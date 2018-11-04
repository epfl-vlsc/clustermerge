
// Entry point for the master server for the distributed
// clustering implementation

#include <iostream>
#include <thread>
#include <fstream>
#include "args.h"
#include "src/dataset/load_dataset.h"
#include "src/common/alignment_environment.h"
#include "src/common/params.h"
#include "src/dist/proto/requests.pb.h"
#include "zmq.hpp"

using namespace std;

constexpr char cluster_format[] =
    "{\n"
    " \"controller\": \"<ip/addr>\",\n"
    " \"push_port\": <port num>,\n"
    " \"pull_port\": <port num>,\n"
    "}\n";

/*
Server cluster format example
{
  "controller": "<ip/addr>",
  "push_port": <port num>,
  "pull_port": <port num>
}
*/

int main(int argc, char* argv[]) {
  args::ArgumentParser parser("ClusterMerge",
                              "Bottom up protein cluster merge.");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<std::string> json_data_dir(
      parser, "data_dir",
      "Directory containing alignment environment matrices in JSON "
      "(logPAM1.json, all_matrices.json) [data/matrices/json]",
      {'d', "data_dir"});
  args::ValueFlag<unsigned int> threads_arg(
      parser, "threads",
      absl::StrCat("Number of threads to use for all-all [",
                   std::thread::hardware_concurrency(), "]"),
      {'t', "threads"});
  args::ValueFlag<unsigned int> cluster_threads_arg(
      parser, "cluster_threads",
      absl::StrCat("Number of threads to use for clustering [",
                   std::thread::hardware_concurrency(), "]"),
      {'c', "cluster-threads"});
  args::ValueFlag<unsigned int> merge_threads_arg(
      parser, "merge_threads",
      absl::StrCat("Number of threads to use for merging [",
                   std::thread::hardware_concurrency(), "]"),
      {'m', "merge-threads"});
  args::ValueFlag<std::string> input_file_list(
      parser, "file_list", "JSON containing list of input AGD datasets.",
      {'i', "input_list"});
  args::ValueFlag<std::string> output_dir(
      parser, "output_dir",
      "Output directory. Will be overwritten if exists."
      "[./output_matches]",
      {'o', "output_dir"});
  args::ValueFlag<std::string> server_config_file(
      parser, "server_config",
      absl::StrCat(
          "JSON file containing the cluster server configuration, example: \n",
          cluster_format),
      {'s', "server_config"});
  args::PositionalList<std::string> datasets_opts(
      parser, "datasets",
      "AGD Protein datasets to cluster. If present, will override `input_list` "
      "argument.");
  args::Flag controller(
      parser, "controller",
      "Designate this process as the cluster controller.",
      {'C', "controller"});

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  } catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // parse the server cluster config file to see if we are a worker or the controller

  bool is_controller = false;
  if (controller) {
    is_controller = args::get(controller);
  }

  json server_config_json;
  string controller_ip;
  int push_port;
  int pull_port;
  if (server_config_file) {
    string server_config_path = args::get(server_config_file);
    std::ifstream server_config_stream(server_config_path);

    if (!server_config_stream.good()) {
      std::cerr << "File " << server_config_path << " not found.\n";
      return 1;
    }
  
    server_config_stream >> server_config_json;

  } else {
    cout << "The server config file (--server_config) is required.\n";
    return 1;
  }

  auto controller_it = server_config_json.find("controller");
  if (controller_it == server_config_json.end() && !is_controller) {
    cout << "This process is not controller and server config lacks controller addr\n";
    return 1;
  } else {
    controller_ip = *controller_it;
  }

  auto push_port_it = server_config_json.find("push_port");
  if (push_port_it == server_config_json.end()) {
    push_port = 5555; // default
  } else {
    push_port = *push_port_it;
  }
  
  auto pull_port_it = server_config_json.find("push_port");
  if (pull_port_it == server_config_json.end()) {
    pull_port = 5555; // default
  } else {
    pull_port = *push_port_it;
  }

  unsigned int threads = std::thread::hardware_concurrency();
  if (threads_arg) {
    threads = args::get(threads_arg);
    // do not allow more than hardware threads
    if (threads > std::thread::hardware_concurrency()) {
      threads = std::thread::hardware_concurrency();
    }
  }
  cout << "Using " << threads << " hardware threads for alignment.\n";

  unsigned int cluster_threads = std::thread::hardware_concurrency();
  if (cluster_threads_arg) {
    cluster_threads = args::get(cluster_threads_arg);
    // do not allow more than hardware threads
    if (cluster_threads > std::thread::hardware_concurrency()) {
      cluster_threads = std::thread::hardware_concurrency();
    }
  }
  cout << "Using " << cluster_threads << " hardware threads for clustering.\n";

  unsigned int merge_threads = std::thread::hardware_concurrency();
  if (merge_threads_arg) {
    merge_threads = args::get(merge_threads_arg);
    // do not allow more than hardware threads
    if (merge_threads > std::thread::hardware_concurrency()) {
      merge_threads = std::thread::hardware_concurrency();
    }
  }
  cout << "Using " << merge_threads << " hardware threads for merging.\n";

  // get output dir to use, only needed if controller
  string dir("output_matches");
  if (output_dir) {
    dir = args::get(output_dir);
  }
  cout << "Using " << dir << " for output.\n";

  // load alignment envs and initialize (this is for SWPS3)
  string json_dir_path = "data/matrices/json/";
  if (json_data_dir) {
    json_dir_path = args::get(json_data_dir);
    // OSX  STILL  doesn't have experimental/filesystem so this isn't portable
    if (json_dir_path.back() != '/') {
      json_dir_path += '/';
    }
  }
  string logpam_json_file = json_dir_path + "logPAM1.json";
  string all_matrices_json_file = json_dir_path + "all_matrices.json";

  std::ifstream logpam_stream(logpam_json_file);
  std::ifstream allmat_stream(all_matrices_json_file);

  if (!logpam_stream.good()) {
    std::cerr << "File " << logpam_json_file << " not found.\n";
    return 1;
  }

  if (!allmat_stream.good()) {
    std::cerr << "File " << all_matrices_json_file << " not found.\n";
    return 1;
  }

  json logpam_json;
  logpam_stream >> logpam_json;

  json all_matrices_json;
  allmat_stream >> all_matrices_json;

  AlignmentEnvironments envs;

  // initializing envs is expensive, so don't copy this
  cout << "Initializing environments from " << json_dir_path << "\n";
  envs.InitFromJSON(logpam_json, all_matrices_json);
  cout << "Done.\n";

  // done init envs

  Parameters params;  // using default params for now


  if (is_controller) {
    // launch controller(push_port, pull_port)
  } else {
    // launch worker(args, controller_ip, push_port, pull_port)  
  }




  cmproto::Cluster c;
  c.add_indexes(1);

  cmproto::MergeBatch batch;
  auto* set = batch.add_sets();
  auto* cluster = set->add_clusters();
  *cluster = c;

  cmproto::MergeRequest r;
  /*auto* b = r.mutable_batch();
   *b = batch;*/

  if (r.has_batch()) {
    cout << "its a batch\n";
  } else if (r.has_partial()) {
    cout << "has partial\n";
  } else {
    cout << "has none\n";
  }

  // parse servers json, determine role
  // load datasets

  // if controller
  // connect to zmq queues
  // push initial batch mergers to ingress queue
  // loop (possibly with threads), pulling from egress queue
  //   construct new batches
  //   if received cluster set is large (some threshold?), construct partial
  //   merge requests

  // else if worker server
  return 0;
}