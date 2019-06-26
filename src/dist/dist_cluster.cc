//
// Entry point for the master server for the distributed
// clustering implementation
//
// Stuart Byma, EPFL
//

#include <limits.h>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>
#include "args.h"
#include "controller.h"
#include "src/common/alignment_environment.h"
#include "src/common/params.h"
#include "src/dataset/load_dataset.h"
#include "src/common/cluster_set.h"
#include "worker.h"

using namespace std;

constexpr char cluster_format[] =
    "{\n"
    " \"controller\": \"<ip/addr>\",\n"
    " \"request_queue_port\": <port num>,\n"
    " \"response_queue_port\": <port num>,\n"
    "}\n";

/*
Server cluster format example
{
  "controller": "<ip/addr>",
  "request_queue_port": <port num>,
  "pull_port": <port num>
}
*/

constexpr char cluster_config_default[] = "data/default_cluster.json";

#define DEFAULT_QUEUE_DEPTH 5
#define DEFAULT_MAX_SET_SIZE \
  500  // max set size parameter for large partial merges
#define DEFAULT_RESPONSE_QUEUE_PORT 5556
#define DEFAULT_REQUEST_QUEUE_PORT 5555
#define DEFAULT_INCOMPLETE_REQUEST_QUEUE_PORT 5554
#define DEFAULT_LARGE_PARTIAL_MERGE_PORT 5557

// captures kill signal and notifies Worker
volatile int signal_num = 0;

void signal_notifier(int sig) {
  signal_num = sig;
  cout << "[signal_num] value changed.\n";
}

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
  args::ValueFlag<unsigned int> queue_depth_arg(
      parser, "queue_depth", "Depth of the local work and response queue.",
      {'q', "queue_depth"});
  args::ValueFlag<unsigned int> batch_size_arg(
      parser, "batch size",
      "How many small clusters should be batched together.",
      {'b', "batch_size"});
  args::ValueFlag<unsigned int> dataset_limit_arg(
      parser, "dataset limit", "Cluster only this many sequences",
      {'D', "dataset_limit"});
  args::ValueFlag<unsigned int> dup_removal_threshold_arg(
      parser, "duplicate removal threshold",
      "How big a set of clusters should be before duplicates are filtered out "
      "[MAX_INT]",
      {'r', "dup_removal_thresh"});
  args::ValueFlag<std::string> input_file_list(
      parser, "file_list", "JSON containing list of input AGD datasets.",
      {'i', "input_list"});
  args::ValueFlag<std::string> aligner_params_arg(
      parser, "aligner parameters",
      "JSON containing alignment and clustering parameters.",
      {'a', "aligner_params"});
  args::ValueFlag<std::string> output_dir(
      parser, "output dir",
      "Output directory. Will be overwritten if exists."
      "[./dist_output_dir]",
      {'o', "output_dir"});
  args::ValueFlag<std::string> server_config_file(
      parser, "server_config",
      absl::StrCat(
          "JSON file containing the cluster server configuration, example: \n",
          cluster_format),
      {'s', "server_config"});
  args::Flag exclude_allall(
      parser, "exclude_allall",
      "Don't perform intra-cluster all-all alignment, just do the clustering.",
      {'x', "exclude_allall"});
  args::PositionalList<std::string> datasets_opts(
      parser, "datasets",
      "AGD Protein datasets to cluster. If present, will override `input_list` "
      "argument.");
  args::Flag controller(parser, "controller",
                        "Designate this process as the cluster controller.",
                        {'C', "controller"});
  args::ValueFlag<unsigned int> checkpoint_interval_arg(
      parser, "checkpoint interval",
      "how long to wait in seconds between checkpointing current computation. "
      "Recommended value 900s"
      " [0 (off)]",
      {'c', "checkpoint_interval"});
  args::ValueFlag<std::string> file_name(
      parser, "file",
      "Adds clustering data from an already clustered json file and merges "
      "with clusters from current dataset",
      {'f', "file"});

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

  long int checkpoint_interval = 0;
  if (checkpoint_interval_arg) {
    checkpoint_interval = args::get(checkpoint_interval_arg);
    cout << "checkpoint interval: " << checkpoint_interval << "\n";
  }

  // parse the server cluster config file to see if we are a worker or the
  // controller

  bool is_controller = false;
  if (controller) {
    is_controller = args::get(controller);
  }

  uint32_t batch_size = 100;
  if (batch_size_arg) {
    batch_size = args::get(batch_size_arg);
  }

  uint32_t dup_removal_threshold = UINT_MAX;
  if (dup_removal_threshold_arg) {
    dup_removal_threshold = args::get(dup_removal_threshold_arg);
  }

  json server_config_json;
  string controller_ip;
  size_t max_set_size;
  int request_queue_port;
  int response_queue_port;
  int incomplete_request_queue_port;
  int set_request_port;
  if (server_config_file) {
    string server_config_path = args::get(server_config_file);
    std::ifstream server_config_stream(server_config_path);

    if (!server_config_stream.good()) {
      std::cerr << "File " << server_config_path << " not found.\n";
      return 1;
    }

    server_config_stream >> server_config_json;

  } else {
    std::ifstream server_config_stream(cluster_config_default);

    if (!server_config_stream.good()) {
      std::cerr << "File " << cluster_config_default << " not found.\n";
      return 1;
    }

    server_config_stream >> server_config_json;
  }

  auto controller_it = server_config_json.find("controller");
  if (controller_it == server_config_json.end() && !is_controller) {
    cout << "This process is not controller and server config lacks controller "
            "addr\n";
    return 1;
  } else {
    controller_ip = *controller_it;
  }
  cout << "controller ip: " << controller_ip << "\n";
  if (is_controller) {
    cout << "The process is controller\n";
  }

  auto request_queue_port_it = server_config_json.find("request_queue_port");
  if (request_queue_port_it == server_config_json.end()) {
    request_queue_port = DEFAULT_REQUEST_QUEUE_PORT;  // default
  } else {
    request_queue_port = *request_queue_port_it;
  }

  auto response_queue_port_it = server_config_json.find("response_queue_port");
  if (response_queue_port_it == server_config_json.end()) {
    response_queue_port = DEFAULT_RESPONSE_QUEUE_PORT;  // default
  } else {
    response_queue_port = *response_queue_port_it;
  }

  auto set_request_port_it = server_config_json.find("set_request_port");
  if (set_request_port_it == server_config_json.end()) {
    set_request_port = DEFAULT_LARGE_PARTIAL_MERGE_PORT;  // default
  } else {
    set_request_port = *set_request_port_it;
  }

  auto incomplete_request_queue_port_it =
      server_config_json.find("incomplete_request_queue_port");
  if (incomplete_request_queue_port_it == server_config_json.end()) {
    incomplete_request_queue_port =
        DEFAULT_INCOMPLETE_REQUEST_QUEUE_PORT;  // default
  } else {
    incomplete_request_queue_port = *incomplete_request_queue_port_it;
  }

  auto max_set_size_it = server_config_json.find("max_set_size");
  if (max_set_size_it == server_config_json.end()) {
    max_set_size = DEFAULT_MAX_SET_SIZE;  // default
  } else {
    max_set_size = *max_set_size_it;
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

  unsigned int queue_depth = DEFAULT_QUEUE_DEPTH;
  if (queue_depth_arg) {
    queue_depth = args::get(queue_depth_arg);
  }

  // get output dir to use, only needed if controller
  string dir("dist_output_dir");
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

  std::vector<unique_ptr<Dataset>> datasets;
  agd::Status s;
  std::string input_file_name_temp = args::get(input_file_list);
  
  if (datasets_opts) {
    // load and parse protein datasets
    // cluster merge sequences are simply string_views over this data
    // so these structures must live until computations complete
    if (input_file_list) {
      cout << "WARNING: ignoring input file list and using positionals!\n";
    }
    s = LoadDatasetsPositional(datasets_opts, &datasets);
  } else if (input_file_list) {
    s = LoadDatasetsJSON(input_file_name_temp, &datasets);
  } else {
    cout << "No datasets provided. See usage: \n";
    std::cerr << parser;
    exit(0);
  }

  if (!s.ok()) {
    cout << s.ToString() << "\n";
    return 0;
  }

  Parameters aligner_params;
  if (aligner_params_arg) {
    json aligner_params_json;
    string aligner_params_file = args::get(aligner_params_arg);
    std::ifstream aligner_params_stream(aligner_params_file);

    if (!aligner_params_stream.good()) {
      std::cerr << "File " << aligner_params_file << " not found.\n";
      return 1;
    }

    aligner_params_stream >> aligner_params_json;

    auto min_score_it = aligner_params_json.find("min_score");
    if (min_score_it != aligner_params_json.end()) {
      aligner_params.min_score = *min_score_it;
    }

    auto max_aa_uncovered_it = aligner_params_json.find("max_aa_uncovered");
    if (max_aa_uncovered_it != aligner_params_json.end()) {
      aligner_params.max_n_aa_not_covered = *max_aa_uncovered_it;
    }

    auto min_full_merge_score_it =
        aligner_params_json.find("min_full_merge_score");
    if (min_full_merge_score_it != aligner_params_json.end()) {
      aligner_params.min_full_merge_score = *min_full_merge_score_it;
    }

    auto blosum_it = aligner_params_json.find("blosum");
    if (blosum_it != aligner_params_json.end()) {
      aligner_params.use_blosum = *blosum_it;
    }
  }  // if not present, aligner params defaults used

  if (is_controller) {
    
    //stuff for handling preexisting merge result
    //parse and build a Datasets object
    std::string input_file_name_temp = args::get(input_file_list);
    std::vector<std::string> dataset_file_names;
    std::vector<unique_ptr<Dataset>> datasets_old;
    json dataset_json_obj;
    ClusterSet old_set;
      
    if (file_name) {
      agd::Status s_old;

      std::ifstream dataset_stream(args::get(file_name));

      if (!dataset_stream.good()) {
        s_old = agd::errors::NotFound("No such file: ", args::get(file_name));
      }

      if (!s_old.ok()) {
        cout << s_old.ToString() << "\n";
        return 0;
      }

      dataset_stream >> dataset_json_obj;

      for (const auto& old_dataset : dataset_json_obj["datasets"]) {
        s_old = LoadDatasetsJSON(old_dataset, &datasets_old);
        dataset_file_names.push_back(old_dataset);

        if (!s_old.ok()) {
          cout << s_old.ToString() << "\n";
          return 0;
        }
      }
     }

    dataset_file_names.push_back(input_file_name_temp);

    // launch controller(push_port, pull_port)
    // TODO put checkpoint dir in params file and ensure exists
    string checkpoint_dir("./");
    
    Controller::Params params;
    params.batch_size = batch_size;
    params.controller_ip = controller_ip;
    params.data_dir_path = json_dir_path;
    params.num_threads = threads;
    params.queue_depth = queue_depth;
    params.request_queue_port = request_queue_port;
    params.response_queue_port = response_queue_port;
    params.incomplete_request_queue_port = incomplete_request_queue_port;
    params.set_request_port = set_request_port;
    params.dup_removal_thresh = dup_removal_threshold;
    params.max_set_size = max_set_size;
    params.exclude_allall = exclude_allall;
    params.checkpoint_interval = checkpoint_interval;
    params.checkpoint_dir = absl::string_view(checkpoint_dir);
    if (dataset_limit_arg) {
      params.dataset_limit = args::get(dataset_limit_arg);
    } else {
      params.dataset_limit = -1;
    }

    Status stat;
    if(file_name) {
      Controller controller(dataset_json_obj, datasets_old);
      stat = controller.Run(params, aligner_params, datasets, dataset_file_names);
    } else {
      Controller controller;
      stat = controller.Run(params, aligner_params, datasets, dataset_file_names);
    }
    if (!stat.ok()) {
      cout << "Error: " << stat.error_message() << "\n";
      return -1;
    }
  } else {
    // load datasets, launch worker
    // launch worker(args, controller_ip, push_port, pull_port)
    Worker worker;
    Worker::Params params;
    params.controller_ip = controller_ip;
    params.data_dir_path = json_dir_path;
    params.num_threads = threads;
    params.queue_depth = queue_depth;
    params.request_queue_port = request_queue_port;
    params.response_queue_port = response_queue_port;
    params.incomplete_request_queue_port = incomplete_request_queue_port;
    params.set_request_port = set_request_port;
    signal(SIGUSR1, signal_notifier);
    // signal(SIGINT, signal_notifier);
    Status stat =
        worker.Run(params, aligner_params, datasets, (int*)&signal_num);
    if (!stat.ok()) {
      cout << "Error: " << stat.error_message() << "\n";
      return -1;
    }
  }

  return 0;
}
