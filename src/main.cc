
#include <iostream>
#include <fstream>
#include "agd/agd_dataset.h"
#include "args.h"
#include "debug.h"
#include "bottom_up_merge.h"
#include "aligner.h"
#include "all_all_executor.h"

using std::cout;
using std::string;
using std::unique_ptr;

agd::Status LoadDatasets(args::PositionalList<std::string>& datasets_opts,
                         std::vector<unique_ptr<agd::AGDDataset>>* datasets) {
    
    agd::Status s = agd::Status::OK();
    for (const auto dataset_opt : args::get(datasets_opts)) {
      unique_ptr<agd::AGDDataset> dataset;
      try {
        s = agd::AGDDataset::Create(dataset_opt, dataset, {"prot"});
      } catch (...) {
        cout << "Error parsing AGD metadata for dataet '" << dataset_opt << "', are you sure it's valid JSON?";
        throw;
      }
      datasets->push_back(std::move(dataset));
    }
    return s;
}

int main(int argc, char** argv) {
  args::ArgumentParser parser("ClusterMerge",
                              "Bottom up protein cluster merge.");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<std::string> json_data_dir(
      parser, "data_dir",
      "Directory containing alignment environment matrices in JSON "
      "(logPAM1.json, all_matrices.json) [data/matrices/json]",
      {"data_dir"});
  args::PositionalList<std::string> datasets_opts(
      parser, "datasets", "AGD Protein datasets to cluster.");

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

  // load alignment envs and initialize
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

  Parameters params;

  // init aligner object
  ProteinAligner aligner(&envs, &params);

  std::vector<unique_ptr<agd::AGDDataset>> datasets;
  agd::Status s;

  if (datasets_opts) {
    // load and parse protein datasets
     s = LoadDatasets(datasets_opts, &datasets);
  } else {
    cout << "No AGD datasets provided. See usage: \n";
    std::cerr << parser;
  }

  if (!s.ok()) {
    cout << s.ToString() << "\n";
    return 0;
  }

  // build initial clustersets
  // one sequence, in one cluster, in one set
  // then, we bottom-up merge them

  cout << "Data loaded, building merger ...\n";
  BottomUpMerge merger(datasets, &aligner);

  AllAllExecutor executor(/*std::thread::hardware_concurrency()*/1, 100, &envs,
                          &params);
  executor.Initialize();

  // pass to Run
  merger.Run(&executor);

  merger.DebugDump();
  // wait and finish call on executor
  // which dumps final results to disk
  executor.FinishAndOutput("outputfoldertest");



  return (0);
}
