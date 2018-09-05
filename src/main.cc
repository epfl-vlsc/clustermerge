
#include <iostream>
#include "agd/agd_dataset.h"
#include "args.h"
#include "debug.h"
#include "bottom_up_merge.h"

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
  args::PositionalList<std::string> datasets_opts(parser, "datasets",
                                         "AGD Protein datasets to cluster.");

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
  BottomUpMerge merger(datasets);
  merger.DebugDump();


  /*agd::AGDDataset::ColumnIterator iter;
  s = dataset->Column("prot", &iter);*/

  /*if (!s.ok()) {
    cout << s.ToString() << "\n";
    return 0;
  }

  const char* data;
  size_t size;
  s = iter.GetNextRecord(&data, &size);
  while (s.ok()) {
    cout << PrintNormalizedProtein(data, size) << "\n\n";
    s = iter.GetNextRecord(&data, &size);
  }*/

  return (0);
}
