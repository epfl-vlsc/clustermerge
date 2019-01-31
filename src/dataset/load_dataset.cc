
#include "load_dataset.h"
#include <fstream>
#include <iostream>

using std::cout;
using std::string;
using std::unique_ptr;

agd::Status LoadDataset(const string& ext, const string& filename,
                        unique_ptr<Dataset>& dataset) {
  // cout << "extension is " << ext<< "\n";
  auto s = agd::Status::OK();
  if (ext == "json") {  // its agd
    try {
      s = AGDProteinDataset::Create(filename, dataset);
    } catch (...) {
      cout << "Error parsing AGD metadata for dataet '" << filename
           << "', are you sure it's valid JSON?";
      throw;
    }
  } else if (ext == "fasta" || ext == "fa") {  // its fasta
    s = FastaDataset::Create(filename, dataset);
  } else {
    s = agd::errors::InvalidArgument("unrecognized dataset file type: ",
                                     filename);
  }

  return s;
}

agd::Status LoadDatasetsPositional(
    args::PositionalList<std::string>& datasets_opts,
    std::vector<unique_ptr<Dataset>>* datasets) {
  agd::Status s = agd::Status::OK();
  for (const auto dataset_opt : args::get(datasets_opts)) {
    unique_ptr<Dataset> dataset;
    auto dot_pos = dataset_opt.find_last_of('.');
    auto extension = dataset_opt.substr(dot_pos + 1);

    s = LoadDataset(extension, dataset_opt, dataset);
    if (!s.ok()) {
      return s;
    }

    /*cout << "data is " << dataset->Name() << " with size " << dataset->Size()
         << "\n";*/
    datasets->push_back(std::move(dataset));
  }
  return s;
}

agd::Status LoadDatasetsJSON(const string& dataset_file,
                             std::vector<unique_ptr<Dataset>>* datasets) {
  std::ifstream dataset_stream(dataset_file);

  if (!dataset_stream.good()) {
    return agd::errors::NotFound("No such file: ", dataset_file);
  }

  json dataset_json_obj;
  dataset_stream >> dataset_json_obj;

  agd::Status s = agd::Status::OK();
  for (const auto dataset_opt : dataset_json_obj) {
    auto dataset_opt_str = dataset_opt.get<std::string>();
    auto dot_pos = dataset_opt_str.find_last_of('.');
    auto extension = dataset_opt_str.substr(dot_pos + 1);
    unique_ptr<Dataset> dataset;

    s = LoadDataset(extension, dataset_opt, dataset);
    if (!s.ok()) {
      return s;
    }

    /*cout << "data is " << dataset->Name() << " with size " << dataset->Size()
         << "\n";*/
    datasets->push_back(std::move(dataset));
  }
  return s;
}