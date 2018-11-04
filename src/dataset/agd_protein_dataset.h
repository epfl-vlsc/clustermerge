
#pragma once
#include "src/agd/agd_dataset.h"
#include "dataset.h"
#include "json.hpp"

using json = nlohmann::json;
using agd::AGDRecordReader;
using agd::Status;

class AGDProteinDataset : public Dataset {
 public:
  // create dataset
  // load only specific columns if given
  static Status Create(const std::string& file_path,
                       std::unique_ptr<Dataset>& dataset);

  Status GetNextRecord(const char** data, size_t* sz,
                       const char** meta = nullptr,
                       size_t* meta_sz = nullptr) override;

  Status GetRecordAt(const size_t index, const char** data, size_t* sz,
                     const char** meta = nullptr,
                     size_t* meta_sz = nullptr) const override;

  void Reset() override {
    prot_iterator_.Reset();
    meta_iterator_.Reset();
  }

 private:
  Status Initialize(const std::string& agd_json_path);

  std::unique_ptr<agd::AGDDataset> dataset_;
  agd::AGDDataset::ColumnIterator prot_iterator_;
  agd::AGDDataset::ColumnIterator meta_iterator_;
};
