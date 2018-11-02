
#include "agd_protein_dataset.h"
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include "absl/strings/str_cat.h"
#include "src/agd/parser.h"

using std::ifstream;
using std::string;
using std::unique_ptr;
using std::vector;

Status AGDProteinDataset::Create(const string& file_path,
                                 unique_ptr<Dataset>& dataset) {
  auto ds = new AGDProteinDataset();
  Status s = ds->Initialize(file_path);
  if (s.ok()) {
    dataset.reset(ds);
  } else {
    delete ds;
  }

  return s;
}

Status AGDProteinDataset::Initialize(const string& agd_json_path) {
  auto s =
      agd::AGDDataset::Create(agd_json_path, dataset_, {"prot", "metadata"});
  if (!s.ok()) {
    return s;
  }

  total_records_ = dataset_->Size();

  s = dataset_->Column("prot", &prot_iterator_);
  if (!s.ok()) {
    return s;
  }
  s = dataset_->Column("metadata", &meta_iterator_);
  if (!s.ok()) {
    return s;
  }

  return Status::OK();
}

Status AGDProteinDataset::GetNextRecord(const char** data, size_t* sz,
                                        const char** meta, size_t* meta_sz) {
  auto s = prot_iterator_.GetNextRecord(data, sz);
  if (!s.ok()) {
    return s;
  }
  if (meta) {
    s = meta_iterator_.GetNextRecord(data, sz);
    if (!s.ok()) {
      return s;
    }
  }
  return Status::OK();
}

Status AGDProteinDataset::GetRecordAt(const size_t index, const char** data,
                                      size_t* sz, const char** meta,
                                      size_t* meta_sz) const {
  return agd::errors::Internal("unimplemented get record at in protein db");
}
