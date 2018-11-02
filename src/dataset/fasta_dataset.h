
#pragma once

#include <vector>
#include "absl/strings/string_view.h"
#include "dataset.h"
#include "src/agd/buffer.h"

// TODO agd and fasta dataset should inherit from an abstract dataset

class FastaDataset : public Dataset {
 public:
  static agd::Status Create(const std::string file_path,
                            std::unique_ptr<Dataset>& dataset);

  agd::Status GetNextRecord(const char** data, size_t* sz,
                            const char** meta = nullptr,
                            size_t* meta_sz = nullptr) override;

  agd::Status GetRecordAt(const size_t index, const char** data, size_t* sz,
                          const char** meta = nullptr,
                          size_t* meta_sz = nullptr) const override;

  void Reset() override { current_record_ = 0; }

 private:
  FastaDataset() = default;
  agd::Status Initialize(const std::string file_path);
  // contiguous buffer to hold all sequences
  std::vector<absl::string_view> seqs_;
  std::vector<std::string> metas_;
  agd::Buffer seqs_buf_;
  std::string genome_;  // the filename basically
};