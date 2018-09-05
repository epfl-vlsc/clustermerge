
#pragma once
#include "agd_record_reader.h"
#include "json.hpp"

namespace agd {
using json = nlohmann::json;

class AGDDataset {
 public:
  // create dataset
  // load only specific columns if given
  static Status Create(const std::string& agd_json_path,
                       std::unique_ptr<AGDDataset>& dataset,
                       std::vector<std::string> columns = {});

  // allows iteration over a complete agd column
  class ColumnIterator {
    friend class AGDDataset;

   public:
    ColumnIterator() = default;
    Status GetNextRecord(const char** data, size_t* size);
    void Reset();

   private:
    ColumnIterator(std::vector<AGDRecordReader>* readers)
        : column_readers_(readers) {}

    std::vector<AGDRecordReader>* column_readers_;
    uint32_t current_reader_ = 0;
  };

  Status Column(const std::string& column, ColumnIterator* iter) {
    if (column_map_.find(column) != column_map_.end()) {
      *iter = ColumnIterator(&column_map_[column]);
      return Status::OK();
    } else {
      return NotFound("column ", column, " not found in column map.");
    }
  }

  uint32_t Size() const { return total_records_; }

  const std::string& Name() const { return name_; }

 private:
  AGDDataset() = default;
  Status Initialize(const std::string& agd_json_path,
                    const std::vector<std::string>& columns);
  std::string name_;

  json agd_metadata_;

  std::vector<Buffer> chunks_;
  // map column name to its set of readers
  std::unordered_map<std::string, std::vector<AGDRecordReader>> column_map_;
  std::vector<size_t> chunk_sizes_;
  size_t total_records_;
};

}  // namespace agd