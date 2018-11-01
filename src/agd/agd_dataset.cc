
#include "agd_dataset.h"
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include "absl/strings/str_cat.h"
#include "parser.h"

namespace agd {

using std::ifstream;
using std::string;
using std::unique_ptr;
using std::vector;

Status AGDDataset::Create(const string& agd_json_path,
                          unique_ptr<AGDDataset>& dataset,
                          vector<string> columns) {

  std::cout << "Loading AGD dataset '" << agd_json_path << "' ...\n";
  dataset.reset(new AGDDataset());
  Status s = dataset->Initialize(agd_json_path, columns);
  if (!s.ok()) {
    dataset.reset();
  }

  return s;
}

Status AGDDataset::Initialize(const string& agd_json_path,
                              const vector<string>& columns) {
  ifstream i(agd_json_path);
  json agd_metadata;
  i >> agd_metadata;
  name_ = agd_metadata["name"];
  agd_metadata_ = agd_metadata;

  vector<string> to_load;
  const auto& cols = agd_metadata["columns"];
  if (!columns.empty()) {
    for (const auto& c : columns) {
      bool found = false;
      for (const auto& json_col : cols) {
        if (json_col == c) {
          found = true;
          break;
        }
      }
      if (!found) {
        return NotFound("column ", c, " was not found in dataset.");
      }
    }
    to_load = columns;
  } else {
    for (const auto& c : cols) {
      to_load.push_back(c);
    }
  }

  string file_path_base =
      agd_json_path.substr(0, agd_json_path.find_last_of('/') + 1);

  RecordParser parser;

  for (const auto& c : to_load) {

    for (const auto& chunk : agd_metadata["records"]) {
      string chunk_name = chunk["path"];
      string path = absl::StrCat(file_path_base, chunk_name, ".", c);

      const int fd = open(path.c_str(), O_RDONLY);
      struct stat st;
      if (stat(path.c_str(), &st) != 0) {
        return Internal("Unable to stat file ", path);
      }
      auto size = st.st_size;
      char* mapped = (char*)mmap(0, size, PROT_READ, MAP_SHARED, fd, 0);
      if (mapped == MAP_FAILED) {
        return Internal("Unable to map file ", path, ", returned ", mapped);
      }

      Buffer chunk_buf;
      uint64_t first_ordinal;
      uint32_t num_records;
      string record_id;
      parser.ParseNew(mapped, size, false, &chunk_buf, &first_ordinal,
                      &num_records, record_id);

      column_map_[c].push_back(AGDRecordReader(chunk_buf.data(), num_records));
      chunks_.push_back(std::move(chunk_buf));
      chunk_sizes_.push_back(num_records);
      total_records_ += num_records;

      munmap(mapped, size);
      close(fd);
    }
  }

  return Status::OK();
}


Status AGDDataset::ColumnIterator::GetNextRecord(const char** data, size_t* size) {

  Status s = (*column_readers_)[current_reader_].GetNextRecord(data, size);
  if (!s.ok()) {
    if (current_reader_ < column_readers_->size() - 1) {
      current_reader_++;
      return (*column_readers_)[current_reader_].GetNextRecord(data, size);
    } else {
      return OutOfRange("Last record in column");
    }
  }
  return s;
}

Status AGDDataset::ColumnIterator::GetNextAt(size_t index, const char** data, size_t* size) {

  if (index >= total_records_) {
    return OutOfRange("index is greater than total records in dataset");
  }

  size_t chunk_index = index / (*column_readers_)[0].NumRecords();
  size_t record_index = index % (*column_readers_)[0].NumRecords();
  Status s = (*column_readers_)[chunk_index].GetRecordAt(record_index, data, size);
  if (!s.ok()) {
    return OutOfRange("Last record in column");
  }
  return s;
}

}  // namespace agd
