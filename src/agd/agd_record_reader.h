#pragma once

#include <vector>
#include "errors.h"
#include "format.h"

namespace agd {

using namespace errors;

/*
 * A class that provides a "view" over the data in the resource container.
 * Does not take ownership of the underlying data
 *
 * This is the class to inherit from if you want another interface to
 * AGD chunk records.
 */
class AGDRecordReader {
 public:
  AGDRecordReader(const char* resource, size_t num_records);

  void Reset();
  int NumRecords() { return num_records_; }

  Status GetNextRecord(const char** data, size_t* size);
  Status PeekNextRecord(const char** data, size_t* size);

  Status GetRecordAt(size_t index, const char** data, size_t* size);

  size_t GetCurrentIndex() { return cur_record_; }

 private:
  const format::RelativeIndex* index_;
  const char *data_, *cur_data_;
  size_t cur_record_ = 0;
  std::vector<size_t> absolute_index_;

  void InitializeIndex();

 protected:
  size_t num_records_;
};

}  // namespace agd
