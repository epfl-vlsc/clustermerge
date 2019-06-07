
#include "fasta_dataset.h"
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include "src/agd/errors.h"

agd::Status FastaDataset::Create(const std::string file_path,
                                 std::unique_ptr<Dataset>& dataset) {
  auto ds = new FastaDataset();
  auto s = ds->Initialize(file_path);
  if (s.ok()) {
    dataset.reset(ds);
  } else {
    delete ds;
  }

  return s;
}

agd::Status FastaDataset::Initialize(const std::string file_path) {
  std::cout << "Loading FASTA dataset " << file_path << " ...\n";
  auto dot_pos = file_path.find_last_of('.');
  auto slash_or_first_pos = file_path.find_last_of('/');
  if (slash_or_first_pos == std::string::npos) {
    slash_or_first_pos = 0;
  }

  genome_ = file_path.substr(slash_or_first_pos + 1,
                             dot_pos - slash_or_first_pos - 1);
  //std::cout << "genome is " << genome_ << "\n";

  const int fd = open(file_path.c_str(), O_RDONLY);
  struct stat st;
  if (stat(file_path.c_str(), &st) != 0) {
    return agd::errors::Internal("Unable to stat file ", file_path);
  }
  auto size = st.st_size;
  char* mapped = (char*)mmap(0, size, PROT_READ, MAP_SHARED, fd, 0);
  if (mapped == MAP_FAILED) {
    return agd::errors::Internal("Unable to map file ", file_path,
                                 ", returned ", mapped);
  }

  seqs_buf_ = agd::Buffer(size, 1024*1024);
  std::vector<size_t> record_sizes;

  char* cur = mapped;
  while (cur - mapped < size) {
    cur++;  // skip the `>`
    auto line_end = cur;
    while (*line_end != '\n') {
      line_end++;
    }
    metas_.push_back(std::string(cur, line_end - cur));
    total_records_++;
    // std::cout << "pushing back meta: " << metas_.back() << "\n";
    cur = line_end + 1;

    size_t cur_rec_size = 0;
    while (*cur != '>' && cur - mapped < size) {
      line_end = cur;
      while (*line_end != '\n') {
        line_end++;
      }
      cur_rec_size += line_end - cur;
      // std::cout << "adding data: " << std::string(cur, line_end - cur) <<
      // "\n";
      seqs_buf_.AppendBuffer(cur, line_end - cur);
      cur = line_end + 1;  // skip the `\n`
    }
    record_sizes.push_back(cur_rec_size);
  }

  // normalize values
  for (size_t i = 0; i < seqs_buf_.size(); i++) {
    seqs_buf_[i] = seqs_buf_[i] - 'A';
  }

  size_t abs_rec = 0;
  for (auto sz : record_sizes) {
    seqs_.push_back(absl::string_view(&seqs_buf_[abs_rec], sz));
    // std::cout << "pushed back sequence: " << seqs_.back() << "\n";
    abs_rec += sz;
  }

  assert(seqs_.size() == metas_.size());
  total_records_ = seqs_.size();
  munmap(mapped, size);
  close(fd);
  return Status::OK();
}

agd::Status FastaDataset::GetNextRecord(const char** data, size_t* sz,
                                        const char** meta, size_t* meta_sz) {
  if (current_record_ >= total_records_) {
    return agd::errors::OutOfRange("Past last fasta record");
  }
  *data = seqs_[current_record_].data();
  *sz = seqs_[current_record_].size();
  if (meta) {
    if (!meta_sz) {
      return agd::errors::Internal("fasta get next record meta sz was null");
    }
    *meta = metas_[current_record_].data();
    *meta_sz = metas_[current_record_].size();
  }

  current_record_++;

  return agd::Status::OK();
}

agd::Status FastaDataset::GetRecordAt(const size_t index, const char** data,
                                      size_t* sz, const char** meta,
                                      size_t* meta_sz) const {
  if (index >= total_records_) {
    return agd::errors::OutOfRange("index out of range");
  }
  *data = seqs_[index].data();
  *sz = seqs_[index].size();
  if (meta) {
    if (!meta_sz) {
      return agd::errors::Internal("fasta get next record meta sz was null");
    }
    *meta = metas_[index].data();
    *meta_sz = metas_[index].size();
  }

  return agd::Status::OK();
}
