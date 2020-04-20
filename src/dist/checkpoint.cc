
#include "src/dist/checkpoint.h"
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include "absl/strings/str_cat.h"
#include "src/agd/errors.h"

using namespace std;

const char* CheckpointFilename = "checkpoint.blob";
const char* TmpCheckpointFilename = "tmp_checkpoint.blob";

// checkpoint file is an indexed binary structure
// [4B uint index size | index (4B uints) | ClusterSets as per index]

bool CheckpointFileExists(const absl::string_view path) {
  std::string checkpoint_file = absl::StrCat(path, CheckpointFilename);
  std::ifstream checkp_stream(checkpoint_file);
  return checkp_stream.good();
}

// path should have a slash on the end
agd::Status WriteCheckpointFile(
    const absl::string_view path,
    const std::unique_ptr<ConcurrentQueue<MarshalledClusterSet>>& queue) {
  std::string tmp_checkpoint_file = absl::StrCat(path, TmpCheckpointFilename);

  std::ofstream checkp_stream(tmp_checkpoint_file);
  if (!checkp_stream.good()) {
    return agd::errors::Internal("Failed to create checkpoint file ",
                                 tmp_checkpoint_file,
                                 ", reason: ", strerror(errno));
  }
  std::vector<uint32_t> index;
  index.reserve(queue->size());
  agd::Buffer buf;

  auto set_iter = queue->begin();
  while (set_iter != queue->end()) {
    buf.AppendBuffer(set_iter->buf.data(), set_iter->buf.size());
    index.push_back(set_iter->TotalSize());
    set_iter++;
  }

  uint32_t sz = index.size();
  checkp_stream.write(reinterpret_cast<char*>(&sz), sizeof(uint32_t));
  checkp_stream.write(reinterpret_cast<char*>(&index[0]),
                      sizeof(uint32_t) * sz);
  checkp_stream.write(buf.data(), buf.size());

  checkp_stream.close();

  // atomically overwrite the old checkpoint
  std::string checkpoint_file = absl::StrCat(path, CheckpointFilename);
  int e = rename(tmp_checkpoint_file.c_str(), checkpoint_file.c_str());
  if (e != 0) {
    return agd::errors::Internal("failed to overwrite checkpoint, returned ", e,
                                 ", reason: ", strerror(errno));
  }

  return agd::Status::OK();
}

agd::Status LoadCheckpointFile(
    const absl::string_view path,
    std::unique_ptr<ConcurrentQueue<MarshalledClusterSet>>& queue) {
  MarshalledClusterSet dummy;
  while (!queue->empty()) {  // empty the queue
    queue->pop(dummy);
  }
  std::string checkpoint_file = absl::StrCat(path, "checkpoint.blob");

  std::ifstream checkp_stream(checkpoint_file, ios::binary);
  if (!checkp_stream.good()) {
    return agd::errors::Internal("unable to open file stream for ",
                                 checkpoint_file, " reason: ", strerror(errno));
  }

  std::vector<uint32_t> index;

  uint32_t index_size;
  checkp_stream.read(reinterpret_cast<char*>(&index_size), sizeof(uint32_t));
  index.resize(index_size);

  checkp_stream.read(reinterpret_cast<char*>(&index[0]),
                     index_size * sizeof(uint32_t));

  for (auto sz : index) {
    MarshalledClusterSet new_set;
    new_set.buf.resize(sz);
    checkp_stream.read(new_set.buf.mutable_data(), sz);
    if (!checkp_stream.good() && !checkp_stream.eof()) {
      return agd::errors::Internal("Unable to read correct number of bytes from stream ",
                                  checkpoint_file, " reason: ", strerror(errno));
    }
    queue->push(std::move(new_set));
  }

  return agd::Status::OK();
}
