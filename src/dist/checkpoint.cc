
#include "src/dist/checkpoint.h"
#include <ifstream>
#include "absl/strings/str_join.h"

bool CheckpointFileExists(const std::string& path) {
  std::string checkpoint_file = absl::StrJoin(path, "checkpoint.blob");
  std::ifstream checkp_stream(checkpoint_file);
  return checkp_stream.good();
}

agd::Status WriteCheckpointFile(const std::string& path,
    const ConcurrentQueue<MarshalledClusterSet>& queue) {
  
  std::string checkpoint_file = absl::StrJoin(path, "checkpoint.blob");
  std::string checkpoint_index = absl::StrJoin(path, "checkpoint.index");
  std::ofstream checkp_stream(checkpoint_file);
  std::ofstream index_stream(checkpoint_index);
  std::vector<size_t> index;
  index.reserve(queue.size());

  const auto& set_iter = queue.begin();
  while (set_iter != queue.end()) {
    checkp_stream.write(set_iter->buf.data(), set_iter->buf.size());
    index.push_back(set_iter->Size());
  }

  checkpoint_index.write(&index[0], sizeof(size_t)*index.size());
}

agd::Status LoadCheckpoinFile(
    const std::string& path,
    const ConcurrentQueue<MarshalledClusterSet>& queue) {

  queue.clear();
  std::string checkpoint_file = absl::StrJoin(path, "checkpoint.blob");
  std::string checkpoint_index = absl::StrJoin(path, "checkpoint.index");
  std::ifstream checkp_stream(checkpoint_file, ios::binary);
  std::ifstream index_stream(checkpoint_index, ios::binary | ios::ate);
  std::vector<size_t> index;

  if (index_stream.tellg % 4 != 0) {
    return errors::Internal("checkpoint index file is not multiple of 4 bytes.");
  }

  index.resize(index_stream.tellg()/sizeof(size_t));

  index_stream.read(&index[0], index_stream.tellg());

  for (auto sz : index) {
    MarshalledClusterSet new_set;
    new_set.buf.resize(sz);
    checkp_stream.read(new_set.buf.mutable_data(), sz);
    queue.push(new_set);
  }

  return Status::OK();
}