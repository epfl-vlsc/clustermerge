
#include "src/dist/checkpoint.h"
#include <fstream>
#include "absl/strings/str_join.h"
#include "src/agd/errors.h"
#include <sys/stat.h>

using namespace std;

bool CheckpointFileExists(const absl::string_view path) {
  struct stat st;
  std::string checkpoint_path = absl::StrJoin(path, "checkpoint/");
  if (stat(checkpoint_path.c_str(), &st) == 0) {
    if (st.st_mode & S_IFDIR == 0) {
      return false;
    }
  } else {
    return false;
  }
  std::string checkpoint_file = absl::StrJoin(path, "checkpoint/checkpoint.blob");
  std::ifstream checkp_stream(checkpoint_file);
  return checkp_stream.good();
}

// path should have a colon on the end
agd::Status WriteCheckpointFile(
    const absl::string_view path,
    const std::unique_ptr<ConcurrentQueue<MarshalledClusterSet>>& queue) {

  std::string checkpoint_tmp_dir = absl::StrJoin(path, "checkpoint_tmp/");
  mkdir(checkpoint_tmp_dir.c_str(), DEFFILEMODE);

  std::string checkpoint_file = absl::StrJoin(checkpoint_tmp_dir, "checkpoint.blob");
  std::string checkpoint_index = absl::StrJoin(checkpoint_tmp_dir, "checkpoint.index");
  std::ofstream checkp_stream(checkpoint_file);
  std::ofstream index_stream(checkpoint_index);
  std::vector<size_t> index;
  index.reserve(queue->size());

  const auto& set_iter = queue->begin();
  while (set_iter != queue->end()) {
    checkp_stream.write(set_iter->buf.data(), set_iter->buf.size());
    index.push_back(set_iter->TotalSize());
  }

  index_stream.write(reinterpret_cast<char*>(&index[0]),
                     sizeof(size_t) * index.size());

  checkp_stream.close();
  index_stream.close();

  // atomically overwrite the old checkpoint
  std::string checkpoint_dir = absl::StrJoin(path, "checkpoint/");
  int e = rename(checkpoint_tmp_dir.c_str(), checkpoint_dir.c_str());
  if (e != 0) {
    return agd::errors::Internal("failed to overwrite checkpoint, returned ", e);
  }

  return agd::Status::OK();
}

agd::Status LoadCheckpoinFile(const absl::string_view path,
                              std::unique_ptr<ConcurrentQueue<MarshalledClusterSet>>& queue) {
  MarshalledClusterSet dummy;
  while (!queue->empty()) {  // empty the queue
    queue->pop(dummy);
  }
  std::string checkpoint_file = absl::StrJoin(path, "checkpoint/checkpoint.blob");
  std::string checkpoint_index = absl::StrJoin(path, "checkpoint/checkpoint.index");
  std::ifstream checkp_stream(checkpoint_file, ios::binary);

  if (!checkp_stream.good()) {
    return agd::errors::Internal("unable to open file stream for ",
                                 checkpoint_file);
  }
  std::ifstream index_stream(checkpoint_index, ios::binary | ios::ate);
  if (!index_stream.good()) {
    return agd::errors::Internal("unable to open file stream for index ",
                                 checkpoint_index);
  }

  std::vector<size_t> index;

  if (index_stream.tellg % 4 != 0) {
    return agd::errors::Internal(
        "checkpoint index file is not multiple of 4 bytes.");
  }

  index.resize(index_stream.tellg() / sizeof(size_t));

  index_stream.read(reinterpret_cast<char*>(&index[0]), index_stream.tellg());

  for (auto sz : index) {
    MarshalledClusterSet new_set;
    new_set.buf.resize(sz);
    checkp_stream.read(new_set.buf.mutable_data(), sz);
    queue->push(new_set);
  }

  return agd::Status::OK();
}