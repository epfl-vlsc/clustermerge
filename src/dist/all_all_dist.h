
#pragma once

#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <vector>
#include "absl/container/node_hash_map.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "src/common/aligner.h"
#include "src/common/all_all_base.h"
#include "src/common/candidate_map.h"
#include "src/common/cluster_set.h"
#include "src/common/concurrent_queue.h"
#include "src/common/params.h"
#include "src/common/sequence.h"
#include "src/comms/requests.h"
#include "zmq.hpp"

// class to manage farming out alignments to remote workers
// we use the existing zmq queues to send requests to existing workers
class AllAllDist : public AllAllBase {
 public:
  AllAllDist(ConcurrentQueue<MarshalledRequest>* req_queue,
             const std::vector<Sequence>& seqs, const std::string& output_dir)
      : request_queue_(req_queue), sequences_(seqs), output_dir_(output_dir) {
  
    struct stat info;
    if (stat(output_dir_.c_str(), &info) != 0) {
      // doesnt exist, create
      std::cout << "creating dir " << output_dir << "\n";
      int e = mkdir(output_dir_.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (e != 0) {
        std::cout << "could not create output dir " << output_dir
             << ", exiting ...\n";
        exit(0);
      }
    } else if (!(info.st_mode & S_IFDIR)) {
      // exists but not dir
      std::cout << "output dir exists but is not dir, exiting ...\n";
      exit(0);
    } else {
      // dir exists, nuke
      // im too lazy to do this the proper way
      std::string cmd = absl::StrCat("rm -rf ", output_dir, "/*");
      std::cout << "dir " << output_dir << " exists, nuking ...\n";
      int nuke_result = system(cmd.c_str());
      if (nuke_result != 0) {
        std::cout << "Could not nuke dir " << output_dir << "\n";
        exit(0);
      }
    }
  }

  // final set alignment scheduling is processed by a single thread
  void EnqueueAlignment(const WorkItem& item) override;

  // may be called concurrently
  void ProcessResult(const char* result_buf);

  // send remaining buffered requests
  // wait until all finished
  // output and close all files
  void Finish();

 private:
  ConcurrentQueue<MarshalledRequest>* request_queue_;

  // to buffer alignments into groups before submitting
  MarshalledRequest req_;
  size_t cur_num_alignments_ = 0;
  // agd::Buffer req_buf_;

  // genome pair denotes unique file to write to
  // we write results directly because storing them all in memory may be too
  // much

  struct LockedStream {
    LockedStream& operator=(LockedStream&& other) {
      out_stream = std::move(other.out_stream);
      return *this;
    }
    std::ofstream out_stream;
    absl::Mutex mu;
  };

  absl::node_hash_map<std::string, std::unique_ptr<std::ofstream>> file_map_;

  // track outstanding alignment requests so we know when we are done
  std::atomic_uint_fast64_t outstanding_{0};
  std::atomic_uint_fast64_t total_alignments_{0};
  std::atomic_uint_fast64_t total_matches_{0};

  const std::vector<Sequence>& sequences_;
  std::string output_dir_;

  size_t num_opened = 0;

  // lock file map
  absl::Mutex mu_;
};