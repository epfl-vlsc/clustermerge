
#pragma once

#include <thread>
#include "absl/container/flat_hash_map.h"
#include "src/common/concurrent_queue.h"
#include "src/common/sequence.h"
#include "src/dataset/dataset.h"
#include "src/proto/requests.pb.h"
#include "src/proto/responses.pb.h"
#include "zmq.hpp"

// one exists in a cluster
// coordinates work among workers, merges results
// currently using zmq push/pull to distribute work, but this may need 
// to be replaced with a load balancing pattern
class Controller {
 public:
  agd::Status Run(size_t num_threads, size_t queue_depth, const std::string& controller_ip,
                  int request_queue_port, int response_queue_port,
                  const std::string& data_dir_path,
                  std::vector<std::unique_ptr<Dataset>>& datasets);

 private:

  zmq::context_t context_;
  std::unique_ptr<zmq::socket_t> zmq_recv_socket_;
  std::unique_ptr<zmq::socket_t> zmq_send_socket_;

  // thread reads from zmq and puts into response queue
  std::unique_ptr<ConcurrentQueue<cmproto::Response>> response_queue_;
  std::thread response_queue_thread_;
  
  // controller work threads
  // read from response queue
  // if is a batch result and is small enough, push to WorkManager
  // if is partial result (ID will be
  // in merge map), lookup and merge with partial_set, if partial set complete,
  //    push ready to merge sets to sets_to_merge_queue 
  std::unique_ptr<ConcurrentQueue<cmproto::ClusterSet>> sets_to_merge_queue_;
  std::thread worker_thread_;

  // merge scheduler thread
  // reads sets from sets to merge queue until a large enough batch is ready
  // if 2 large sets, split into multiple requests
  // form MergeRequest, push to request queue
  // might be bottleneck, we'll see
  std::thread merge_scheduler_thread_;

  // thread reads from request queue and pushes to zmq
  std::unique_ptr<ConcurrentQueue<cmproto::MergeRequest>> request_queue_;
  // may need more than one thread .... we'll see
  std::thread request_queue_thread_;

  std::vector<Sequence> sequences_;  // abs indexable sequences

  // map structure to track incomplete partial merges
  struct PartialMergeItem {
    cmproto::ClusterSet partial_set;
    uint32_t num_expected;
    uint32_t num_received;
    uint32_t original_size;
  };

  absl::flat_hash_map<uint32_t, PartialMergeItem> partial_merge_map_;

  //uint32_t current_request_id_ = 0;
  uint32_t batch_size_ = 100;

  volatile bool run_ = true;
  uint32_t outstanding_merges_ = 0;
  absl::Mutex mu_; // for partial_merge_map_

};