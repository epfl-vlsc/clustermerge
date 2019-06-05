
#pragma once

#include <thread>
#include <atomic>
#include "absl/container/node_hash_map.h"
#include "src/common/concurrent_queue.h"
#include "src/common/sequence.h"
#include "src/common/params.h"
#include "src/comms/requests.h"
#include "src/dataset/dataset.h"
#include "src/dist/partial_merge.h"
#include "zmq.hpp"

// one exists in a cluster
// coordinates work among workers, merges results
// currently using zmq push/pull to distribute work, but this may need
// to be replaced with a load balancing pattern
class Controller {
 public:
  struct Params {
    size_t num_threads;
    size_t queue_depth;
    absl::string_view controller_ip;
    int request_queue_port;
    int response_queue_port;
    int incomplete_request_queue_port;
    absl::string_view data_dir_path;
    uint32_t batch_size;
    uint32_t dup_removal_thresh;
    bool exclude_allall;
    int dataset_limit;
  };

  agd::Status Run(const Params& params, const Parameters& aligner_params,
                  std::vector<std::unique_ptr<Dataset>>& datasets);

 private:
  zmq::context_t context_;
  zmq::context_t context_sink_;
  zmq::context_t context_ir_sink_;  //ir stands for incomplete requests
  std::unique_ptr<zmq::socket_t> zmq_recv_socket_;
  std::unique_ptr<zmq::socket_t> zmq_send_socket_;
  std::unique_ptr<zmq::socket_t> zmq_incomp_req_socket_;

  // thread reads from zmq and puts into response queue
  std::unique_ptr<ConcurrentQueue<MarshalledResponse>> response_queue_;
  std::thread response_queue_thread_;

  // controller work threads
  // read from response queue
  // if is a batch result, push to sets_to_merge_queue_
  // if is partial result (ID will be
  // in merge map), lookup and merge with partial_set, if partial set complete,
  //    push ready to merge sets to sets_to_merge_queue
  std::unique_ptr<ConcurrentQueue<MarshalledClusterSet>> sets_to_merge_queue_;
  std::vector<std::thread> worker_threads_;

  // thread reads from request queue and pushes to zmq
  std::unique_ptr<ConcurrentQueue<MarshalledRequest>> request_queue_;
  // may need more than one thread .... we'll see
  std::thread request_queue_thread_;

  // controller reads incomplete requests into this queue
  std::unique_ptr<ConcurrentQueue<MarshalledRequest>> incomplete_request_queue_;
  std::thread incomplete_request_queue_thread_;

  std::vector<Sequence> sequences_;  // abs indexable sequences

  // indexed cluster and partial merge set to facilitate efficient 
  // parallel merging of remotely executed partial merge items

  // map structure to track incomplete partial merges
  struct PartialMergeItem {
    //cmproto::ClusterSet partial_set;
    PartialMergeItem() = default;
    PartialMergeItem(const PartialMergeItem& other) = delete;
    PartialMergeItem(PartialMergeItem&& other) {
      partial_set = std::move(other.partial_set);
      num_expected = other.num_expected;
      num_received.store(other.num_received.load());
    }
    PartialMergeItem& operator=(PartialMergeItem&& other) {
      partial_set = std::move(other.partial_set);
      num_expected = other.num_expected;
      num_received.store(other.num_received.load());
      return *this;
    }
    PartialMergeSet partial_set;
    uint32_t num_expected;
    std::atomic_uint_fast32_t num_received;
    //uint32_t original_size;
  };

  // we use a node map so that pointers remain stable, and we can reduce the
  // time spent in critical sections
  absl::node_hash_map<uint32_t, PartialMergeItem> partial_merge_map_;

  // uint32_t current_request_id_ = 0;

  volatile bool run_ = true;
  uint32_t outstanding_merges_ = 0;
  volatile bool outstanding_partial_ = false;
  absl::Mutex mu_;  // for partial_merge_map_
};
