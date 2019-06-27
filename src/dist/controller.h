
#pragma once

#include <atomic>
#include <thread>
#include "absl/container/node_hash_map.h"
#include "src/common/concurrent_queue.h"
#include "src/common/params.h"
#include "src/common/sequence.h"
#include "src/comms/requests.h"
#include "src/common/cluster_set.h"
#include "src/dataset/dataset.h"
#include "src/dist/partial_merge.h"
#include "zmq.hpp"

// one exists in a cluster
// coordinates work among workers, merges results
// currently using zmq push/pull to distribute work, but this may need
// to be replaced with a load balancing pattern
class Controller {
 public:
  Controller(nlohmann::json dataset_json_obj,
    std::vector<std::unique_ptr<Dataset>>& datasets_old) {
    const char* data_old;
    size_t size_old;
    uint32_t id_old = 0;

    //std::vector<Sequence> old_sequences;

    for (auto& dataset_old : datasets_old) {
      //std::cout << "Parsing dataset " << dataset_old->Name() << " ...\n";

      auto s_old = dataset_old->GetNextRecord(&data_old, &size_old);
      uint32_t genome_index_old = 0;
      while (s_old.ok()) {
        // cout << "Adding sequence id " << id << "\n";
        if (size_old > 60000) {
          std::cout << "over size " << size_old << "\n";
          exit(0);
        }

        Sequence seq(absl::string_view(data_old, size_old), dataset_old->Name(),
                    dataset_old->Size(), genome_index_old++, id_old++);

        sequences_.push_back(std::move(seq));
        s_old = dataset_old->GetNextRecord(&data_old, &size_old);
      }
    }

    for (const auto& cluster : dataset_json_obj["clusters"]) {
      Cluster c;
      for (const auto& seq : cluster) {
        int abs_index = seq["AbsoluteIndex"];
        c.AddSequence(sequences_[abs_index]);
      }
      old_set_.AddCluster(c);
    }

    absolute_id_ = id_old;  
    std::cout << "Existing clusters: " << old_set_.Size() << "\n";
  }

  Controller() = default;

  struct Params {
    size_t num_threads;
    size_t queue_depth;
    size_t max_set_size;
    absl::string_view controller_ip;
    int request_queue_port;
    int response_queue_port;
    int incomplete_request_queue_port;
    int set_request_port;
    absl::string_view data_dir_path;
    uint32_t batch_size;
    uint32_t dup_removal_thresh;
    bool exclude_allall;
    int dataset_limit;
    long int checkpoint_interval;
    absl::string_view checkpoint_dir;
  };

  agd::Status Run(const Params& params, const Parameters& aligner_params,
                  std::vector<std::unique_ptr<Dataset>>& datasets, std::vector<std::string>& dataset_names);

 private:
  zmq::context_t context_;
  zmq::context_t context_sink_;
  zmq::context_t context_ir_sink_;  // ir stands for incomplete requests
  zmq::context_t context_set_request_;
  std::unique_ptr<zmq::socket_t> zmq_recv_socket_;
  std::unique_ptr<zmq::socket_t> zmq_send_socket_;
  std::unique_ptr<zmq::socket_t> zmq_incomplete_request_socket_;
  std::unique_ptr<zmq::socket_t> zmq_set_request_socket_;

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

  // thread to send partial merge sets
  std::thread set_request_thread_;

  std::vector<Sequence> sequences_;  // abs indexable sequences

  long int checkpoint_timer_;
  // indexed cluster and partial merge set to facilitate efficient

  // parallel merging of remotely executed partial merge items

  // map structure to track incomplete partial merges
  struct PartialMergeItem {
    // cmproto::ClusterSet partial_set;
    PartialMergeItem() = default;
    PartialMergeItem(const PartialMergeItem& other) = delete;
    PartialMergeItem(PartialMergeItem&& other) {
      partial_set = std::move(other.partial_set);
      num_expected = other.num_expected;
      num_received.store(other.num_received.load());
      marshalled_set_buf = std::move(other.marshalled_set_buf);
    }
    PartialMergeItem& operator=(PartialMergeItem&& other) {
      partial_set = std::move(other.partial_set);
      num_expected = other.num_expected;
      num_received.store(other.num_received.load());
      marshalled_set_buf = std::move(other.marshalled_set_buf);
      return *this;
    }
    PartialMergeSet partial_set;
    uint32_t num_expected;
    std::atomic_uint_fast32_t num_received;
    // holds the MarshalledClusterSet responding to set requests
    agd::Buffer marshalled_set_buf = agd::Buffer(2048, 512);
  };

  // we use a node map so that pointers remain stable, and we can reduce the
  // time spent in critical sections
  absl::node_hash_map<uint32_t, PartialMergeItem> partial_merge_map_;

  // uint32_t current_request_id_ = 0;

  //holds the preclustered result if json is passed
  ClusterSet old_set_;
  size_t absolute_id_ = 0;  //for syncing with preclustered sequences 
  
  volatile bool run_ = true;
  uint32_t outstanding_merges_ = 0;
  volatile bool outstanding_partial_ = false;
  absl::Mutex mu_;  // for partial_merge_map_
};
