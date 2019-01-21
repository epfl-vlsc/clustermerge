
#pragma once

#include <thread>
#include "src/agd/status.h"
#include "src/common/concurrent_queue.h"
#include "src/common/sequence.h"
#include "src/dataset/dataset.h"
#include "src/proto/requests.pb.h"
#include "src/proto/responses.pb.h"
#include "zmq.hpp"

// interface to manage local worker
class Worker {
 public:
  Worker() {}
  // main worker entry point
  // reuse agd status for error passing
  // requires num_threads, data_dir_path, controller_ip, push_port, pull_port
  agd::Status Run(size_t num_threads, size_t queue_depth,
                  const std::string& data_dir_path,
                  const std::string& controller_ip, int push_port,
                  int pull_port,
                  std::vector<std::unique_ptr<Dataset>>& datasets);

 private:
  zmq::context_t context_;
  std::unique_ptr<zmq::socket_t> zmq_recv_socket_;
  std::unique_ptr<zmq::socket_t> zmq_send_socket_;
  // local input buffer
  // one thread receives zmq messages, decodes, puts work in work queue
  // all other threads do work and push to output buffer
  std::unique_ptr<ConcurrentQueue<cmproto::MergeRequest>> work_queue_;
  std::thread work_queue_thread_;

  std::vector<std::thread> worker_threads_;

  // local output buffer
  // one thread encodes, sends results to controller
  std::unique_ptr<ConcurrentQueue<cmproto::Response>> result_queue_;
  std::thread result_queue_thread_;

  volatile bool run_ = true;

  std::vector<Sequence> sequences_; // abs indexable sequences
};