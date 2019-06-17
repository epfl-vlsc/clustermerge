
#pragma once

#include <thread>
#include "src/agd/status.h"
#include "src/common/concurrent_queue.h"
#include "src/common/params.h"
#include "src/common/sequence.h"
#include "src/comms/requests.h"
#include "src/dataset/dataset.h"
#include "zmq.hpp"

// interface to manage local worker
class Worker {
 public:
  Worker() {}

  struct Params {
    size_t num_threads;
    size_t queue_depth;
    absl::string_view controller_ip;
    int request_queue_port;
    int response_queue_port;
    int incomplete_request_queue_port;
    absl::string_view data_dir_path;
  };

  // main worker entry point
  // reuse agd status for error passing
  // requires num_threads, data_dir_path, controller_ip, push_port, pull_port
  agd::Status Run(const Params& params, const Parameters& aligner_params,
                  std::vector<std::unique_ptr<Dataset>>& datasets,
                  int* const signal_num);

  // signal handler
  agd::Status SignalHandler(int signal_num);

 private:
  zmq::context_t context_;
  std::unique_ptr<zmq::socket_t> zmq_recv_socket_;
  std::unique_ptr<zmq::socket_t> zmq_send_socket_;
  std::unique_ptr<zmq::socket_t> zmq_incomplete_request_socket_;
  // local input buffer
  // one thread receives zmq messages, decodes, puts work in work queue
  // all other threads do work and push to output buffer
  std::unique_ptr<ConcurrentQueue<zmq::message_t>> work_queue_;
  std::thread work_queue_thread_;
  bool wqt_signal_ = false;

  std::vector<std::thread> worker_threads_;
  bool worker_signal_ = false;

  // local output buffer
  // one thread encodes, sends results to controller
  std::unique_ptr<ConcurrentQueue<MarshalledResponse>> result_queue_;
  std::thread result_queue_thread_;
  bool rqt_signal_ = false;

  // local output buffer for incomplete requests
  // one thread which encodes and sends to the controller
  std::unique_ptr<ConcurrentQueue<MarshalledRequest>> incomplete_request_queue_;
  std::thread incomplete_request_queue_thread_;
  bool irqt_signal_ = false;

  std::vector<Sequence> sequences_;  // abs indexable sequences

  std::thread queue_measure_thread_;
  std::vector<long int> timestamps_;
  std::vector<size_t> queue_sizes_;

  // addresses
  std::string response_queue_address;
  std::string request_queue_address;
  std::string incomplete_request_queue_address;
};