
#include "controller.h"
#include <iostream>
#include "src/agd/errors.h"

using std::cout;
using std::thread;

agd::Status Controller::Run(size_t num_threads, size_t queue_depth,
                            const std::string& controller_ip,
                            int request_queue_port, int response_queue_port,
                            std::vector<std::unique_ptr<Dataset>>& datasets) {
  // index all sequences
  agd::Status s = Status::OK();
  const char* data;
  size_t length;
  size_t id = 0;  // absolute ID
  for (auto& dataset : datasets) {
    size_t dataset_index = 0;
    s = dataset->GetNextRecord(&data, &length);
    while (s.ok()) {
      auto seq = Sequence(absl::string_view(data, length), dataset->Name(),
                          dataset->Size(), dataset_index, id);
      sequences_.push_back(std::move(seq));
      id++;
      dataset_index++;
      s = dataset->GetNextRecord(&data, &length);
    }
  }

  auto total_merges = sequences_.size() - 1;

  // connect to zmq queues
  auto address = std::string("tcp://*:");
  auto response_queue_address = absl::StrCat(address, response_queue_port);
  auto request_queue_address = absl::StrCat(address, request_queue_port);

  context_ = zmq::context_t(1);
  try {
    zmq_recv_socket_.reset(new zmq::socket_t(context_, ZMQ_PULL));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PULL socket ");
  }

  try {
    zmq_send_socket_.reset(new zmq::socket_t(context_, ZMQ_PUSH));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PUSH socket ");
  }

  try {
    zmq_recv_socket_->bind(response_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 response_queue_address);
  }

  try {
    zmq_send_socket_->bind(request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 request_queue_address);
  }

  request_queue_.reset(new ConcurrentQueue<cmproto::MergeRequest>(queue_depth));
  response_queue_.reset(new ConcurrentQueue<cmproto::Response>(queue_depth));

  request_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    cmproto::MergeRequest merge_request;
    while (run_) {
      if (!request_queue_->pop(merge_request)) {
        continue;
      }
      auto size = merge_request.ByteSizeLong();
      zmq::message_t msg(size);
      auto success = merge_request.SerializeToArray(msg.data(), size);
      if (!success) {
        cout << "Thread failed to serialize request protobuf!\n";
      }

      success = zmq_send_socket_->send(std::move(msg));
      if (!success) {
        cout << "Thread failed to send request over zmq!\n";
      }
    }

    cout << "Work queue thread ending.\n";
  });

  response_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    cmproto::Response response;
    zmq::message_t msg;
    while (run_) {
      zmq_recv_socket_->recv(&msg);

      if (!response.ParseFromArray(msg.data(), msg.size())) {
        cout << "Failed to parse merge request protobuf!!\n";
        return;
      }

      assert(outstanding_requests_.load() > 0);
      outstanding_requests_--;
      response_queue_->push(response);
    }

    cout << "Work queue thread ending.\n";
  });

  auto worker_func = [this]() {
    // read from result queue
    // if is a batch result and is small enough, push to WorkManager
    // if is partial result (ID will be
    // in merge map), lookup and merge with partial_set, if partial set
    // complete,
    //    push to WorkManager
    cmproto::Response response;
    while (run_) {
      if (!response_queue_->pop(response)) {
        continue;
      }

      auto id = response.id();
      auto partial_it = partial_merge_map_.find(id);
      if (partial_it != partial_merge_map_.end()) {
        // its part of a larger merger
        auto& partial_item = partial_it->second;
        // merge response_set into partial_item.partial_set
        // if complete now, manager.ProcessResponse(partial_item.partial_set)
        // remove partial_it
      } else {
        // it's a full result
        sets_to_merge_queue_->push(std::move(response.set()));
      }
      // outstanding_requests_--;
    }
  };

  auto scheduler_func = [this] () {
    cmproto::MergeRequest request;
    std::vector<cmproto::ClusterSet> sets;
    sets.resize(2);
    while (run_) {
      

    }
  };

  return agd::Status::OK();
}
