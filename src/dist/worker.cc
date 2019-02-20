
#include "worker.h"
#include <deque>
#include <fstream>
#include "absl/strings/str_cat.h"
#include "json.hpp"
#include "merge_batch.h"
#include "src/common/aligner.h"
#include "src/common/alignment_environment.h"
#include "src/common/params.h"

using std::cout;
using std::string;
using std::thread;

agd::Status Worker::Run(const Params& params, const Parameters& aligner_params,
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

  // create envs, params
  string logpam_json_file = absl::StrCat(params.data_dir_path, "logPAM1.json");
  string all_matrices_json_file =
      absl::StrCat(params.data_dir_path, "all_matrices.json");

  std::ifstream logpam_stream(logpam_json_file);
  std::ifstream allmat_stream(all_matrices_json_file);

  if (!logpam_stream.good()) {
    return agd::errors::Internal("File ", logpam_json_file, " not found.");
  }

  if (!allmat_stream.good()) {
    return agd::errors::Internal("File ", all_matrices_json_file,
                                 " not found.");
  }

  json logpam_json;
  logpam_stream >> logpam_json;

  json all_matrices_json;
  allmat_stream >> all_matrices_json;

  AlignmentEnvironments envs;

  // initializing envs is expensive, so don't copy this
  cout << "Worker initializing environments from " << params.data_dir_path
       << "\n";
  envs.InitFromJSON(logpam_json, all_matrices_json);
  cout << "Done.\n";

  // done init envs

  // connect to zmq queues
  auto address = absl::StrCat("tcp://", params.controller_ip, ":");
  auto response_queue_address =
      absl::StrCat(address, params.response_queue_port);
  auto request_queue_address = absl::StrCat(address, params.request_queue_port);

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
    zmq_recv_socket_->connect(request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 request_queue_address);
  }

  try {
    zmq_send_socket_->connect(response_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 response_queue_address);
  }

  work_queue_.reset(
      new ConcurrentQueue<cmproto::MergeRequest>(params.queue_depth));
  result_queue_.reset(
      new ConcurrentQueue<cmproto::Response>(params.queue_depth));

  work_queue_thread_ = thread([this]() {
    // get msg from zmq
    // decode
    // put in work queue
    // repeat
    cmproto::MergeRequest merge_request;
    while (run_) {
      zmq::message_t msg;
      zmq_recv_socket_->recv(&msg);

      if (!merge_request.ParseFromArray(msg.data(), msg.size())) {
        cout << "Failed to parse merge request protobuf!!\n";
        return;
      }

      work_queue_->push(merge_request);
    }

    cout << "Work queue thread ending.\n";
  });

  auto worker_func = [this, &envs, &aligner_params]() {
    ProteinAligner aligner(&envs, &aligner_params);

    cmproto::MergeRequest request;
    std::deque<ClusterSet> sets_to_merge;
    while (run_) {
      if (!work_queue_->pop(request)) {
        continue;
      }

      // build cluster(s)
      // merge (do work)
      // encode result, put in queue

      if (request.has_batch()) {
        // cout << "its a batch\n";
        auto& batch = request.batch();
        for (size_t i = 0; i < batch.sets_size(); i++) {
          auto& set_proto = batch.sets(i);
          // construct cluster set from proto
          // TODO add ClusterSet constructor
          ClusterSet cs(set_proto, sequences_);
          sets_to_merge.push_back(std::move(cs));
        }
        // sets are constructed and in the queue, proceed to merge them
        // MergeSets(sets_to_merge)
        // cout << "merging a batch...\n";
        MergeBatch(sets_to_merge, &aligner);
        // queue now has final set
        auto& final_set = sets_to_merge[0];
        // encode to protobuf, push to queue
        // cout << "the merged set has " << final_set.Size() << " clusters\n";

        cmproto::Response response;
        auto* new_cs_proto = response.mutable_set();
        final_set.ConstructProto(new_cs_proto);

        response.set_id(request.id());
        result_queue_->push(response);
        sets_to_merge.clear();

      } else if (request.has_partial()) {
        // cout << "has partial\n";
        auto& partial = request.partial();
        // execute a partial merge, merge cluster into cluster set
        // do not remove any clusters, simply mark fully merged so
        // the controller can merge other partial merge requests

        // cout << "set has " << partial.set().clusters_size() << " clusters\n";
        ClusterSet cs(partial.set(), sequences_);
        Cluster c(partial.cluster(), sequences_);

        // cout << "merging cluster set with cluster\n";
        auto new_cs = cs.MergeCluster(c, &aligner);
        // cout << "cluster set now has " << new_cs.Size() << " clusters\n";

        cmproto::Response response;
        response.set_id(request.id());
        auto* new_cs_proto = response.mutable_set();
        new_cs.ConstructProto(new_cs_proto);

        result_queue_->push(response);

      } else {
        cout << "request was empty!!\n";
        return;
      }
    }
  };

  worker_threads_.reserve(params.num_threads);
  cout << "spinning up " << params.num_threads << " threads. \n";
  for (size_t i = 0; i < params.num_threads; i++) {
    worker_threads_.push_back(std::thread(worker_func));
  }

  result_queue_thread_ = thread([this]() {
    cmproto::Response response;
    while (run_) {
      if (!result_queue_->pop(response)) {
        continue;
      }
      auto size = response.ByteSizeLong();
      zmq::message_t msg(size);
      auto success = response.SerializeToArray(msg.data(), size);
      if (!success) {
        cout << "Thread failed to serialize response protobuf!\n";
      }

      success = zmq_send_socket_->send(std::move(msg));
      if (!success) {
        cout << "Thread failed to send response over zmq!\n";
      }
    }

    cout << "Work queue thread ending.\n";
  });

  timestamps_.reserve(100000);
  queue_sizes_.reserve(100000);

  queue_measure_thread_ = std::thread([this](){
      // get timestamp, queue size
      //cout << "queue measure thread starting ...\n";
      while(run_) {
        time_t result = std::time(nullptr);
        timestamps_.push_back(static_cast<long int>(result));
        queue_sizes_.push_back(work_queue_->size());
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        if (queue_sizes_.size() >= 1000000) {
          break;  // dont run forever ... 
        }
      }
      //cout << "queue measure thread finished\n";
    }
  );

  cout << "Worker running, press button to exit\n";
  //std::cin.get();

  run_ = false;
  result_queue_thread_.join();
  work_queue_thread_.join();
  queue_measure_thread_.join();
  for (auto& t : worker_threads_) {
    t.join();
  }
  
  // queue size stats
  std::vector<std::pair<size_t, size_t>> values;
  for (size_t i = 0; i < timestamps_.size(); i++) {
    values.push_back(std::make_pair(timestamps_[i], queue_sizes_[i]));
  }

  nlohmann::json j(values);

  std::cout << "dumping queue sizes ...\n";
  std::ofstream o("queue.json");

  o << std::setw(2) << j << std::endl;

  return agd::Status::OK();
}
