
#include <fstream>
#include <deque>
#include "absl/strings/str_cat.h"
#include "json.hpp"
#include "src/common/alignment_environment.h"
#include "src/common/aligner.h"
#include "src/common/params.h"
#include "worker.h"
#include "merge_batch.h"

using std::cout;
using std::string;
using std::thread;

agd::Status Worker::Run(size_t num_threads, size_t queue_depth,
                        const std::string& data_dir_path,
                        const std::string& controller_ip, int request_queue_port,
                        int response_queue_port, std::vector<std::unique_ptr<Dataset>>& datasets) {
  // index all sequences
  agd::Status s = Status::OK();
  const char* data;
  size_t length;
  size_t id = 0; // absolute ID
  for (auto& dataset : datasets) {
    size_t dataset_index = 0;
    s = dataset->GetNextRecord(&data, &length);
    while (s.ok()) {
      auto seq = Sequence(absl::string_view(data, length), dataset->Name(), dataset->Size(), dataset_index, id);
      sequences_.push_back(std::move(seq));
      id++;
      dataset_index++;
      s = dataset->GetNextRecord(&data, &length);
    }
  }

  // create envs, params
  string logpam_json_file = data_dir_path + "logPAM1.json";
  string all_matrices_json_file = data_dir_path + "all_matrices.json";

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
  cout << "Worker initializing environments from " << data_dir_path << "\n";
  envs.InitFromJSON(logpam_json, all_matrices_json);
  cout << "Done.\n";

  // done init envs

  Parameters params;  // using default params for now

  // connect to zmq queues
  auto address = absl::StrCat("tcp://", controller_ip, ":");
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
    zmq_recv_socket_->connect(request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ", request_queue_address);
  }
  
  try {
    zmq_send_socket_->connect(response_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ", response_queue_address);
  }

  work_queue_.reset(new ConcurrentQueue<cmproto::MergeRequest>(queue_depth));
  result_queue_.reset(new ConcurrentQueue<cmproto::Response>(queue_depth));

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

  auto worker_func = [this, &envs, &params]() {

    ProteinAligner aligner(&envs, &params);

    cmproto::MergeRequest request;
    std::deque<ClusterSet> sets_to_merge;
    while(run_) {
      if (!work_queue_->pop(request)) {
        continue;
      }

      // build cluster(s)
      // merge (do work)
      // encode result, put in queue

      if (request.has_batch()) {
        cout << "its a batch\n";
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
        MergeBatch(sets_to_merge, &aligner);
        // queue now has final set
        auto& final_set = sets_to_merge[0];
        // encode to protobuf, push to queue
        
        cmproto::Response response;
        auto* new_cs_proto = response.mutable_set();
        final_set.ConstructProto(new_cs_proto);

        result_queue_->push(response);

      } else if (request.has_partial()) {
        cout << "has partial\n";
        auto& partial = request.partial();
        // execute a partial merge, merge cluster into cluster set
        // do not remove any clusters, simply mark fully merged so 
        // the controller can merge other partial merge requests
         ClusterSet cs(partial.set(), sequences_);
         Cluster c(partial.cluster(), sequences_);

         cs.MergeCluster(c, &aligner);

         cmproto::Response response;
         auto* new_cs_proto = response.mutable_set();
         cs.ConstructProto(new_cs_proto);

         result_queue_->push(response);

      } else {
        cout << "request was empty!!\n";
        return;
      }
    }
  };

  worker_threads_.reserve(num_threads);
  for (size_t i = 0; i < num_threads; i++) {
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

  cout << "Worker running, press ctrl-C to exit\n";

  result_queue_thread_.join();
  work_queue_thread_.join();
  for (auto& t : worker_threads_) {
    t.join();
  }

  return agd::Status::OK();

}