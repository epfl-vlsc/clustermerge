
#include "worker.h"
#include <sys/types.h>
#include <unistd.h>
#include <csignal>
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

agd::Status Worker::SignalHandler(int signal_num) {
  zmq::message_t msg;

  // terminate wqt and drain work_queue_
  wqt_signal_ = true;
  work_queue_thread_.join();
  while (work_queue_->size() != 0) {
    work_queue_->pop(msg);
    MarshalledRequest rq;
    rq.buf.AppendBuffer(reinterpret_cast<const char*>(msg.data()), msg.size());
    incomplete_request_queue_->push(std::move(rq));
  }
  assert(work_queue_->size() == 0);
  cout << "Work queue emptied.\n";

  // terminate all worker threads
  worker_signal_ = true;
  work_queue_->unblock();
  for (auto& t : worker_threads_) {
    t.join();
  }
  cout << "All worker threads have joined.\n";

  // terminating rqt
  // expecting size of set request queue == 0 (since all workers will have quit)
  srt_signal_ = true;
  set_request_queue_->unblock();
  set_request_thread_.join();
  assert(set_request_queue_->size() == 0);
  cout << "Set request queue emptied and thread terminated.\n";

  // terminate irqt and drain incomplete_request_queue_
  irqt_signal_ = true;
  incomplete_request_queue_->unblock();
  incomplete_request_queue_thread_.join();
  assert(incomplete_request_queue_->size() == 0);
  cout << "Incomplete request queue emptied\n";

  // terminate rqt and drain result_queue_
  rqt_signal_ = true;
  result_queue_->unblock();
  result_queue_thread_.join();
  assert(result_queue_->size() == 0);
  cout << "Resuly queue emptied.\n";

  try {
    zmq_recv_socket_->disconnect(request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not disconnect to zmq at ",
                                 request_queue_address);
  }

  try {
    zmq_send_socket_->disconnect(response_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not disconnect to zmq at ",
                                 response_queue_address);
  }

  try {
    zmq_set_request_socket_->disconnect(set_request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not disconnect to zmq at ",
                                 set_request_queue_address);
  }

  try {
    zmq_incomplete_request_socket_->disconnect(
        incomplete_request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not disconnect to zmq at ",
                                 incomplete_request_queue_address);
  }

  cout << "Zmq's disconnected.\n";

  cout << "Quitting nicely..!\n";
  // exit(signal_num);
  return agd::Status::OK();
}

agd::Status Worker::Run(const Params& params, const Parameters& aligner_params,
                        std::vector<std::unique_ptr<Dataset>>& datasets,
                        int* const signal_num) {
  // print the process id
  pid_t pid = getpid();
  cout << "Process id: " << pid << std::endl;

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
  string blosum_json_file = absl::StrCat(params.data_dir_path, "BLOSUM62.json");

  std::ifstream logpam_stream(logpam_json_file);
  std::ifstream blosum_stream(blosum_json_file);
  std::ifstream allmat_stream(all_matrices_json_file);

  if (!logpam_stream.good()) {
    return agd::errors::Internal("File ", logpam_json_file, " not found.");
  }

  if (!blosum_stream.good()) {
    return agd::errors::Internal("File ", blosum_json_file, " not found.");
  }

  if (!allmat_stream.good()) {
    return agd::errors::Internal("File ", all_matrices_json_file,
                                 " not found.");
  }

  json logpam_json;
  logpam_stream >> logpam_json;

  json blosum_json;
  blosum_stream >> blosum_json;

  json all_matrices_json;
  allmat_stream >> all_matrices_json;

  AlignmentEnvironments envs;

  // initializing envs is expensive, so don't copy this
  cout << "Worker initializing environments from " << params.data_dir_path
       << "\n";
  envs.InitFromJSON(logpam_json, all_matrices_json, aligner_params.min_score);
  if (aligner_params.use_blosum) {
    envs.UseBlosum(blosum_json, aligner_params.min_score);
  }
  cout << "Done.\n";

  // done init envs

  // connect to zmq queues
  auto address = absl::StrCat("tcp://", params.controller_ip, ":");
  response_queue_address = absl::StrCat(address, params.response_queue_port);
  request_queue_address = absl::StrCat(address, params.request_queue_port);
  incomplete_request_queue_address =
      absl::StrCat(address, params.incomplete_request_queue_port);
  set_request_queue_address = absl::StrCat(address, params.set_request_port);

  context_ = zmq::context_t(1);
  // worker is the requester --> requests for work items
  try {
    zmq_recv_socket_.reset(new zmq::socket_t(context_, ZMQ_REQ));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq REQ socket ");
  }

  try {
    zmq_send_socket_.reset(new zmq::socket_t(context_, ZMQ_PUSH));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PUSH socket ");
  }

  try {
    zmq_set_request_socket_.reset(new zmq::socket_t(context_, ZMQ_REQ));
  } catch (...) {
    return agd::errors::Internal(
        "Could not create zmq REQ socket -- large pm ");
  }

  try {
    zmq_incomplete_request_socket_.reset(new zmq::socket_t(context_, ZMQ_PUSH));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq INCOMPLETE REQ socket ");
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

  try {
    zmq_set_request_socket_->connect(set_request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 set_request_queue_address);
  }

  try {
    zmq_incomplete_request_socket_->connect(
        incomplete_request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 incomplete_request_queue_address);
  }

  zmq_recv_socket_->setsockopt(ZMQ_SNDHWM, 5);
  int val = zmq_recv_socket_->getsockopt<int>(ZMQ_SNDHWM);
  cout << "snd hwm value is " << val << " \n";

  work_queue_.reset(new ConcurrentQueue<zmq::message_t>(params.queue_depth));
  result_queue_.reset(
      new ConcurrentQueue<MarshalledResponse>(params.queue_depth));
  incomplete_request_queue_.reset(
      new ConcurrentQueue<MarshalledRequest>(params.queue_depth));
  set_request_queue_.reset(
      new ConcurrentQueue<std::pair<int, MultiNotification*>>(
          params.num_threads));

  work_queue_thread_ = thread([this]() {
    // get msg from zmq
    // decode
    // put in work queue
    // repeat
    // cmproto::MergeRequest merge_request;
    while (!wqt_signal_) {
      zmq::message_t msg;
      memcpy(msg.data(), "Send-Work", 9);
      // std::cout << "Requesting for work..\n";
      zmq_recv_socket_->send(msg);
      bool msg_received = zmq_recv_socket_->recv(&msg);  //, ZMQ_NOBLOCK);
      if (!msg_received) {
        continue;
      }
      // cout << "Received work item.\n";
      work_queue_->push(std::move(msg));
      // cout << "Pushing into work_queue_\n";
    }

    cout << "Work queue thread ending.\n";
  });

  set_request_thread_ = thread([this]() {
    while (!srt_signal_) {
      std::pair<int, MultiNotification*> pr;
      if (!set_request_queue_->pop(pr)) {
        continue;
      }
      int id = pr.first;
      // if its not in the cache
      mu_.Lock();
      if (set_map_.find(id) == set_map_.end()) {
        mu_.Unlock();

        zmq::message_t msg(&id, sizeof(int));
        bool success = zmq_set_request_socket_->send(msg);
        if (!success) {
          cout << "Failed to send the id over zmq -- worker func.\n";
        }
        success = zmq_set_request_socket_->recv(&msg);
        if (!success) {
          cout << "Requested set was not received.\n";
        }

        // cout << "Received message size --> " << msg.size() << std::endl;
        // copy to put in the map
        agd::Buffer buf(msg.size());
        buf.AppendBuffer(reinterpret_cast<char*>(msg.data()), msg.size());
        
        //build the offset vector
        MarshalledClusterSetView set(buf.data());
        std::vector<size_t> offsets;
        set.Offsets(offsets);
        std::pair<agd::Buffer, std::vector<size_t>> set_and_offsets  = {std::move(buf), offsets};

        mu_.Lock();
        set_map_[id] = std::move(set_and_offsets);
        mu_.Unlock();
        cout << "Fetched set: [" << id << "]\n";
      } else {
        mu_.Unlock();
        // cout << "Communicating thread --> set already present: " << id <<
        // "\n";
      }

      pr.second->Notify();
    }
  });

  int nPartial = 0;
  int nJoined = 0;
  auto worker_func = [this, &envs, &aligner_params, &nPartial, &nJoined]() {
    ProteinAligner aligner(&envs, &aligner_params);
    // cmproto::MergeRequest request;
    std::deque<ClusterSet> sets_to_merge;
    while (!worker_signal_) {
      zmq::message_t msg;
      if (!work_queue_->pop(msg)) {
        continue;
      }
      MarshalledRequestView request(reinterpret_cast<char*>(msg.data()),
                                    msg.size());
      // build cluster(s)
      // merge (do work)
      // encode result, put in queue
      // cout << "Got a request: " << request.ID() << std::endl;
      if (request.Type() == RequestType::Batch) {
        // cout << "its a batch, request ID is " << request.ID() << " \n";
        // auto& batch = request.batch();
        MarshalledClusterSetView clusterset;
        while (request.NextClusterSet(&clusterset)) {
          // construct cluster set from proto
          // TODO add ClusterSet constructor
          // cout << "marshalled cluster set has " << clusterset.NumClusters()
          // << " clusters\n";
          ClusterSet cs(clusterset, sequences_);
          sets_to_merge.push_back(std::move(cs));
        }
        // sets are constructed and in the queue, proceed to merge them
        // MergeSets(sets_to_merge)
        // cout << "merging a batch...\n";
        MergeBatch(sets_to_merge, &aligner);
        // queue now has final set
        // cout << "--\n";
        auto& final_set = sets_to_merge[0];
        // encode to protobuf, push to queue
        // cout << "the merged set has " << final_set.Size() << " clusters\n";

        // final_set.ConstructProto(new_cs_proto);
        MarshalledResponse response;
        final_set.BuildMarshalledResponse(request.ID(), &response);
        MarshalledClusterSetView view;
        view = response.Set();
        // cout << "final set has " << view.NumClusters() << " clusters.\n";

        // cout << "Pushing to result queue.\n";
        result_queue_->push(std::move(response));
        sets_to_merge.clear();
        // cout << "Pushed to result queue.\n";
      } else if (request.Type() == RequestType::Partial) {
        nPartial++;
        // cout << "Its a partial request\n";
        // auto& partial = request.partial();
        // execute a partial merge, merge cluster into cluster set
        // do not remove any clusters, simply mark fully merged so
        // the controller can merge other partial merge requests
        MarshalledClusterSetView set;
        MarshalledClusterView cluster;
        request.ClusterAndSet(&set, &cluster);

        // cout << "set has " << set.NumClusters() << " clusters\n";
        ClusterSet cs(set, sequences_);
        Cluster c(cluster, sequences_);

        // cout << "merging cluster set with cluster\n";
        auto new_cs = cs.MergeCluster(c, &aligner, worker_signal_);
        if (worker_signal_) {
          MarshalledRequest rq;
          rq.buf.AppendBuffer(reinterpret_cast<const char*>(msg.data()),
                              msg.size());
          incomplete_request_queue_->push(std::move(rq));
          std ::cout << "Pushing into incomplete queue with ID: "
                     << request.ID() << std::endl;
          continue;
        }
        // cout << "cluster set now has " << new_cs.Size() << " clusters\n";
        assert(set.NumClusters() <= new_cs.Size());
        MarshalledResponse response;
        new_cs.BuildMarshalledResponse(request.ID(), &response);
        assert(response.Set().NumClusters() == new_cs.Size());
        result_queue_->push(std::move(response));

      } else if (request.Type() == RequestType::LargePartial) {
        // extract id and cluster
        int id = request.ID();
        // cout << "Received a large partial request with id " << id << "\n";
        MarshalledClusterView cluster;
        request.Cluster(&cluster);

        MarshalledClusterSet set;
        // search in map
        mu_.Lock();
        auto it = set_map_.find(id);
        if (it == set_map_.end()) {
          mu_.Unlock();
          MultiNotification n;
          // cout << "Worker thread --> [" << id << "] Set not cached.
          // Requesting...\n";
          set_request_queue_->push(std::make_pair(id, &n));
          n.SetMinNotifies(1);
          n.WaitForNotification();
          mu_.Lock();
          it = set_map_.find(id);
          assert(it != set_map_.end());
          // cout << "Worker thread --> [" << id << "] Received set.\n";
        }

        set.buf.AppendBuffer(it->second.first.data(), it->second.first.size());
        mu_.Unlock();

        ClusterSet cs(set, sequences_);
        Cluster c(cluster, sequences_);

        // same code as Partial Merge
        auto new_cs = cs.MergeCluster(c, &aligner, worker_signal_);
        if (worker_signal_) {
          MarshalledRequest rq;
          rq.buf.AppendBuffer(reinterpret_cast<const char*>(msg.data()),
                              msg.size());
          incomplete_request_queue_->push(std::move(rq));
          std ::cout << "Pushing into incomplete queue with ID: "
                     << request.ID() << std::endl;
          continue;
        }
        // cout << "cluster set now has " << new_cs.Size() << " clusters\n";
        assert(set.NumClusters() <= new_cs.Size());
        MarshalledResponse response;
        new_cs.BuildMarshalledResponse(request.ID(), &response);
        assert(response.Set().NumClusters() == new_cs.Size());
        result_queue_->push(std::move(response));

      } else if(request.Type() == RequestType::SubLargePartial) {
        int start_index, end_index, cluster_index, id;
        id = request.ID();
        MarshalledClusterView cluster;
        request.IndexesAndCluster(&start_index, &end_index, &cluster_index, &cluster);
        cout << "Got request with: " << start_index << " " << end_index << " " << cluster_index << " \n";      

        // search in map
        mu_.Lock();
        auto it = set_map_.find(id);
        if (it == set_map_.end()) {
          mu_.Unlock();
          MultiNotification n;
          // cout << "Worker thread --> [" << id << "] Set not cached.
          // Requesting...\n";
          set_request_queue_->push(std::make_pair(id, &n));
          n.SetMinNotifies(1);
          n.WaitForNotification();
          mu_.Lock();
          it = set_map_.find(id);
          assert(it != set_map_.end());
          // cout << "Worker thread --> [" << id << "] Received set.\n";
        }

        MarshalledClusterSetView set(it->second.first.data());
        ClusterSet cs(set, it->second.second, start_index, end_index, sequences_);
        mu_.Unlock();

        Cluster c(cluster, sequences_);

        // same code as Partial Merge
        auto new_cs = cs.MergeCluster(c, &aligner, worker_signal_);
        if (worker_signal_) {
          MarshalledRequest rq;
          rq.buf.AppendBuffer(reinterpret_cast<const char*>(msg.data()),
                              msg.size());
          incomplete_request_queue_->push(std::move(rq));
          std ::cout << "Pushing into incomplete queue with ID: "
                     << request.ID() << std::endl;
          continue;
        }
        // cout << "cluster set now has " << new_cs.Size() << " clusters\n";
        assert(set.NumClusters() <= new_cs.Size());
        
        //this stuff needs some change, we need to a way to know what response it is
        MarshalledResponse response;
        new_cs.BuildMarshalledResponse(request.ID(), &response);
        assert(response.Set().NumClusters() == new_cs.Size());
        result_queue_->push(std::move(response));

      } else {
        cout << "request was not any type!!!!!!\n";
        return;
      }
    }
    nJoined++;
    cout << "Ending thread [" << nJoined << "]\n";
  };

  worker_threads_.reserve(params.num_threads);
  cout << "spinning up " << params.num_threads << " threads. \n";
  for (size_t i = 0; i < params.num_threads; i++) {
    worker_threads_.push_back(std::thread(worker_func));
  }

  result_queue_thread_ = thread([this, &params]() {
    while (!(rqt_signal_ && result_queue_->size() == 0)) {
      MarshalledResponse response;

      if (!result_queue_->pop(response)) {
        continue;
      }
      bool success = zmq_send_socket_->send(std::move(response.msg));
      if (!success) {
        cout << "Thread failed to send response over zmq!\n";
      }
    }

    cout << "Result queue thread ending.\n";
  });

  incomplete_request_queue_thread_ = thread([this, &params]() {
    auto free_func = [](void* data, void* hint) {
      delete[] reinterpret_cast<char*>(data);
    };
    while (!(irqt_signal_ && incomplete_request_queue_->size() == 0)) {
      MarshalledRequest request;

      if (!incomplete_request_queue_->pop(request)) {
        continue;
      }

      auto size = request.buf.size();
      zmq::message_t msg(request.buf.release_raw(), size, free_func, NULL);
      bool success = zmq_incomplete_request_socket_->send(std::move(msg));
      if (!success) {
        cout << "INCOMP REQ Thread failed to send response over zmq!\n";
      }
    }

    cout << "Incomplete request queue thread ending.\n";
  });

  /*timestamps_.reserve(100000);
  queue_sizes_.reserve(100000);*/

  // queue_measure_thread_ = std::thread([this](){
  //     // get timestamp, queue size
  //     //cout << "queue measure thread starting ...\n";
  //     while(run_) {
  //       time_t result = std::time(nullptr);
  //       timestamps_.push_back(static_cast<long int>(result));
  //       queue_sizes_.push_back(work_queue_->size());
  //       std::this_thread::sleep_for(std::chrono::milliseconds(500));
  //       if (queue_sizes_.size() >= 1000000) {
  //         break;  // dont run forever ...
  //       }
  //     }
  //     cout << "queue measure thread finished\n";
  //   }
  // );

  while (!(*signal_num)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10 * 100));
  }
  cout << "Response [signal_num] value change received.\n";
  return SignalHandler(*signal_num);

  /*cout << "Worker running, press button to exit\n";
  // std::cin.get();

  cout << "joining threads ...\n";

  // incomplete_request_queue_thread_.join();
  // result_queue_thread_.join();

  // queue size stats
  std::vector<std::pair<size_t, size_t>> values;
  for (size_t i = 0; i < timestamps_.size(); i++) {
    values.push_back(std::make_pair(timestamps_[i], queue_sizes_[i]));
  }

  nlohmann::json j(values);

  std::cout << "dumping queue sizes ...\n";
  std::ofstream o("queue.json");

  o << std::setw(2) << j << std::endl;*/

  return agd::Status::OK();
}
