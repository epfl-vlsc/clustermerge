
#include "controller.h"
#include <google/protobuf/text_format.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "src/agd/errors.h"
#include "src/common/all_all_executor.h"
#include "src/common/cluster_set.h"

using std::cout;
using std::string;
using std::thread;

void RemoveDuplicates(cmproto::ClusterSet& set) {
  absl::flat_hash_set<std::vector<size_t>> set_map;

  auto cluster_it = set.mutable_clusters()->begin();
  std::vector<size_t> cluster_set;
  while (cluster_it != set.mutable_clusters()->end()) {
    for (const auto& s : cluster_it->indexes()) {
      cluster_set.push_back(s);
    }
    std::sort(cluster_set.begin(), cluster_set.end());

    auto result = set_map.insert(std::move(cluster_set));
    if (!result.second) {
      cluster_it = set.mutable_clusters()->erase(cluster_it);
    } else {
      cluster_it++;
    }
    cluster_set.clear();
  }
}

// merge other into
void MergePartials(cmproto::ClusterSet& set, const cmproto::ClusterSet& other,
                   uint32_t original_size) {
  // relies on clusters in the sets being in the same order
  // where any new cluster is the last element
  // we do not "fully merge" any of the clusters in `set`

  // cout << "merging partial clusters ...\n";
  for (uint32_t i = 0; i < original_size; i++) {
    auto* mut_cluster = set.mutable_clusters(i);
    auto& cluster = other.clusters(i);
    /*string s;
    google::protobuf::TextFormat::PrintToString(*mut_cluster, &s);
    cout << "Merging partial: " << s << "\n";
    google::protobuf::TextFormat::PrintToString(cluster, &s);
    cout << "with " << s << "\n";*/
    if (cluster.fully_merged()) {
      mut_cluster->set_fully_merged(true);
    }

    for (auto index : cluster.indexes()) {
      if (std::find(mut_cluster->indexes().begin(),
                    mut_cluster->indexes().end(),
                    index) == mut_cluster->indexes().end()) {
        // doesnt exist, add
        // cout << "merger adding index to cluster\n";
        mut_cluster->add_indexes(index);
      }
    }
  }
  // if the partial merge generated a new cluster, add it to the set
  if (other.clusters_size() > original_size) {
    assert(original_size == other.clusters_size() - 1);
    auto* c = set.add_clusters();
    c->CopyFrom(other.clusters(other.clusters_size() - 1));
  }
  // cout << "done merging\n";
}

agd::Status Controller::Run(const Params& params,
                            const Parameters& aligner_params,
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
  outstanding_merges_ = total_merges;
  cout << "outstanding merges to complete: " << outstanding_merges_ << "\n";
  cout << "dup removal thresh is " << params.dup_removal_thresh << "\n";

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
  auto address = std::string("tcp://*:");
  auto response_queue_address =
      absl::StrCat(address, params.response_queue_port);
  auto request_queue_address = absl::StrCat(address, params.request_queue_port);

  context_ = zmq::context_t(1);
  context_sink_ = zmq::context_t(2);
  try {
    zmq_recv_socket_.reset(new zmq::socket_t(context_sink_, ZMQ_PULL));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PULL socket ");
  }

  try {
    zmq_send_socket_.reset(new zmq::socket_t(context_, ZMQ_PUSH));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PUSH socket ");
  }

  zmq_send_socket_->setsockopt(ZMQ_SNDHWM, 1);
  int val = zmq_send_socket_->getsockopt<int>(ZMQ_SNDHWM);
  cout << "snd hwm value is " << val << " \n";

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

  request_queue_.reset(
      new ConcurrentQueue<cmproto::MergeRequest>(params.queue_depth));
  response_queue_.reset(
      new ConcurrentQueue<cmproto::Response>(params.queue_depth));
  sets_to_merge_queue_.reset(
      new ConcurrentQueue<cmproto::ClusterSet>(sequences_.size()));

  request_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    cmproto::MergeRequest merge_request;
    int total_sent = 0;
    while (run_) {
      if (!request_queue_->pop(merge_request)) {
        continue;
      }
      auto size = merge_request.ByteSizeLong();
      zmq::message_t msg(size);
      /*cout << "pushing request of size " << size << " of type "
           << (merge_request.has_batch() ? "batch " : "partial ") << "\n";*/
      auto success = merge_request.SerializeToArray(msg.data(), size);
      if (!success) {
        cout << "Thread failed to serialize request protobuf!\n";
      }

      success = zmq_send_socket_->send(std::move(msg));
      if (!success) {
        cout << "Thread failed to send request over zmq!\n";
      }
      total_sent++;
    }

    cout << "Work queue thread ending. Total sent: " << total_sent << "\n";
  });

  response_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    cmproto::Response response;
    zmq::message_t msg;
    int total_received = 0;
    while (run_) {
      bool msg_received = zmq_recv_socket_->recv(&msg, ZMQ_NOBLOCK);
      // basically implements polling to avoid blocking recv calls
      // not ideal but whatever, this is a research project! :-D
      if (!msg_received) {
        continue;
      }
      total_received++;

      if (!response.ParseFromArray(msg.data(), msg.size())) {
        cout << "Failed to parse merge request protobuf!!\n";
        return;
      }
      /*cout << "parsed a zmq response with " << msg.size() << " bytes and "
           << response.set().clusters_size() << " clusters\n";*/

      response_queue_->push(response);
    }

    cout << "Work queue thread ending. Total received: " << total_received << "\n";
  });

  // partial mergers in this thread may need to be more fully parallelized to
  // prevent bottlenecks. May be required to use a different structure for
  // tracking partial mergers rather than the current map, which needs to be
  // locked
  worker_thread_ = thread([this, &params]() {
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
      // cout << "repsonse id is " << id << "\n";
      PartialMergeItem* partial_item;

      {
        absl::MutexLock l(&mu_);
        auto partial_it = partial_merge_map_.find(id);
        if (partial_it == partial_merge_map_.end()) {
          if (id != -1) {
            cout << "error thing in map was not -1\n";
            exit(0);
          }
          // cout << "pushing full result \n";
          if (response.set().clusters_size() > params.dup_removal_thresh) {
            RemoveDuplicates(*response.mutable_set());
          }
          sets_to_merge_queue_->push(std::move(response.set()));
          continue;
        } else {
          // do assign
          partial_item = &partial_it->second;
        }
        // cout << "preparing to merge partial result... \n";
        // its part of a larger merger
        // auto& partial_item = partial_it->second;
        // merge response_set into partial_item.partial_set
        MergePartials(partial_item->partial_set, response.set(),
                      partial_item->original_size);
      }

      // cout << "done\n";
      // check if this was the last one
      partial_item->num_received++;
      if (partial_item->num_expected == partial_item->num_received) {
        // go through  and delete any fully merged clusters
        auto cluster_it = partial_item->partial_set.mutable_clusters()->begin();
        while (cluster_it !=
               partial_item->partial_set.mutable_clusters()->end()) {
          if (cluster_it->fully_merged()) {
            cluster_it =
                partial_item->partial_set.mutable_clusters()->erase(cluster_it);
            continue;
          }
          cluster_it++;
        }

        if (partial_item->partial_set.clusters_size() >
            params.dup_removal_thresh) {
          RemoveDuplicates(partial_item->partial_set);
        }
        sets_to_merge_queue_->push(std::move(partial_item->partial_set));
        // remove partial it, its done now
        cout << "partial id " << id << " is complete\n";
        {
          absl::MutexLock l(&mu_);
          partial_merge_map_.erase(id);
        }
      }
    }
  });

  // dump all sequences in single cluster sets into the queue
  // for now assumes all of this will fit in memory
  // even a million sequences would just be a few MB
  for (const auto& s : sequences_) {
    cmproto::ClusterSet set;
    auto* c = set.add_clusters();
    c->add_indexes(s.ID());
    sets_to_merge_queue_->push(std::move(set));
  }

  auto t0 = std::chrono::high_resolution_clock::now();
  // just use 'this' thread to schedule outgoing work
  // take stuff from sets_to_merge_ and schedule the work in
  // batches or split partial merges for larger sets
  cmproto::MergeRequest request;
  std::vector<cmproto::ClusterSet> sets;
  sets.resize(2);
  while (outstanding_merges_ > 0) {
    sets[0].Clear();
    sets[1].Clear();
    if (!sets_to_merge_queue_->pop(sets[0])) {
      continue;
    }
    // get next so we have two to merge
    if (!sets_to_merge_queue_->pop(sets[1])) {
      return agd::errors::Internal(
          "error: did not get set for second to merge with.");
    }
    /*cout << "processing two sets ...\n";
    cout << "set one size: " << sets[0].clusters_size()
         << ", set two size: " << sets[1].clusters_size() << "\n\n";*/

    // form request, push to queue
    if (sets[0].clusters_size() < params.batch_size ||
        sets[1].clusters_size() < params.batch_size) {
      // create a batch, they are small
      // cout << "two sets are small, batching ...\n";
      cmproto::MergeBatch* batch = request.mutable_batch();
      uint32_t total_clusters =
          sets[0].clusters_size() + sets[1].clusters_size();

      auto* c = batch->add_sets();
      c->CopyFrom(sets[0]);
      c = batch->add_sets();
      c->CopyFrom(sets[1]);
      outstanding_merges_--;
      while (total_clusters < params.batch_size && outstanding_merges_ > 0) {
        // add more cluster sets to batch
        if (!sets_to_merge_queue_->pop(sets[0])) {
          return agd::errors::Internal(
              "ERROR: did not get set for second to merge with.");
        }

        c = batch->add_sets();
        c->CopyFrom(sets[0]);
        outstanding_merges_--;
        total_clusters += sets[0].clusters_size();
      }
      // cout << "batched " << batch->sets_size() << " sets\n";
      request.set_id(-1);
      // if the queue uses copy semantics im not sure how protobufs
      // with submessages will behave
      request_queue_->push(std::move(request));
      request.Clear();

    } else {
      // either set is large enough, split the computation into multiple
      // requests
      // cout << "splitting merger of two large sets into partial mergers\n";
      outstanding_merges_--;
      /*if (outstanding_merges_ == 0) {
        cout << "final two merging\n\n\n";
        for (int i = 0; i < sets[0].clusters_size(); i++) {
          cout << "cluster has " << sets[0].clusters(i).indexes_size() << "
      seqs\n";
        }
        cout << "\n";
        for (int i = 0; i < sets[1].clusters_size(); i++) {
          cout << "cluster has " << sets[1].clusters(i).indexes_size() << "
      seqs\n";
        }
      }*/
      // make a map entry for this multi-part request
      PartialMergeItem item;
      item.num_received = 0;
      // use the outstanding merges as id
      if (sets[0].clusters_size() < sets[1].clusters_size()) {
        sets[0].Swap(&sets[1]);
      }
      // each work item does a partial merge of one cluster in sets[0] to
      // all clusters in sets[1]
      item.num_expected = sets[0].clusters_size();
      item.partial_set.CopyFrom(sets[1]);
      item.original_size = sets[1].clusters_size();
      {
        absl::MutexLock l(&mu_);
        partial_merge_map_.insert_or_assign(outstanding_merges_, item);
      }

      cout << "pushing id " << outstanding_merges_ << "\n";
      for (const auto& c : sets[0].clusters()) {
        request.set_id(outstanding_merges_);
        cmproto::MergePartial* partial_request = request.mutable_partial();
        auto* cluster = partial_request->mutable_cluster();
        cluster->CopyFrom(c);
        auto* cluster_set = partial_request->mutable_set();
        cluster_set->CopyFrom(sets[1]);
        // cout << "pushing partial request with " <<
        // partial_request->set().clusters_size() << " clusters in set and ID: "
        // << request.id() << "\n";
        request_queue_->push(std::move(request));
        request.Clear();
      }
    }
    cout << "outstanding merges: " << outstanding_merges_ << "\n";
  }

  // we must now wait for the last results to come in
  // wait for worker thread to push last merged set
  // TODO add a timeout or something?
  cout << "done and waiting for final result...\n";
  // while (sets_to_merge_queue_->size() != 1);;

  cout << "scheduling final alignments on controller...\n";
  cmproto::ClusterSet final_set;
  sets_to_merge_queue_->peek(final_set);
  cout << "final set size is " << final_set.clusters_size() << " clusters\n";
  cout << "partial merge map size is " << partial_merge_map_.size() << "\n";
  auto t1 = std::chrono::high_resolution_clock::now();

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time: " << sec.count() << " seconds.\n";

  ClusterSet set(final_set, sequences_);
  set.DumpJson("dist_clusters.json");
  AllAllExecutor executor(std::thread::hardware_concurrency(), 500, &envs,
                          &aligner_params);
  executor.Initialize();
  set.ScheduleAlignments(&executor);
  executor.FinishAndOutput("dist_output_dir");

  cout << "clustering complete!! Joining threads ...\n";

  run_ = false;
  response_queue_->unblock();
  request_queue_->unblock();
  sets_to_merge_queue_->unblock();
  worker_thread_.join();
  request_queue_thread_.join();
  response_queue_thread_.join();

  cout << "All threads joined.\n";

  return agd::Status::OK();
}
