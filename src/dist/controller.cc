
#include "controller.h"
//#include <google/protobuf/text_format.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "src/agd/errors.h"
#include "src/common/all_all_executor.h"
#include "src/common/cluster_set.h"
#include <unistd.h>

using std::cout;
using std::string;
using std::thread;

/*class PartialMergeSet {
  public:
    void MergeClusterSet(const cmproto::ClusterSet& set);
    void BuildClusterSetProto(cmproto::ClusterSet* set);
  private:
    std::vector<IndexedCluster> clusters_;
    // lock to add new clusters
    absl::Mutex mu_;
};*/

/*void RemoveDuplicates(cmproto::ClusterSet& set) {
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
}*/

// merge other into
/*void MergePartials(cmproto::ClusterSet& set, const cmproto::ClusterSet& other,
                   uint32_t original_size) {
  // relies on clusters in the sets being in the same order
  // where any new cluster is the last element
  // we do not "fully merge" any of the clusters in `set`

  // cout << "merging partial clusters ...\n";
  for (uint32_t i = 0; i < original_size; i++) {
    auto* mut_cluster = set.mutable_clusters(i);
    auto& cluster = other.clusters(i);
    //string s;
    //google::protobuf::TextFormat::PrintToString(*mut_cluster, &s);
    //cout << "Merging partial: " << s << "\n";
    //google::protobuf::TextFormat::PrintToString(cluster, &s);
    //cout << "with " << s << "\n";
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
}*/

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
  if (params.dataset_limit > 0) {
    outstanding_merges_ = params.dataset_limit - 1;
  }
  cout << "outstanding merges to complete: " << outstanding_merges_ << "\n";
  cout << "dup removal thresh is " << params.dup_removal_thresh << "\n";
  cout << "Using " << params.num_threads << " threads to merge partials\n";

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
  envs.InitFromJSON(logpam_json, all_matrices_json, aligner_params.min_score);
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

  zmq_send_socket_->setsockopt(ZMQ_SNDHWM, 5);
  int val = zmq_send_socket_->getsockopt<int>(ZMQ_SNDHWM);
  cout << "snd hwm value is " << val << " \n";

  request_queue_.reset(
      new ConcurrentQueue<MarshalledRequest>(params.queue_depth));
  response_queue_.reset(
      new ConcurrentQueue<MarshalledResponse>(params.queue_depth));
  sets_to_merge_queue_.reset(
      new ConcurrentQueue<MarshalledClusterSet>(sequences_.size()));

  request_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    MarshalledRequest merge_request;
    int total_sent = 0;
    auto free_func = [](void* data, void* hint) { delete [] reinterpret_cast<char*>(data); };

    while (run_) {
      if (!request_queue_->pop(merge_request)) {
        continue;
      }
      auto size = merge_request.buf.size();
      // inline message_t(void *data_, size_t size_, free_fn *ffn_, void *hint_
      // = NULL) release the buf pointer directly to avoid an additional copy
      // here
      zmq::message_t msg(merge_request.buf.release_raw(), size, free_func,
                         NULL);
      /*cout << "pushing request of size " << size << " of type "
           << (merge_request.has_batch() ? "batch " : "partial ") << "\n";*/

      bool success = zmq_send_socket_->send(std::move(msg));
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
    int total_received = 0;
    while (run_) {
      MarshalledResponse response;
      bool msg_received = zmq_recv_socket_->recv(&response.msg, ZMQ_NOBLOCK);
      // basically implements polling to avoid blocking recv calls
      // not ideal but whatever, this is a research project! :-D
      if (!msg_received) {
        continue;
      }
      total_received++;

      //cout << "response set size "  << response.Set().NumClusters() << " clusters\n";

      response_queue_->push(std::move(response));
    }

    cout << "Work queue thread ending. Total received: " << total_received
         << "\n";
  });

  // partial mergers in this thread may need to be more fully parallelized to
  // prevent bottlenecks. May be required to use a different structure for
  // tracking partial mergers rather than the current map, which needs to be
  // locked
  std::ofstream return_times("return_times.txt");

  auto worker_func = [this, &return_times]() {
    // read from result queue
    // if is a batch result and is small enough, push to WorkManager
    // if is partial result (ID will be
    // in merge map), lookup and merge with partial_set, if partial set
    // complete,
    //    push to WorkManager
    MarshalledResponse response;
    while (run_) {
      if (!response_queue_->pop(response)) {
        continue;
      }

      auto id = response.ID();
      // cout << "repsonse id is " << id << "\n";
      PartialMergeItem* partial_item;

      {
        absl::MutexLock l(&mu_);
        auto partial_it = partial_merge_map_.find(id);
        if (partial_it == partial_merge_map_.end()) {
          if (id != -1) {
            cout << "error thing in map was not -1, was " << id << "\n";
            exit(0);
          }
          //cout << "pushing full result \n";
          /*if (response.set().clusters_size() > params.dup_removal_thresh) {
            RemoveDuplicates(*response.mutable_set());
          }*/
          MarshalledClusterSet new_set(response);

          sets_to_merge_queue_->push(std::move(new_set));
          continue;
        } else {
          // do assign
          partial_item = &partial_it->second;
        }
        if (outstanding_merges_ == 1) {
          return_times << static_cast<long int>(std::time(0)) << "\n";
        }
      }

      /*MergePartials(partial_item->partial_set, response.set(),
                    partial_item->original_size);*/

      partial_item->partial_set.MergeClusterSet(response.Set());
      // cout << "done\n";
      // check if this was the last one
      auto val = partial_item->num_received++;
      if (partial_item->num_expected - 1 == val) {
        // TODO are we sure that all other potential threads are finished
        // merging their partials with this one?

        // go through  and delete any fully merged clusters
        // do this when constructing the new full set
        // partial_item->partial_set.RemoveFullyMerged();

        MarshalledClusterSet set;
        partial_item->partial_set.BuildMarshalledSet(&set);
        /*if (set.clusters_size() >
            params.dup_removal_thresh) {
          RemoveDuplicates(set);
        }*/
        // remove partial it, its done now
        //cout << "partial id " << id << " is complete, with " << set.NumClusters() << " clusters\n";
        {
          if (outstanding_merges_ == 1) {
            cout << "last request complete, 1 merge left, time: " << static_cast<long int>(std::time(0)) << "\n";
          }
          absl::MutexLock l(&mu_);
          partial_merge_map_.erase(id);
        }
        sets_to_merge_queue_->push(std::move(set));
        // set partial not outstanding
        // outstanding_partial_ = false;
      }
    }
  };

  worker_threads_.reserve(params.num_threads);
  for (int i = 0; i < params.num_threads; i++) {
    worker_threads_.push_back(std::thread(worker_func));
  }

  cout << "loading to marshalled sets\n";
  // dump all sequences in single cluster sets into the queue
  // for now assumes all of this will fit in memory
  // even a million sequences would just be a few MB
  int i = 0;
  for (const auto& s : sequences_) {
    /*cmproto::ClusterSet set;
    auto* c = set.add_clusters();
    c->add_indexes(s.ID());*/
    if (params.dataset_limit > 0) {
      if (i == params.dataset_limit) break;
    }
    MarshalledClusterSet set(s.ID());
    sets_to_merge_queue_->push(std::move(set));
    i++;
  }

  cout << "done\n";

  auto t0 = std::chrono::high_resolution_clock::now();
  // just use 'this' thread to schedule outgoing work
  // take stuff from sets_to_merge_ and schedule the work in
  // batches or split partial merges for larger sets
  std::vector<MarshalledClusterSet> sets;
  sets.resize(2);
  while (outstanding_merges_ > 0) {
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
    if (sets[0].NumClusters() < params.batch_size ||
        sets[1].NumClusters() < params.batch_size) {
      // create a batch, they are small
      //cout << "two sets are small, batching ...\n";

      uint32_t total_clusters = sets[0].NumClusters() + sets[1].NumClusters();
      //cout << "total clusters: " << total_clusters << "\n";
      // marshal id, type, num sets
      // int id = -1; auto t = RequestType::Batch; int num_sets = 2;
      // sets[1].MarshalToBuffer(request.buf)
      MarshalledRequest request;
      request.CreateBatchRequest(-1);
      request.AddSetToBatch(sets[0]);
      request.AddSetToBatch(sets[1]);

      outstanding_merges_--;
      while (total_clusters < params.batch_size && outstanding_merges_ > 0) {
        // add more cluster sets to batch
        if (!sets_to_merge_queue_->pop(sets[0])) {
          return agd::errors::Internal(
              "ERROR: did not get set for second to merge with.");
        }

        // c = batch->add_sets();
        // c->CopyFrom(sets[0]);
        total_clusters += sets[0].NumClusters();
        request.AddSetToBatch(sets[0]);
        outstanding_merges_--;
      }
      // cout << "batched " << batch->sets_size() << " sets\n";
      // if the queue uses copy semantics im not sure how protobufs
      // with submessages will behave
      request_queue_->push(std::move(request));

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

      //cout << "pushing partial ...\n";
      /*cout << "waiting for outstanding partial\n";
      while(outstanding_partial_) {
        ;;
      }
      cout << "outstanding finished\n";*/
      PartialMergeItem item;
      item.num_received = 0;
      // use the outstanding merges as id
      if (sets[0].NumClusters() < sets[1].NumClusters()) {
        // sets[0].Swap(&sets[1]);
        std::swap(sets[0], sets[1]);
      }
      // each work item does a partial merge of one cluster in sets[0] to
      // all clusters in sets[1]
      item.num_expected = sets[0].NumClusters();
      // item.partial_set.CopyFrom(sets[1]);
      item.partial_set.Init(sets[1]);
      // item.original_size = sets[1].clusters_size();
      {
        absl::MutexLock l(&mu_);
        //cout << "pushing id " << outstanding_merges_ << " to map\n";
        partial_merge_map_.insert_or_assign(outstanding_merges_, std::move(item));
      }

      MarshalledClusterView cluster;
      //cout << "pushing id " << outstanding_merges_ << "\n";
      uint32_t total_cluster = sets[0].NumClusters();
      uint32_t i = 0;
      while (sets[0].NextCluster(&cluster)) {
        i++;
        // for (const auto& c : sets[0].clusters()) {
        MarshalledRequest request;
        request.CreatePartialRequest(outstanding_merges_, cluster, sets[1]);
        // cout << "pushing partial request with " <<
        // partial_request->set().clusters_size() << " clusters in set and ID: "
        // << request.id() << "\n";
        request_queue_->push(std::move(request));
      }
      if (outstanding_merges_ == 1) {
        cout << "last request sent, 1 merge left, time: " << static_cast<long int>(std::time(0)) << "\n";
      }
      assert(i == total_cluster);
      // set partial outstanding
      // outstanding_partial_ = true;
    }
    //cout << "outstanding merges: " << outstanding_merges_ << "\n";
  }

  // we must now wait for the last results to come in
  // wait for worker thread to push last merged set
  // TODO add a timeout or something?
  cout << "done and waiting for final result...\n";
  // while (sets_to_merge_queue_->size() != 1);;

  cout << "scheduling final alignments on controller...\n";
  MarshalledClusterSet final_set;
  sets_to_merge_queue_->peek(final_set);
  auto t1 = std::chrono::high_resolution_clock::now();
  cout << "final set size is " << final_set.NumClusters() << " clusters\n";
  cout << "partial merge map size is " << partial_merge_map_.size() << "\n";
  cout << "sets to merge size is " << sets_to_merge_queue_->size() << "\n";
  cout << "request queue size " << request_queue_->size() << "\n";
  cout << "response queue size " << response_queue_->size() << "\n";

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time: " << sec.count() << " seconds.\n";

  std::ofstream timing_file("dist_timing.txt", std::ofstream::out);
  timing_file << sec.count() << "\n";

  ClusterSet set(final_set, sequences_);
  set.DumpJson("dist_clusters.json");

  if (!params.exclude_allall) {
    AllAllExecutor executor(std::thread::hardware_concurrency(), 500, &envs,
        &aligner_params);
    executor.Initialize();
    set.ScheduleAlignments(&executor);
    executor.FinishAndOutput("dist_output_dir");
  } else {
    cout << "Skipping all-all alignments ...\n";
  }

  cout << "clustering complete!! Joining threads ...\n";

  run_ = false;
  response_queue_->unblock();
  request_queue_->unblock();
  sets_to_merge_queue_->unblock();
  // worker_thread_.join();
  for (auto& t : worker_threads_) {
    t.join();
  }
  request_queue_thread_.join();
  response_queue_thread_.join();

  cout << "All threads joined.\n";

  return agd::Status::OK();
}
