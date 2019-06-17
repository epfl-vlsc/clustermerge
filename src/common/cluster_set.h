
#pragma once

#include <vector>
#include "all_all_executor.h"
#include "cluster.h"
#include "src/comms/requests.h"

class MergeExecutor;

class ClusterSet {
 public:
  ClusterSet() = default;
  ClusterSet(ClusterSet&& other) { clusters_ = std::move(other.clusters_); }
  ClusterSet(Sequence& seed) {
    // construct from a single sequence
    Cluster c(seed);
    clusters_.push_back(std::move(c));
  }
  ClusterSet(Cluster& c) {
    // construct from a single sequence
    clusters_.push_back(std::move(c));
  }

  ClusterSet(size_t num) { clusters_.reserve(num); }

  // construct from protobuf (for dist version)
  ClusterSet(MarshalledClusterSet& marshalled_set,
             const std::vector<Sequence>& sequences);

  ClusterSet(MarshalledClusterSetView& marshalled_set,
             const std::vector<Sequence>& sequences);

  void BuildMarshalledResponse(int id, MarshalledResponse* response);

  // void ConstructProto(cmproto::ClusterSet* set_proto);

  void Swap(ClusterSet* other) { clusters_.swap(other->clusters_); }

  // merge two cluster sets by building a new one
  ClusterSet MergeClusters(ClusterSet& other, ProteinAligner* aligner);

  // execute a PartialMerge for the distributed runtime
  // merge cluster into cluster set
  // keep any fully merged cluster, but mark it
  // the central controller will eliminate fully merged clusters
  // after merging partial results
  ClusterSet MergeCluster(Cluster& c_other, ProteinAligner* aligner,
                          bool& worker_signal_);

  // merge two cluster sets by building a new one, in parallel (uses std::async)
  ClusterSet MergeClustersParallel(ClusterSet& other, MergeExecutor* executor);

  //Add by akash
  void AddCluster(Cluster& c) { clusters_.push_back(std::move(c)); }

  // merge `this` with `cluster`, called from parallel merge executor
  void MergeClusterLocked(Cluster* cluster, ProteinAligner* aligner);

  // schedule all-all alignments onto the executor threadpool
  void ScheduleAlignments(AllAllExecutor* executor);

  // remove duplicate clusters from this set
  void RemoveDuplicates();

  void DebugDump() const;
  void DumpJson(const std::string& filename, std::vector <std::string>& dataset_file_names) const;

  void MarshalToBuffer(agd::Buffer* buf);

  size_t Size() { return clusters_.size(); }

 private:
  std::vector<Cluster> clusters_;
};
