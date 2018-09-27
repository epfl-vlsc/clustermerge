
#pragma once

#include <vector>
#include "all_all_executor.h"
#include "cluster.h"

class MergeExecutor;

class ClusterSet {
 public:
  ClusterSet() = default;
  ClusterSet(ClusterSet&& other) {
    clusters_ = std::move(other.clusters_);
  }
  ClusterSet(Sequence& seed) {
    Cluster c(seed);
    clusters_.push_back(std::move(c));
  }

  ClusterSet(size_t num) { clusters_.reserve(num); }

  // merge two cluster sets by building a new one
  ClusterSet MergeClusters(ClusterSet& other, ProteinAligner* aligner);
  
  // merge two cluster sets by building a new one, in parallel (uses std::async)
  ClusterSet MergeClustersParallel(ClusterSet& other, MergeExecutor* executor);

  // merge `this` with `cluster`, called from parallel merge executor
  void MergeClusterLocked(Cluster* cluster, ProteinAligner* aligner);

  // schedule all-all alignments onto the executor threadpool
  void ScheduleAlignments(AllAllExecutor* executor);

  void DebugDump() const;
  void DumpJson() const;

  size_t Size() { return clusters_.size(); }

 private:
  std::vector<Cluster> clusters_;
};
