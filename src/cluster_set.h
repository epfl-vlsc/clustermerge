
#pragma once

#include <vector>
#include "all_all_executor.h"
#include "cluster.h"

class ClusterSet {
 public:
  ClusterSet(Sequence& seed) {
    Cluster c(seed);
    clusters_.push_back(std::move(c));
  }

  ClusterSet(size_t num) { clusters_.reserve(num); }

  // merge two cluster sets by building a new one
  ClusterSet MergeClusters(ClusterSet& other, ProteinAligner* aligner);

  // schedule all-all alignments onto the executor threadpool
  void ScheduleAlignments(AllAllExecutor* executor);

  void DebugDump() const;

 private:
  ClusterSet() = default;
  std::vector<Cluster> clusters_;
};