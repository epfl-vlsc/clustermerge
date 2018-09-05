
#pragma once

#include "cluster.h"
#include <vector> 

class ClusterSet {
 public:
  ClusterSet(Sequence& seed) {
    Cluster c(seed);
    clusters_.push_back(std::move(c));
  }

  ClusterSet(size_t num) {
    clusters_.resize(num);
  }

  // merge two cluster sets by building a new one
  ClusterSet MergeClusters(ClusterSet& other);

  void DebugDump() const;

 private:
  ClusterSet() = default;
  std::vector<Cluster> clusters_;
};