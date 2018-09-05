
#include "cluster_set.h"
#include "aligner.h"
#include <iostream>
#include "debug.h"

ClusterSet ClusterSet::MergeClusters(ClusterSet& other) {
  // merge clusters, clusters can "disappear" from either
  // set, so we just create a new one and resize its internal
  // cluster vector for a single alloc

  ClusterSet new_cluster_set(clusters_.size() + other.clusters_.size());

  ProteinAligner::Alignment alignment;
  for (auto& c : clusters_) {
    for (auto& c_other : other.clusters_) {
      // s = c.AlignReps(c_other, alignment);
      // analyze 4 different merge cases
    }
  }
}

void ClusterSet::DebugDump() const {
  std::cout << "Dumping " << clusters_.size() << " clusters in set... \n";
  for (const auto& cluster : clusters_) {
    std::cout << "\tCluster seqs:\n";
    for (const auto& seq : cluster.Sequences()) {
      std::cout << "\t\tGenome: " << seq.Genome() << ", sequence: "
                << PrintNormalizedProtein(seq.Seq().data(), seq.Seq().length())
                << "\n\n";
    }
  }
}