#pragma once

#include <vector>
#include "agd/errors.h"
#include "sequence.h"

class Cluster {
 public:
  Cluster(Sequence& seed) { seqs_.push_back(seed); }

  Cluster(Cluster&& other) {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
  }

  agd::Status CompareCluster(const Cluster& other);

  bool IsFullyMerged() const { return fully_merged_; }

  // mark that this cluster has been fully merged
  // with another, and will go away
  void SetFullyMerged() { fully_merged_ = true; }

  const std::vector<Sequence>& Sequences() const { return seqs_; }

 private:
  // representative is first seq
  std::vector<Sequence> seqs_;
  bool fully_merged_ = false;
};