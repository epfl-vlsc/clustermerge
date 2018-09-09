#pragma once

#include <vector>
#include "agd/errors.h"
#include "aligner.h"
#include "sequence.h"

class Cluster {
 public:
  Cluster() = default;  // an empty cluster
  Cluster(Sequence& seed) { seqs_.push_back(seed); }

  Cluster(Cluster&& other) {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
  }

  void Merge(const Cluster& other, ProteinAligner* aligner);

  agd::Status AlignReps(const Cluster& other,
                        ProteinAligner::Alignment* alignment,
                        ProteinAligner* aligner);

  bool PassesThreshold(const Cluster& other, ProteinAligner* aligner);

  bool IsFullyMerged() const { return fully_merged_; }

  const Sequence& Rep() { return seqs_[0]; }

  // add (move) seq into seqs_
  void AddSequence(const Sequence& seq);

  // mark that this cluster has been fully merged
  // with another, and will go away. seqs_ may not be valid anymore
  void SetFullyMerged() { fully_merged_ = true; }

  const std::vector<Sequence>& Sequences() const { return seqs_; }

 private:
  // representative is first seq
  std::vector<Sequence> seqs_;
  bool fully_merged_ = false;
};