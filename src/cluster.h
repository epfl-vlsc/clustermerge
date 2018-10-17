#pragma once

#include <list>
#include "agd/errors.h"
#include "aligner.h"
#include "sequence.h"
#include "absl/synchronization/mutex.h"

class Cluster {
 public:
  Cluster() = default;  // an empty cluster
  Cluster(Sequence& seed) { seqs_.push_back(seed); }

  Cluster(Cluster&& other) noexcept {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    residue_total = other.residue_total;
  }

  Cluster& operator=(Cluster&& other) {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    residue_total = other.residue_total;
    return *this;
  }

  void Merge(Cluster* other, ProteinAligner* aligner);
  void MergeOther(Cluster* other, ProteinAligner* aligner);

  agd::Status AlignReps(const Cluster& other,
                        ProteinAligner::Alignment* alignment,
                        ProteinAligner* aligner);

  bool PassesThreshold(const Cluster& other, ProteinAligner* aligner);

  bool IsFullyMerged() const { return fully_merged_; }

  const Sequence& Rep() { return seqs_.front(); }

  // add (move) seq into seqs_
  void AddSequence(const Sequence& seq);

  // mark that this cluster has been fully merged
  // with another, and will go away. seqs_ may not be valid anymore
  void SetFullyMerged() { fully_merged_ = true; }

  const std::list<Sequence>& Sequences() const { return seqs_; }

  uint64_t Residues() const { return residue_total; }

  void Lock() { mu_.Lock(); }
  void Unlock() { mu_.Unlock(); }

 private:
  // representative is first seq
  // use a list so refs aren't invalidated
  std::list<Sequence> seqs_;
  bool fully_merged_ = false;
  uint64_t residue_total = 0;

  absl::Mutex mu_; // protects seqs_ and fully_merged_

};
