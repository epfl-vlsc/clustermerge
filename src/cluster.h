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
    residue_total_ = other.residue_total_;
    longest_ = other.longest_;
    duplicate_ = other.duplicate_;
  }

  Cluster& operator=(Cluster&& other) noexcept {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    residue_total_ = other.residue_total_;
    longest_ = other.longest_;
    duplicate_ = other.duplicate_;
    return *this;
  }

  void Merge(Cluster* other, ProteinAligner* aligner);
  void MergeOther(Cluster* other, ProteinAligner* aligner);

  agd::Status AlignReps(const Cluster& other,
                        ProteinAligner::Alignment* alignment,
                        ProteinAligner* aligner);

  bool PassesThreshold(const Cluster& other, ProteinAligner* aligner);

  const Sequence& Rep() { return seqs_.front(); }

  // add (move) seq into seqs_
  void AddSequence(const Sequence& seq);

  // mark that this cluster has been fully merged
  // with another, and will go away. seqs_ may not be valid anymore
  void SetFullyMerged() { fully_merged_ = true; }
  bool IsFullyMerged() const { return fully_merged_; }
  
  void SetDuplicate() { duplicate_ = true; }
  bool IsDuplicate() const { return duplicate_; }

  const std::list<Sequence>& Sequences() const { return seqs_; }

  uint64_t Residues() const { return residue_total_; }
  
  uint64_t LongestLength() const { return longest_; }

  void Lock() { mu_.Lock(); }
  void Unlock() { mu_.Unlock(); }

 private:
  // representative is first seq
  // use a list so refs aren't invalidated
  std::list<Sequence> seqs_;
  bool fully_merged_ = false;
  uint64_t residue_total_ = 0;
  uint64_t longest_ = 0;
  bool duplicate_ = false;

  absl::Mutex mu_; // protects seqs_ and fully_merged_

};
