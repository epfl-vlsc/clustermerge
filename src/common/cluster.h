#pragma once

#include <list>
#include "src/agd/errors.h"
#include "aligner.h"
#include "sequence.h"
#include "absl/synchronization/mutex.h"
#include "src/comms/requests.h"

class Cluster {
 public:
  Cluster() = default;  // an empty cluster
  Cluster(const Sequence& seed) { seqs_.push_back(seed); }

  Cluster(Cluster&& other) noexcept {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    duplicate_ = other.duplicate_;
  }

  Cluster& operator=(Cluster&& other) noexcept {
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    duplicate_ = other.duplicate_;
    return *this;
  }

  Cluster(const MarshalledClusterView& cluster, const std::vector<Sequence>& sequences) {
    // construct a cluster object from a protobuf representation
    uint32_t num_seqs = cluster.NumSeqs();
    for (size_t seq_i = 0; seq_i < num_seqs; seq_i++) {
      AddSequence(sequences[cluster.SeqIndex(seq_i)]);
    }
    fully_merged_ = cluster.IsFullyMerged();
  }


  void Merge(Cluster* other, ProteinAligner* aligner);
  void MergeOther(Cluster* other, ProteinAligner* aligner);

  agd::Status AlignReps(const Cluster& other,
                        ProteinAligner::Alignment* alignment,
                        ProteinAligner* aligner);

  bool PassesThreshold(const Cluster& other, ProteinAligner* aligner);

  const Sequence& Rep() { return seqs_.front(); }

  // add seq into seqs_
  void AddSequence(const Sequence& seq);

  // mark that this cluster has been fully merged
  // with another, and will go away. seqs_ may not be valid anymore
  void SetFullyMerged() { fully_merged_ = true; }
  bool IsFullyMerged() const { return fully_merged_; }
  
  void SetDuplicate() { duplicate_ = true; }
  bool IsDuplicate() const { return duplicate_; }

  const std::list<Sequence>& Sequences() const { return seqs_; }

  void Lock() { mu_.Lock(); }
  void Unlock() { mu_.Unlock(); }

  void MarshalToBuffer(agd::Buffer* buf);
  // marshalled cluster is [fully_merged, num_idx, (cluster indexes)]
  uint32_t ByteSize() { return sizeof(bool) + sizeof(int) + sizeof(int)*seqs_.size(); }

 private:
  // representative is first seq
  // use a list so refs aren't invalidated
  std::list<Sequence> seqs_;
  bool fully_merged_ = false;
  bool duplicate_ = false;

  absl::Mutex mu_; // protects seqs_ and fully_merged_

};
