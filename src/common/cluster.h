#pragma once

#include <list>
#include "src/agd/errors.h"
#include "aligner.h"
#include "sequence.h"
#include "absl/synchronization/mutex.h"
#include "src/comms/requests.h"

class Cluster {
 public:
  Cluster(const std::vector<Sequence>& sequences) : all_seqs_(&sequences) {};  // an empty cluster
  Cluster(uint32_t seed, const std::vector<Sequence>& sequences) : all_seqs_(&sequences) { seqs_.push_back(seed); }

  Cluster(Cluster&& other) noexcept {
    all_seqs_ = other.all_seqs_;
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    duplicate_ = other.duplicate_;
  }

  Cluster& operator=(Cluster&& other) noexcept {
    all_seqs_ = other.all_seqs_;
    seqs_ = std::move(other.seqs_);
    fully_merged_ = other.fully_merged_;
    duplicate_ = other.duplicate_;
    return *this;
  }

  Cluster(const MarshalledClusterView& cluster, const std::vector<Sequence>& sequences) : all_seqs_(&sequences) {
    // construct a cluster object from a protobuf representation
    uint32_t num_seqs = cluster.NumSeqs();
    seqs_.reserve(num_seqs);
    for (size_t seq_i = 0; seq_i < num_seqs; seq_i++) {
      AddSequence(seq_i);
    }
    fully_merged_ = cluster.IsFullyMerged();
  }


  void Merge(Cluster* other, ProteinAligner* aligner);
  void MergeOther(Cluster* other, ProteinAligner* aligner);

  agd::Status AlignReps(const Cluster& other,
                        ProteinAligner::Alignment* alignment,
                        ProteinAligner* aligner);

  bool PassesThreshold(const Cluster& other, ProteinAligner* aligner);

  uint32_t Rep() { return seqs_.front(); }
  const Sequence& SeqRep() const { return all_seqs_->at(seqs_.front()); }
  const Sequence& SeqAt(uint32_t idx) const { return all_seqs_->at(idx); }

  // add seq into seqs_
  void AddSequence(uint32_t seq);

  // mark that this cluster has been fully merged
  // with another, and will go away. seqs_ may not be valid anymore
  void SetFullyMerged() { fully_merged_ = true; }
  bool IsFullyMerged() const { return fully_merged_; }
  
  void SetDuplicate() { duplicate_ = true; }
  bool IsDuplicate() const { return duplicate_; }

  void Reserve(size_t num_seqs) {
    seqs_.reserve(num_seqs);
  }

  const std::vector<uint32_t>& Sequences() const { return seqs_; }
  const std::vector<Sequence>& AllSequences() const { return *all_seqs_; }

  void Lock() { mu_.Lock(); }
  void Unlock() { mu_.Unlock(); }

  // marshalled cluster is [fully_merged, num_idx, (cluster indexes)]
  uint32_t ByteSize() { return sizeof(bool) + sizeof(int) + sizeof(int)*seqs_.size(); }

 private:
  // representative is first seq
  // use a list so refs aren't invalidated
  // NOTE testing vector here for dist version mem consumption
  // TODO have a base cluster class and inherit versions for dist and local?
  std::vector<uint32_t> seqs_;
  const std::vector<Sequence>* all_seqs_ = nullptr;
  bool fully_merged_ = false;
  bool duplicate_ = false;

  absl::Mutex mu_; // protects seqs_ and fully_merged_

};
