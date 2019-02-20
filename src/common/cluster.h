#pragma once

#include <list>
#include "src/agd/errors.h"
#include "aligner.h"
#include "sequence.h"
#include "absl/synchronization/mutex.h"
#include "src/proto/cluster.pb.h"

class Cluster {
 public:
  Cluster() = default;  // an empty cluster
  Cluster(const Sequence& seed) { seqs_.push_back(seed); }

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

  Cluster(const cmproto::Cluster& cluster_proto, const std::vector<Sequence>& sequences) {
    // construct a cluster object from a protobuf representation
    for (size_t seq_i = 0; seq_i < cluster_proto.indexes_size(); seq_i++) {
      AddSequence(sequences[cluster_proto.indexes(seq_i)]);
    }
    fully_merged_ = cluster_proto.fully_merged();
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

  uint64_t Residues() const { return residue_total_; }
  
  uint64_t LongestLength() const { return longest_; }

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
  uint64_t residue_total_ = 0;
  uint64_t longest_ = 0;
  bool duplicate_ = false;

  absl::Mutex mu_; // protects seqs_ and fully_merged_

};
