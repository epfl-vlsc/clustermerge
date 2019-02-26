
#include "absl/container/flat_hash_set.h"
#include "absl/synchronization/mutex.h"
#include "src/comms/requests.h"
#include <iostream>

class IndexedCluster {
 public:
  IndexedCluster() = default;
  IndexedCluster(const IndexedCluster& other) {
    seq_indexes_ = other.seq_indexes_;
    fully_merged_ = other.fully_merged_;
    respresentative_idx_ = other.respresentative_idx_;
    orig_seqs_ = other.orig_seqs_;
  }
  IndexedCluster& operator=(const IndexedCluster& other) {
    seq_indexes_ = other.seq_indexes_;
    fully_merged_ = other.fully_merged_;
    respresentative_idx_ = other.respresentative_idx_;
    orig_seqs_ = other.orig_seqs_;
    return *this;
  }
  IndexedCluster(const MarshalledClusterView& c) {
    respresentative_idx_ = c.SeqIndex(0);
    uint32_t num_seqs = c.NumSeqs();
    for (uint32_t i = 0; i < num_seqs; i++) {
      seq_indexes_.insert(c.SeqIndex(i));
    }
    orig_seqs_ = num_seqs;
    if (num_seqs != seq_indexes_.size()) {
      std::cout <<  "there were dups!!!!!!!!!!!!!\n";
    }
  }
  void Insert(uint32_t seq_index) {
    absl::MutexLock l(&mu_);
    seq_indexes_.insert(seq_index);
  }
  void SetFullyMerged() { fully_merged_ = true; }
  bool IsFullyMerged() const { return fully_merged_; }

  const absl::flat_hash_set<uint32_t>& SeqIndexes() const { return seq_indexes_; }
  uint32_t Representative() const { return respresentative_idx_; }
  uint32_t NumOrigSeqs() const { return orig_seqs_; }

 private:
  absl::flat_hash_set<uint32_t> seq_indexes_;
  bool fully_merged_ = false;
  // track explicitly the rep because the set is not ordered
  uint32_t respresentative_idx_ = 0;
  uint32_t orig_seqs_;
  // lock to insert new set indexes
  absl::Mutex mu_;
};

class PartialMergeSet {
 public:
  PartialMergeSet& operator=(PartialMergeSet&& other) {
    clusters_ = std::move(other.clusters_);
    new_clusters_ = std::move(other.new_clusters_);
    return *this;
  };
  void MergeClusterSet(MarshalledClusterSetView set);
  // build final set, not including fully merged clusters
  void BuildMarshalledSet(MarshalledClusterSet* set);
  void Init(MarshalledClusterSet& set);
  //void RemoveFullyMerged();
  const std::vector<IndexedCluster>& Clusters() const { return clusters_; }

 private:
  std::vector<IndexedCluster> clusters_;  // does not change after construction
  std::vector<MarshalledCluster> new_clusters_;
  // lock to add new clusters
  absl::Mutex mu_;
};
