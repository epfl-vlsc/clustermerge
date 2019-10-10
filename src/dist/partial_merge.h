
#include <iostream>
#include "absl/container/flat_hash_set.h"
#include "absl/synchronization/mutex.h"
#include "src/comms/requests.h"

/*template <typename T>
struct Node {
  T t;
  Node* next;
};

template <typename T>
class AtomicList {
  std::atomic<Node<T>*> head_{nullptr};
  std::atomic_uint_fast32_t size_{0};

 public:
  // move ops are NOT threadsafe
  AtomicList() = default;
  AtomicList(AtomicList&& other) {
    head_ = other.head_;
    other.head_ = nullptr;
    size_ = other.size_.load();
  }
  AtomicList& operator=(AtomicList&& other) {
    head_ = other.head_.load();
    other.head_ = nullptr;
    size_ = other.size_.load();
    return *this;
  }
  ~AtomicList() {
    auto* current = head_.load();
    Node<T>* next;
    while (current != nullptr) {
      next = current->next;
      delete current;
      current = next;
    }
  }
  void push_atomic(T t) {
    auto p = new Node<T>();
    p->t = t;
    p->next = head_;
    while (!head_.compare_exchange_weak(p->next, p)) {
    }
    size_++;
  }
  uint32_t Size() const { return size_.load(); }
  // this is bad practice but im in a hurry and too lazy to implement iterators
  Node<T>* Head() const { return head_.load(); }
};*/

class IndexedCluster {
 public:
  IndexedCluster() = default;
  IndexedCluster(IndexedCluster&& other) {
    seq_indexes_ = std::move(other.seq_indexes_);
    fully_merged_ = other.fully_merged_;
    respresentative_idx_ = other.respresentative_idx_;
    new_seqs_hash_ = std::move(other.new_seqs_hash_);
    orig_seqs_ = other.orig_seqs_;
  }
  IndexedCluster& operator=(IndexedCluster&& other) {
    seq_indexes_ = std::move(other.seq_indexes_);
    fully_merged_ = other.fully_merged_;
    respresentative_idx_ = other.respresentative_idx_;
    new_seqs_hash_ = std::move(other.new_seqs_hash_);
    orig_seqs_ = other.orig_seqs_;
    return *this;
  }
  IndexedCluster(const MarshalledClusterView& c) {
    respresentative_idx_ = c.SeqIndex(0);
    uint32_t num_seqs = c.NumSeqs();
    seq_indexes_.reserve(num_seqs);
    for (uint32_t i = 0; i < num_seqs; i++) {
      seq_indexes_.push_back(c.SeqIndex(i));
    }
    orig_seqs_ = num_seqs;
    if (num_seqs != seq_indexes_.size()) {
      std::cout << "there were dups!!!!!!!!!!!!!\n";
    }
  }
  void Insert(uint32_t seq_index) {
    absl::MutexLock l(&mu_);
    new_seqs_hash_.insert(seq_index);
    // seq_indexes_.insert(seq_index);
    //new_seqs_.push_atomic(seq_index);
  }
  void SetFullyMerged() { fully_merged_ = true; }
  bool IsFullyMerged() const { return fully_merged_; }

  const std::vector<uint32_t>& SeqIndexes() const { return seq_indexes_; }
  uint32_t Representative() const { return respresentative_idx_; }
  uint32_t NumOrigSeqs() const { return orig_seqs_; }
  uint32_t NumNewSeqs() const { return new_seqs_hash_.size(); }

  // add new seqs to set to facilitate duplicate removal
  void AddNewSeqs(absl::flat_hash_set<uint32_t>* set) const {
    for (const auto& x : new_seqs_hash_) {
      set->insert(x);
    }
  }

 private:
  std::vector<uint32_t> seq_indexes_;
  //AtomicList<uint32_t> new_seqs_;
  absl::flat_hash_set<uint32_t> new_seqs_hash_;
  bool fully_merged_ = false;
  uint32_t orig_seqs_;
  // track explicitly the rep because the set is not ordered
  uint32_t respresentative_idx_ = 0;
  // lock to insert new set indexes
  absl::Mutex mu_;
};

class PartialMergeSet {
 public:
  PartialMergeSet& operator=(PartialMergeSet&& other) {
    clusters_set1_ = std::move(other.clusters_set1_);
    clusters_set2_ = std::move(other.clusters_set2_);
    return *this;
  };
  void MergeClusterSet(MarshalledClusterSetView set, int start_index,
                       int end_index, int cluster_index);
  // build final set, not including fully merged clusters
  void BuildMarshalledSet(MarshalledClusterSet* set);
  void Init(MarshalledClusterSet& set1, MarshalledClusterSet& set2);
  // void RemoveFullyMerged();
  const std::vector<IndexedCluster>& Clusters() const { return clusters_set2_; }

 private:
  std::vector<IndexedCluster> clusters_set1_;
  std::vector<IndexedCluster>
      clusters_set2_;  // does not change after construction
};
