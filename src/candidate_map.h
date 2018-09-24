
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include "absl/synchronization/mutex.h"

namespace {
template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}  // namespace

struct PairHash {
  template <class T1, class T2>
  size_t operator()(const std::pair<T1, T2>& p) const {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T1>{}(p.second);
    hash_combine(h1, h2);
    return h1;
  }
};

typedef std::pair<std::string, std::string>
    GenomePair;  // could be replaced by ints that map to genome strings
typedef std::pair<int, int> SequencePair;

class CandidateMap {
  typedef std::unordered_map<
      GenomePair, std::unordered_map<SequencePair, bool, PairHash>, PairHash>
      GenomeSequenceMap;

 public:
  // returns true, or false and inserts key
  bool ExistsOrInsert(const GenomePair& g, const SequencePair& s) {
    absl::MutexLock l(&mu_);

    auto genome_pair_it = map_.find(g);
    if (genome_pair_it != map_.end()) {
      auto seq_pair_it = genome_pair_it->second.find(s);
      if (seq_pair_it != genome_pair_it->second.end()) {
        return true;
      }
    }

    map_[g][s] = true;
    return false;
  }

 private:
  GenomeSequenceMap map_;
  absl::Mutex mu_;
};

// class for associating sequence ID pairs
// as passed threshold or not
class SequenceIDMap {
 public:
  // determine if a sequencepair has been aligned, and
  // if it passed the threshold
  bool Exists(const SequencePair& s, bool* passed_threshold)
      LOCKS_EXCLUDED(mu_) {
    absl::MutexLock l(&mu_);

    auto it = map_.find(s);
    if (it != map_.end()) {
      *passed_threshold = it->second;
      return true;
    } else {
      return false;
    }
  }

  void Insert(const SequencePair& s, const bool& passed_threshold)
      LOCKS_EXCLUDED(mu_) {
    absl::MutexLock l(&mu_);

    map_[s] = passed_threshold;
  }

 private:
  absl::Mutex mu_;
  std::unordered_map<SequencePair, bool, PairHash> map_ GUARDED_BY(mu_);
};