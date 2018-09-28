
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include "absl/synchronization/mutex.h"
#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"

/*namespace {
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}
}

struct PairHash {
  template <class T1, class T2>
  size_t operator()(const std::pair<T1, T2>& p) const {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T1>{}(p.second);
    hash_combine(h1, h2);
    return h1;
  }
};*/

typedef std::pair<absl::string_view, absl::string_view>
    GenomePair;  // could be replaced by ints that map to genome strings
typedef std::pair<int, int> SequencePair;

class CandidateMap {
  typedef absl::flat_hash_map<
      GenomePair, absl::flat_hash_map<SequencePair, bool>>
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
