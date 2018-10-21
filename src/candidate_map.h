
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include "absl/synchronization/mutex.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"


typedef std::pair<absl::string_view, absl::string_view>
    GenomePair;  // could be replaced by ints that map to genome strings
typedef std::pair<int, int> SequencePair;
typedef std::pair<uint32_t, uint32_t> AbsSequencePair;

class CandidateMap {
  typedef absl::flat_hash_set<AbsSequencePair> AbsSequenceMap;

 public:
  CandidateMap() = default;
  CandidateMap(size_t size) {
    map_.reserve(size);
  }
  // do not copy or move
  CandidateMap(CandidateMap&& other) = delete;
  CandidateMap(const CandidateMap& other) = delete;

  // returns true, or false and inserts key
  bool ExistsOrInsert(const AbsSequencePair& s) {
    if (map_.contains(s)) {
      return true;
    } else {
      map_.insert(s);
      return false;
    }
  }

 private:
  AbsSequenceMap map_;
  absl::Mutex mu_;
};
