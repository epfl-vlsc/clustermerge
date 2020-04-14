
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <chrono>
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"


typedef std::pair<absl::string_view, absl::string_view>
    GenomePair;
typedef std::pair<int, int> SequencePair;
typedef std::pair<uint32_t, uint32_t> AbsSequencePair;

class CandidateMap {
  typedef absl::flat_hash_set<AbsSequencePair> AbsSequenceMap;

 public:
  CandidateMap() = default;
  CandidateMap(size_t size) {
    map_.reserve(size);
  }
  ~CandidateMap() {
    //std::cout << "Candidate map was rehashed " << times_rehashed << " times.\n";
    //std::cout << "The longest insert time was: " << longest_ << " ms.\n";
  }
  // do not copy or move
  CandidateMap(CandidateMap&& other) = delete;
  CandidateMap(const CandidateMap& other) = delete;
  CandidateMap& operator=(const CandidateMap& other) = delete;

  // returns true, or false and inserts key
  bool ExistsOrInsert(const AbsSequencePair& s) {
    if (map_.contains(s)) {
      return true;
    } else {
      auto start = std::chrono::steady_clock::now();

      auto old_load_factor = map_.load_factor();
      map_.insert(s);
      auto new_load_factor = map_.load_factor();
      if (new_load_factor < old_load_factor) {
        times_rehashed++;
      }

      auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - start)
                     .count();
      if (diff > longest_) {
        longest_ = diff;
      }
      return false;
    }
  }

  size_t size(){
    return map_.size();
  }

 private:
  int times_rehashed = 0;
  int longest_ = 0;
  AbsSequenceMap map_;
};
