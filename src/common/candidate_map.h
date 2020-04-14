
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <chrono>
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "bloom_filter/bloom_filter.hpp"


typedef std::pair<absl::string_view, absl::string_view>
    GenomePair;
typedef std::pair<int, int> SequencePair;
typedef std::pair<uint32_t, uint32_t> AbsSequencePair;

class CandidateMap {
  typedef absl::flat_hash_set<AbsSequencePair> AbsSequenceMap;

 public:
  CandidateMap() = default;
  CandidateMap(size_t size) {
    //map_.reserve(size);
    // BLOOM FILTER PARAMS
    bl_params_.projected_element_count = 10000000;
    bl_params_.false_positive_probability = 0.0001f;
    if (!bl_params_) {
      std::cout << "invalid BLFilter params!\n";
      exit(0);
    }
    bl_params_.compute_optimal_parameters();
    
    bl_filter_.reset(new bloom_filter(bl_params_));
    std::cout << "[CandMap] the table size is " << bl_filter_->size() << "\n";
  }

  ~CandidateMap() {
    //std::cout << "Candidate map was rehashed " << times_rehashed << " times.\n";
    //std::cout << "The longest insert time was: " << longest_ << " ms.\n";
    std::cout << "[CandMap] There were " << lookups_ << " lookups.\n";
    std::cout << "[CandMap] There were " << disagreements_ << " disagreements between Map and Filter.\n";
  }
  // do not copy or move
  CandidateMap(CandidateMap&& other) = delete;
  CandidateMap(const CandidateMap& other) = delete;
  CandidateMap& operator=(const CandidateMap& other) = delete;

  // returns true, or false and inserts key
  /*bool ExistsOrInsertMap(const AbsSequencePair& s) {
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
  }*/
  
  bool ExistsOrInsert(const AbsSequencePair& s) {
    lookups_++;
    bool ret_filter;
    if (bl_filter_->contains(s)) {
      ret_filter = true;
    } else {
      bl_filter_->insert(s);
      ret_filter = false;
    }

    // instrumentation to guage how well the filter works
    /*bool ret_map = ExistsOrInsertMap(s);
    if (ret_filter != ret_map) disagreements_++;*/
    
    return ret_filter;
  }

 private:
  uint64_t disagreements_ = 0;
  uint64_t lookups_ = 0;
  int times_rehashed = 0;
  int longest_ = 0;
  //AbsSequenceMap map_;
  bloom_parameters bl_params_;
  std::unique_ptr<bloom_filter> bl_filter_;
};
