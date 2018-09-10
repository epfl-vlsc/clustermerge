
#pragma once

#include <atomic>
#include <string>
#include <unordered_map>
#include "absl/strings/string_view.h"
#include "alignment_environment.h"
#include "candidate_map.h"
#include "concurrent_queue.h"
#include "params.h"
#include "sequence.h"
#include "threadpool/ThreadPool.h"

class AllAllExecutor {
 public:
  typedef std::tuple<const Sequence*, const Sequence*> WorkItem;

  AllAllExecutor(size_t num_threads, size_t capacity,
                 AlignmentEnvironments* envs, Parameters* params);

  void EnqueueAlignment(const WorkItem& item);

  void FinishAndOutput(absl::string_view output_folder);

  void Initialize();

  int Worker();
 private:
  CandidateMap candidate_map_;
  std::unique_ptr<ConcurrentQueue<WorkItem>> work_queue_;

  std::vector<std::thread> threads_;

  std::atomic_uint_fast32_t num_active_threads_, id_{0};
  std::atomic<bool> run_{true};
  AlignmentEnvironments* envs_;
  Parameters* params_;
  size_t num_threads_;

  struct Match {
    int seq1_min;
    int seq1_max;
    int seq2_min;
    int seq2_max;
    double score;
    double distance;
    double variance;
    inline bool operator==(const Match& rhs) {
      return seq1_min == rhs.seq1_min && seq1_max == rhs.seq1_max &&
             seq2_min == rhs.seq2_min && seq2_max == rhs.seq2_max &&
             score == rhs.score && distance == rhs.distance &&
             variance == rhs.variance;
    }
    inline bool operator!=(const Match& rhs) { return !(*this == rhs); }

    std::string ToString() {
      std::ostringstream s;
      s << "s1m: " << seq1_min << ", s1M: " << seq1_max << ", s2m: " << seq2_min
        << ", s2M: " << seq2_max << ", score: " << score
        << ", dist: " << distance << ", var: " << variance;
      return s.str();
    }
  };

  typedef std::unordered_map<
      GenomePair, std::unordered_map<SequencePair, Match, PairHash>, PairHash>
      ResultMap;
  // candidate map prevents dups, we store the actual matches here
  // each thread gets its own map, to avoid any sync here
  std::vector<ResultMap> matches_per_thread_;

};