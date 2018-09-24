
#pragma once

#include <algorithm>
#include <atomic>
#include <iostream>
#include <string>
#include <unordered_map>
#include <thread>
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "aligner.h"
#include "alignment_environment.h"
#include "candidate_map.h"
#include "concurrent_queue.h"
#include "params.h"
#include "sequence.h"

class AllAllExecutor {
 public:
  AllAllExecutor() = delete;
  typedef std::tuple<const Sequence*, const Sequence*> WorkItem;

  AllAllExecutor(size_t num_threads, size_t capacity,
                 AlignmentEnvironments* envs, Parameters* params);

  void EnqueueAlignment(const WorkItem& item);

  void FinishAndOutput(absl::string_view output_folder);

  void Initialize();

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

  bool PassesLengthConstraint(const ProteinAligner::Alignment& alignment,
                              int seq1_len, int seq2_len) {
    float min_alignment_len =
        std::min(float(alignment.seq1_length), float(alignment.seq2_length));
    float max_min_seq_len =
        std::max(30.0f, 0.3f * float(std::min(seq1_len, seq2_len)));
    return min_alignment_len >= max_min_seq_len;
  }

  bool PassesScoreConstraint(const Parameters* params, int score) {
    return score >= params->min_score;
  }

  int Worker() {
    // int test = id_.load();
    // cout << string("cur val of test is ") + std::to_string(test);

    int my_id = id_.fetch_add(1, std::memory_order_relaxed);
    auto& matches = matches_per_thread_[my_id];

    /*std::cout << string("Alignment thread spinning up with id ") +
                     std::to_string(my_id) + "\n";*/

    // int capacity = work_queue_->capacity();

    ProteinAligner aligner(envs_, params_);

    while (run_.load()) {
      // read from queue, and align work item
      WorkItem item;
      if (!work_queue_->pop(item)) {
        continue;
      }

      // std::cout << "aligner thread got work\n";

      auto seq1 = std::get<0>(item);
      auto seq2 = std::get<1>(item);

      if (seq1->GenomeSize() > seq2->GenomeSize() ||
          ((seq1->GenomeSize() == seq2->GenomeSize()) &&
           seq1->Genome() > seq2->Genome())) {
        // seq1 = sequence;
        std::swap(seq1, seq2);
      } else {
        // seq2 = sequence;
      }

      ProteinAligner::Alignment alignment;
      alignment.score = 0;  // 0 score will signify not to create candidate

      if (seq1->Genome() == seq2->Genome() &&
          seq1->GenomeIndex() == seq2->GenomeIndex()) {
        // not sure if this can actually happen yet, but no need to align
        // against self
        continue;
      }

      if (seq1->Genome() == seq2->Genome() &&
          seq1->GenomeIndex() > seq2->GenomeIndex()) {
        std::swap(seq1, seq2);
      }

      auto genome_pair = std::make_pair(seq1->Genome(), seq2->Genome());
      auto seq_pair = std::make_pair(seq1->GenomeIndex(), seq2->GenomeIndex());
      if (!candidate_map_.ExistsOrInsert(genome_pair, seq_pair)) {
        if (aligner.PassesThreshold(seq1->Seq().data(), seq2->Seq().data(),
                                    seq1->Seq().size(), seq2->Seq().size())) {
          agd::Status s = aligner.AlignLocal(
              seq1->Seq().data(), seq2->Seq().data(), seq1->Seq().size(),
              seq2->Seq().size(), alignment);

          if (PassesLengthConstraint(alignment, seq1->Seq().size(),
                                     seq2->Seq().size()) &&
              PassesScoreConstraint(params_, alignment.score)) {
            Match new_match;
            new_match.seq1_min = alignment.seq1_min;
            new_match.seq1_max = alignment.seq1_max;
            new_match.seq2_min = alignment.seq2_min;
            new_match.seq2_max = alignment.seq2_max;
            new_match.score = alignment.score;
            new_match.variance = alignment.pam_variance;
            new_match.distance = alignment.pam_distance;
            matches[genome_pair][seq_pair] = new_match;
          }
        }
      }
      // else, we already aligned these two seqs, done
    }

    std::cout << "aligner executor thread ending.\n";
    num_active_threads_.fetch_sub(1, std::memory_order_relaxed);
    return 0;
  }
};
