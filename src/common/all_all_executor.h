
#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include "absl/container/flat_hash_map.h"
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
  typedef std::tuple<const Sequence*, const Sequence*, size_t> WorkItem;

  AllAllExecutor(size_t num_threads, size_t capacity,
                 AlignmentEnvironments* envs, const Parameters* params);

  void EnqueueAlignment(const WorkItem& item);

  void FinishAndOutput(const std::string& output_dir);

  void Initialize();

 private:
  std::unique_ptr<ConcurrentQueue<WorkItem>> work_queue_;

  std::vector<std::thread> threads_;
  std::thread queue_measure_thread_;

  std::atomic_uint_fast32_t num_active_threads_, id_{0};
  std::atomic<bool> run_{true};
  AlignmentEnvironments* envs_;
  const Parameters* params_;
  size_t num_threads_;

  // statistics
  std::atomic<uint64_t> num_full_alignments_{0};
  std::atomic<uint64_t> num_pass_threshold_{0};
  std::atomic<uint64_t> num_avoided_{0};
  std::vector<long int> timestamps_;
  std::vector<size_t> queue_sizes_;

  struct Match {
    int seq1_min;
    int seq1_max;
    int seq2_min;
    int seq2_max;
    double score;
    double distance;
    double variance;
    size_t cluster_size;
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

  typedef absl::flat_hash_map<GenomePair,
                              absl::flat_hash_map<SequencePair, Match>>
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
    int my_id = id_.fetch_add(1, std::memory_order_relaxed);
    auto& matches = matches_per_thread_[my_id];

    /*std::cout << string("Alignment thread spinning up with id ") +
                     std::to_string(my_id) + "\n";*/

    ProteinAligner aligner(envs_, params_);
    // std::vector<size_t> alignment_times;

    WorkItem item;
    size_t ms_wait = 0;
    int longest_wait = 0;
    bool first = true;
    while (run_.load()) {
      // read from queue, and align work item
      auto start = std::chrono::steady_clock::now();
      if (!work_queue_->pop(item)) {
        continue;
      }
      if (!first) {
        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - start)
                        .count();
        ms_wait += diff;
        if (diff > longest_wait) {
          longest_wait = diff;
        }
      } else {
        first = false;
      }

      // std::cout << "aligner thread got work\n";

      auto seq1 = std::get<0>(item);
      auto seq2 = std::get<1>(item);
      auto cluster_size = std::get<2>(item);

      ProteinAligner::Alignment alignment;
      alignment.score = 0;  // 0 score will signify not to create candidate

      auto genome_pair = std::make_pair(absl::string_view(seq1->Genome()),
                                        absl::string_view(seq2->Genome()));
      auto seq_pair = std::make_pair(seq1->GenomeIndex(), seq2->GenomeIndex());
      num_pass_threshold_++;
      if (aligner.PassesThreshold(seq1->Seq().data(), seq2->Seq().data(),
                                  seq1->Seq().size(), seq2->Seq().size())) {
        // auto t0 = std::chrono::high_resolution_clock::now();
        agd::Status s = aligner.AlignLocal(
            seq1->Seq().data(), seq2->Seq().data(), seq1->Seq().size(),
            seq2->Seq().size(), alignment);
        /*auto t1 = std::chrono::high_resolution_clock::now();
        auto duration = t1 - t0;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(duration);
        alignment_times.push_back(msec.count());*/
        num_full_alignments_++;

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
          new_match.cluster_size = cluster_size;
          matches[genome_pair][seq_pair] = new_match;
        }
      }
    }

    /*auto longest_time =
        max_element(alignment_times.begin(), alignment_times.end());
    std::cout << absl::StrCat("aligner executor thread ending, max time is ",
    *longest_time, " ms \n");*/
    /*std::cout << absl::StrCat("spent ", float(ms_wait) / 1000.0f,
        " waiting on queue\n\t and the longest wait was ", longest_wait, " ms
       \n");*/
    num_active_threads_.fetch_sub(1, std::memory_order_relaxed);
    return 0;
  }
};
