
#include "all_all_executor.h"
#include <iostream>
#include "aligner.h"

using std::cout;
using std::get;
using std::make_pair;

AllAllExecutor::AllAllExecutor(size_t num_threads, size_t capacity,
                               AlignmentEnvironments* envs, Parameters* params)
    : envs_(envs), params_(params) {
  thread_pool_.reset(new ThreadPool(num_threads));
  work_queue_.reset(new ConcurrentQueue<WorkItem>(capacity));
  matches_per_thread_.resize(num_threads);

  auto aligner_func = [this]() {
    // std::chrono::high_resolution_clock::time_point start_time =
    // std::chrono::high_resolution_clock::now();
    int my_id = id_.fetch_add(1, std::memory_order_relaxed);
    auto& matches = matches_per_thread_[my_id];

    std::cout << "Alignment thread spinning up with id " << my_id;

    // int capacity = work_queue_->capacity();

    ProteinAligner aligner(envs_, params_);

    while (run_.load()) {
      // read from queue, and align work item
      WorkItem item;
      if (!work_queue_->pop(item)) {
        continue;
      }

      auto seq1 = get<0>(item);
      auto seq2 = get<1>(item);

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

      auto genome_pair = make_pair(seq1->Genome(), seq2->Genome());
      auto seq_pair = make_pair(seq1->GenomeIndex(), seq2->GenomeIndex());
      if (!candidate_map_.ExistsOrInsert(genome_pair, seq_pair)) {
        if (aligner.PassesThreshold(seq1->Seq().data(), seq2->Seq().data(),
                                    seq1->Seq().size(), seq2->Seq().size())) {
          agd::Status s = aligner.AlignLocal(
              seq1->Seq().data(), seq2->Seq().data(), seq1->Seq().size(),
              seq2->Seq().size(), alignment);
        }

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
      // else, we already aligned these two seqs, done
    }

    cout << "aligner executor thread ending.\n";
    num_active_threads_.fetch_sub(1, std::memory_order_relaxed);
  };

  num_active_threads_ = num_threads;
  for (int i = 0; i < num_threads; i++) thread_pool_->enqueue(aligner_func);
}