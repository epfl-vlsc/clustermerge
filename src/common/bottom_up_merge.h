#pragma once

#include <atomic>
#include <deque>
#include <list>
#include "aligner.h"
#include "all_all_executor.h"
#include "cluster_set.h"
#include "merge_executor.h"
#include "src/dataset/dataset.h"

class BottomUpMerge {
 public:
  // build one sequence, put in one cluster, put cluster in one set
  BottomUpMerge(std::vector<std::unique_ptr<Dataset>>& datasets,
                ProteinAligner* aligner);

  // single threaded mode
  // without mutltithread sync overhead
  agd::Status Run(AllAllExecutor* executor, size_t dup_removal_threshold,
                  bool do_allall);

  // use multiple threads to merge clusters in parallel
  agd::Status RunMulti(size_t num_threads, size_t dup_removal_threshold,
                       AllAllExecutor* executor, MergeExecutor* merge_executor,
                       bool do_allall);

  void DebugDump();

 private:
  std::deque<ClusterSet> sets_;

  // threads to run cluster mergers in parallel
  std::vector<std::thread> threads_;

  // aligner object
  ProteinAligner* aligner_;

  // mutex and sync vars
  absl::Mutex queue_mu_;
  std::atomic<uint32_t> cluster_sets_left_;
  mutable absl::CondVar queue_pop_cv_;
};