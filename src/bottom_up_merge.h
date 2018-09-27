#pragma once

#include <atomic>
#include <deque>
#include <list>
#include "agd/agd_dataset.h"
#include "aligner.h"
#include "all_all_executor.h"
#include "cluster_set.h"
#include "merge_executor.h"

class BottomUpMerge {
 public:
  // build one sequence, put in one cluster, put cluster in one set
  BottomUpMerge(std::vector<std::unique_ptr<agd::AGDDataset>>& datasets,
                ProteinAligner* aligner);

  // single threaded mode
  // without mutltithread sync overhead
  agd::Status Run(AllAllExecutor* executor);

  // use multiple threads to merge clusters in parallel
  agd::Status RunMulti(size_t num_threads, AllAllExecutor* executor,
                       MergeExecutor* merge_executor);

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