#pragma once

#include <list>
#include <deque>
#include "agd/agd_dataset.h"
#include "cluster_set.h"
#include "aligner.h"
#include "all_all_executor.h"

class BottomUpMerge {
 public:
  // build one sequence, put in one cluster, put cluster in one set
  BottomUpMerge(std::vector<std::unique_ptr<agd::AGDDataset>>& datasets, ProteinAligner* aligner);

  agd::Status Run(AllAllExecutor* executor);

  void DebugDump();

 private:
  std::deque<ClusterSet> sets_;

  //std::list<std::string> coverages_;
    
  // aligner object
  ProteinAligner* aligner_;
};