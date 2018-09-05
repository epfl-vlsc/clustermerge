#pragma once

#include <list>
#include <deque>
#include "agd/agd_dataset.h"
#include "cluster_set.h"

class BottomUpMerge {
 public:
  // build one sequence, put in one cluster, put cluster in one set
  BottomUpMerge(std::vector<std::unique_ptr<agd::AGDDataset>>& datasets);

  agd::Status Run();

  void DebugDump();

 private:
  std::deque<ClusterSet> sets_;

  std::list<std::string> coverages_;
    
  // aligner object
};