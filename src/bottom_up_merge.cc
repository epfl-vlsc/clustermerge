
#include "bottom_up_merge.h"
#include <iostream>

using std::cout;
using std::string;

BottomUpMerge::BottomUpMerge(
    std::vector<std::unique_ptr<agd::AGDDataset>>& datasets,
    ProteinAligner* aligner) {
  aligner_ = aligner;

  const char* data;
  size_t size;
  uint32_t id = 0;
  for (auto& dataset : datasets) {
    cout << "Parsing dataset " << dataset->Name() << " ...\n";
    agd::AGDDataset::ColumnIterator iter;
    auto s = dataset->Column("prot", &iter);

    if (!s.ok()) {
      cout << "dataset " << dataset->Name()
           << " had no prot column, skipping ...";
      continue;
    }

    s = iter.GetNextRecord(&data, &size);
    uint32_t genome_index = 0;
    while (s.ok()) {
      coverages_.push_back(string());
      coverages_.back().resize(size);

      cout << "Adding sequence id " << id << "\n";
      Sequence seq(absl::string_view(data, size), coverages_.back(),
                   dataset->Name(), dataset->Size(), genome_index++, id++);

      ClusterSet cs(seq);

      sets_.push_back(std::move(cs));
      s = iter.GetNextRecord(&data, &size);
    }
  }
}

void BottomUpMerge::DebugDump() {
  cout << "Dumping merger ... \n";
  for (const auto& cs : sets_) {
    cs.DebugDump();
  }
}

agd::Status BottomUpMerge::Run(AllAllExecutor* executor) {
  while (sets_.size() > 1) {
    // dequeue 2 sets
    // merge the sets into one
    // push onto queue

    auto& s1 = sets_[0];
    auto& s2 = sets_[1];

    auto merged_set = s1.MergeClusters(s2, aligner_);

    sets_.pop_front();
    sets_.pop_front();
    sets_.push_back(std::move(merged_set));
  };

  auto& final_set = sets_[0];

  // for all clusters in final set, schedule all-all alignments with executor

  return agd::Status::OK();
}