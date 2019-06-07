
#include "merge_batch.h"

void MergeBatch(std::deque<ClusterSet> &sets_to_merge,
                ProteinAligner *aligner) {
  while (sets_to_merge.size() > 1) {
    // dequeue 2 sets
    // merge the sets into one
    // push onto queue

    auto &s1 = sets_to_merge[0];
    auto &s2 = sets_to_merge[1];

    /*cout << "Merging cluster sets of size " << s1.Size()
      << " and " << s2.Size() << "\n";
    cout << "Cluster sets remaining: " << sets_.size() - 2 << "\n";*/
    // s1.DebugDump();
    // cout << "\nand\n";
    // s2.DebugDump();
    auto merged_set = s1.MergeClusters(s2, aligner);

    sets_to_merge.pop_front();
    sets_to_merge.pop_front();
    sets_to_merge.push_back(std::move(merged_set));
  };
}