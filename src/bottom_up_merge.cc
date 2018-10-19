
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
      // coverages_.push_back(string());
      // coverages_.back().resize(size);

      // cout << "Adding sequence id " << id << "\n";
      if (size > 60000) {
        cout << "over size " << size << "\n";
        exit(0);
      }
      Sequence seq(absl::string_view(data, size),  // coverages_.back(),
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

agd::Status BottomUpMerge::RunMulti(size_t num_threads,
                                    AllAllExecutor* executor,
                                    MergeExecutor* merge_executor, bool do_allall) {
  cluster_sets_left_ = sets_.size();

  // launch threads, join threads
  auto cluster_worker = [this, &merge_executor]() {
    cout << "cluster worker starting\n";
    // need own aligner per thread
    ProteinAligner aligner(aligner_->Envs(), aligner_->Params());

    // is atomic actually needed here?
    while (cluster_sets_left_.load() > 1) {
      queue_mu_.Lock();
      if (sets_.size() > 1 && cluster_sets_left_.load() > 1) {  // enough to pop
        auto s1 = std::move(sets_.front());
        sets_.pop_front();
        auto s2 = std::move(sets_.front());
        sets_.pop_front();
        cluster_sets_left_.fetch_sub(1);
        //cout << "cluster sets left: " << cluster_sets_left_.load() << "\n";
        queue_mu_.Unlock();

        // this part takes a while for larger sets
        auto t0 = std::chrono::high_resolution_clock::now();

        // eventually we may want to call single thread mergeClusters for
        // small cluster sets as it may be more efficient

        // swap so we have the larger set first, this results
        // in a larger number of smaller work items
        if (s1.Size() < s2.Size()) {
          s1.Swap(&s2);
          //cout << "swapped89\n";
        }

        auto merged_set = s1.MergeClustersParallel(s2, merge_executor);
        auto t1 = std::chrono::high_resolution_clock::now();

        auto duration = t1 - t0;
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
        //cout << "merge took: " << sec.count() << " seconds.\n";

        queue_mu_.Lock();
        sets_.push_back(std::move(merged_set));
        queue_pop_cv_.SignalAll();
        queue_mu_.Unlock();

      } else if (sets_.size() <= 1) {  // wait until enough
        while (sets_.size() <= 1 && cluster_sets_left_.load() > 1) {
          /*cout << "waiting, set size: " << sets_.size()
               << " sets left: " << cluster_sets_left_.load() << "\n";*/
          queue_pop_cv_.Wait(&queue_mu_);
        }
        if (cluster_sets_left_.load() <= 1) {
          queue_mu_.Unlock();
          break;
        }

        auto s1 = std::move(sets_.front());
        sets_.pop_front();
        auto s2 = std::move(sets_.front());
        sets_.pop_front();
        cluster_sets_left_.fetch_sub(1);
        //cout << "cluster sets left: " << cluster_sets_left_.load() << "\n";
        queue_mu_.Unlock();

        // swap so we have the larger set first, this results
        // in a larger number of smaller work items
        if (s1.Size() < s2.Size()) {
          s1.Swap(&s2);
          //cout << "swapped127\n";
        }

        // this part takes a while for larger sets
        auto t0 = std::chrono::high_resolution_clock::now();
        auto merged_set = s1.MergeClustersParallel(s2, merge_executor);
        auto t1 = std::chrono::high_resolution_clock::now();

        auto duration = t1 - t0;
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
        //cout << "merge took: " << sec.count() << " seconds.\n";

        queue_mu_.Lock();
        sets_.push_back(std::move(merged_set));
        queue_pop_cv_.SignalAll();
        queue_mu_.Unlock();

      } else {  // only one cluster set left, done.
        queue_mu_.Unlock();
        break;
      }
    }

    cout << "Cluster eval thread finishing...\n";
  };

  cout << "scheduling cluster threads, sets remaining: "
       << cluster_sets_left_.load();
  auto t0 = std::chrono::high_resolution_clock::now();
  threads_.reserve(num_threads);
  for (size_t i = 0; i < num_threads; i++) {
    threads_.push_back(std::thread(cluster_worker));
  }

  for (auto& t : threads_) {
    t.join();
  }

  auto t1 = std::chrono::high_resolution_clock::now();

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time: " << sec.count() << " seconds.\n";

  assert(sets_.size() == 1);
  // now we are all finished clustering
  auto& final_set = sets_[0];

  final_set.DumpJson();
  cout << "Total clusters: " << final_set.Size() << "\n";

  // for all clusters in final set, schedule all-all alignments with executor
  if (do_allall) {
    cout << "Scheduling all-all alignments ...\n";
    final_set.ScheduleAlignments(executor);
    cout << "Finished alignment scheduling. \n";
  } else {
    cout << "Skipping all all computation ...\n";
  }

  return agd::Status::OK();
}

agd::Status BottomUpMerge::Run(AllAllExecutor* executor, bool do_allall) {
  auto t0 = std::chrono::high_resolution_clock::now();
  while (sets_.size() > 1) {
    // dequeue 2 sets
    // merge the sets into one
    // push onto queue

    auto& s1 = sets_[0];
    auto& s2 = sets_[1];

    /*cout << "Merging cluster sets of size " << s1.Size()
      << " and " << s2.Size() << "\n";
    cout << "Cluster sets remaining: " << sets_.size() - 2 << "\n";*/
    // s1.DebugDump();
    // cout << "\nand\n";
    // s2.DebugDump();
    auto merged_set = s1.MergeClusters(s2, aligner_);

    sets_.pop_front();
    sets_.pop_front();
    sets_.push_back(std::move(merged_set));
  };

  auto t1 = std::chrono::high_resolution_clock::now();

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time: " << sec.count() << " seconds.\n";

  auto& final_set = sets_[0];

  final_set.DumpJson();

  // for all clusters in final set, schedule all-all alignments with executor
  if (do_allall) {
    final_set.ScheduleAlignments(executor);
  } else {
    cout << "Skipping all all computation ...\n";
  }

  cout << "Total clusters: " << final_set.Size() << "\n";

  return agd::Status::OK();
}
