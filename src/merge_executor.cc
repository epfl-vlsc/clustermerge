#include "merge_executor.h"

using std::cout;
using namespace std::literals::chrono_literals;

MergeExecutor::MergeExecutor(size_t num_threads, size_t capacity,
                             AlignmentEnvironments* envs, Parameters* params)
    : envs_(envs), params_(params), num_threads_(num_threads) {
  work_queue_.reset(new ConcurrentQueue<WorkItem>(capacity));

  num_active_threads_ = num_threads;
  threads_.reserve(num_threads_);
  for (size_t i = 0; i < num_threads_; i++) {
    threads_.push_back(std::thread(&MergeExecutor::Worker, this));
  }
  cout << "merge executor started threads\n";
}

MergeExecutor::~MergeExecutor() {
  cout << "MergeExecutor waiting for work queue to empty\n";
  while (!work_queue_->empty()) {
    std::this_thread::sleep_for(1s);
    cout << "MergeExecutor work queue has " << work_queue_->size() << " entries left\n";
  }
  cout << "MergeExecutor Queue emptied, unblocking...\n";
  run_ = false;
  work_queue_->unblock();
  for (auto& f : threads_) {
    f.join();
  }
  cout << "All threads finished.\n";
  cout << "Num pass threshold alignments: " << num_alignments_ << "\n";
}

void MergeExecutor::EnqueueMerge(const WorkItem& item) {
  if (!work_queue_->push(item)) {
    cout << "failed to push work item to queue.\n";
  }
}

void MergeExecutor::Worker() {
  int my_id = id_.fetch_add(1, std::memory_order_relaxed);

  std::cout << absl::StrCat("merger thread spinning up with id ",
                   my_id, "\n");

  ProteinAligner aligner(envs_, params_);

  while (run_.load()) {
    // read from queue, and align work item
    WorkItem item;
    if (!work_queue_->pop(item)) {
      continue;
    }

    auto* cluster = std::get<0>(item);
    auto* cluster_set = std::get<1>(item);
    auto* notification = std::get<2>(item);

    //cout << "merge thread merging ...\n";

    cluster_set->MergeClusterLocked(cluster, &aligner);

    notification->Notify();

  }

  num_alignments_ += aligner.NumAlignments();

}
