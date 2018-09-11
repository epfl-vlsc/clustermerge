
#include "all_all_executor.h"
#include <iostream>
#include "aligner.h"

using std::cout;
using std::get;
using std::make_pair;
using std::string;

void AllAllExecutor::EnqueueAlignment(const WorkItem& item) {
  if (!work_queue_->push(item)) {
    cout << "failed to push work item to queue.\n";
  }
}

void AllAllExecutor::FinishAndOutput(absl::string_view output_folder) {
  while (!work_queue_->empty()) {
    ;
    ;
  }
  cout << "Queue emptied, unblocking...\n";
  run_ = false;
  work_queue_->unblock();
  for (auto& f : threads_) {
    f.join();
  }
  cout << "All threads finished.\n";

  cout << "Output to dir " << output_folder << "\n";
}

AllAllExecutor::AllAllExecutor(size_t num_threads, size_t capacity,
                               AlignmentEnvironments* envs, Parameters* params)
    : envs_(envs), params_(params), num_threads_(num_threads) {
  work_queue_.reset(new ConcurrentQueue<WorkItem>(capacity));
  matches_per_thread_.resize(num_threads);

  cout << "Start executor, id is " << id_.load() << "\n";

  num_active_threads_ = num_threads;
}

void AllAllExecutor::Initialize() {
  threads_.reserve(num_threads_);
  for (size_t i = 0; i < num_threads_; i++) {
    threads_.push_back(std::thread(&AllAllExecutor::Worker, this));
  }
}