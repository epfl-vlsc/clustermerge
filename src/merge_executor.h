
#pragma once

#include <atomic>
#include <utility>
#include "absl/synchronization/mutex.h"
#include "absl/synchronization/notification.h"
#include "alignment_environment.h"
#include "multi_notification.h"
#include "cluster_set.h"
#include "params.h"

class MergeExecutor {
 public:
  typedef std::tuple<Cluster*, ClusterSet*, MultiNotification*> WorkItem;

  MergeExecutor() = delete;
  ~MergeExecutor();

  MergeExecutor(size_t num_threads, size_t capacity,
                AlignmentEnvironments* envs, Parameters* params);

  void EnqueueMerge(const WorkItem& item);

  void Initialize();

 private:
  std::unique_ptr<ConcurrentQueue<WorkItem>> work_queue_;

  std::vector<std::thread> threads_;

  std::atomic_uint_fast32_t num_active_threads_, id_{0};
  std::atomic<bool> run_{true};
  AlignmentEnvironments* envs_;
  Parameters* params_;
  size_t num_threads_;

  // thread worker func
  void Worker();
};