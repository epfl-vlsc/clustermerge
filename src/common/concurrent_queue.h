#pragma once

#include <queue>
#include <utility>
#include "absl/synchronization/mutex.h"

// a class wrapping STL queue
// thread-safe, limited buffer capacity, blocks on push()
// to a full queue.

template <typename T>
class ConcurrentQueue {
 public:
  ConcurrentQueue(size_t capacity);
  ~ConcurrentQueue() = default;

  // return true if pushed, false otherwise
  // will block until pushed if block_ is true
  bool push(const T& item);
  bool push(T&& item);
  // return true if success and item is valid, false otherwise
  bool pop(T& item);

  bool drop_if_equal(T& item);

  bool peek(T& item);

  // unblock the queue, notify all threads
  void unblock();
  // set blocking behavior
  void set_block();

  bool empty() const;
  size_t capacity() const;
  size_t size() const;

  uint64_t num_pop_waits();
  uint64_t num_push_waits();
  uint64_t num_peek_waits();

  // these are for iterating the queue to do the checkpointing
  // these are not threadsafe
  typename std::deque<T>::const_iterator begin() const { return queue_.begin(); }
  typename std::deque<T>::const_iterator end() const { return queue_.end(); }

 private:
  // mutex to protect the queue
  mutable absl::Mutex mu_;
  // cond vars for block/wait/notify on queue push/pop
  mutable absl::CondVar queue_pop_cv_;
  mutable absl::CondVar queue_push_cv_;
  std::deque<T> queue_;
  size_t capacity_;
  // block on calls to push, pop
  bool block_ = true;
  uint64_t num_pop_waits_ = 0;
  uint64_t num_push_waits_ = 0;
  uint64_t num_peek_waits_ = 0;
  uint64_t num_push_ = 0;

};

template <typename T>
class ScopeDropIfEqual {
 public:
  ScopeDropIfEqual(ConcurrentQueue<T>& queue, T& item)
      : queue_(queue), item_(item) {}
  ~ScopeDropIfEqual() { queue_.drop_if_equal(item_); }

 private:
  ConcurrentQueue<T>& queue_;
  T& item_;
};

template <typename T>
bool ConcurrentQueue<T>::peek(T& item) {
  bool popped = false;
  {
    absl::MutexLock l(&mu_);

    if (queue_.empty() && block_) {
      num_peek_waits_++;
      while (queue_.empty() && block_) {
        // queue_pop_cv_.wait(l);
        queue_pop_cv_.Wait(&mu_);
      }
    }

    if (!queue_.empty()) {
      item = std::move(queue_.front());
      popped = true;
    }
  }
  if (popped) queue_pop_cv_.Signal();
  return popped;
}

template <typename T>
bool ConcurrentQueue<T>::drop_if_equal(T& item) {
  bool ret = false;
  {
    absl::MutexLock l(&mu_);
    if (!queue_.empty() && queue_.front() == item) {
      queue_.pop_front();
      ret = true;
    }
  }
  queue_push_cv_.Signal();
  return ret;
}

template <typename T>
bool ConcurrentQueue<T>::pop(T& item) {
  bool popped = false;
  {
    absl::MutexLock l(&mu_);

    if (queue_.empty() && block_) {
      num_pop_waits_++;
      while (queue_.empty() && block_) {
        queue_pop_cv_.Wait(&mu_);
      }
    }

    if (!queue_.empty()) {
      item = std::move(queue_.front());
      queue_.pop_front();
      popped = true;
    }
  }
  if (popped) {
    // tell someone blocking on write they can now write to the queue
    queue_push_cv_.Signal();
    return true;
  } else
    return false;
}

template <typename T>
bool ConcurrentQueue<T>::push(const T& item) {
  bool pushed = false;
  {
    absl::MutexLock l(&mu_);
    // we block until something pops and makes room for us
    // unless blocking is set to false
    if (queue_.size() == capacity_ && block_) {
      num_push_waits_++;
      while (queue_.size() == capacity_ && block_) {
        queue_push_cv_.Wait(&mu_);
      }
    }

    if (queue_.size() < capacity_) {
      queue_.push_back(item);
      pushed = true;
    }
  }

  if (pushed) {
    // tell someone blocking on read they can now read from the queue
    // TODO maybe notify_all is better here? If so, good to drop the notify_one
    // in peek
    queue_pop_cv_.Signal();
    num_push_++;
    return true;
  } else
    return false;
}

template <typename T>
bool ConcurrentQueue<T>::push(T&& item) {
  bool pushed = false;
  {
    absl::MutexLock l(&mu_);
    // we block until something pops and makes room for us
    // unless blocking is set to false
    if (queue_.size() == capacity_ && block_) {
      num_push_waits_++;
      while (queue_.size() == capacity_ && block_) {
        queue_push_cv_.Wait(&mu_);
      }
    }

    if (queue_.size() < capacity_) {
      queue_.push_back(std::move(item));
      pushed = true;
    }
  }

  if (pushed) {
    // tell someone blocking on read they can now read from the queue
    // TODO maybe notify_all is better here? If so, good to drop the notify_one
    // in peek
    queue_pop_cv_.Signal();
    num_push_++;
    return true;
  } else
    return false;
}

template <typename T>
void ConcurrentQueue<T>::unblock() {
  {
    absl::MutexLock l(&mu_);
    block_ = false;
  }

  queue_push_cv_.SignalAll();
  queue_pop_cv_.SignalAll();
}

template <typename T>
void ConcurrentQueue<T>::set_block() {
  absl::MutexLock l(&mu_);
  block_ = true;
}

template <typename T>
ConcurrentQueue<T>::ConcurrentQueue(size_t capacity) : capacity_(capacity) {}

template <typename T>
bool ConcurrentQueue<T>::empty() const {
  return queue_.empty();
}

template <typename T>
size_t ConcurrentQueue<T>::capacity() const {
  return capacity_;
}

template <typename T>
size_t ConcurrentQueue<T>::size() const {
  return queue_.size();
}

template <typename T>
uint64_t ConcurrentQueue<T>::num_pop_waits() {
  return num_pop_waits_;
}

template <typename T>
uint64_t ConcurrentQueue<T>::num_push_waits() {
  return num_push_waits_;
}

template <typename T>
uint64_t ConcurrentQueue<T>::num_peek_waits() {
  return num_peek_waits_;
}
