#pragma once
#include <cstddef>
#include <memory>
#include "data.h"
#include "status.h"

namespace agd {

class Buffer : public Data {
 private:
  std::unique_ptr<char[]> buf_{nullptr};
  std::size_t size_ = 0, allocation_ = 0, extend_extra_ = 0;

 public:
  ~Buffer() override = default;
  Buffer(const Buffer& other) = delete;
  void operator=(const Buffer& other) = delete;

  Buffer(Buffer&&) = default;
  Buffer& operator=(Buffer&&) = default;
  Buffer(decltype(size_) initial_size = 2 * 1024 * 1024,
         decltype(size_) extend_extra = 8 * 1024 * 1024);

  Status WriteBuffer(const char* content, std::size_t content_size);
  Status AppendBuffer(const char* content, std::size_t content_size);

  // ensures that the total capacity is at least `capacity`
  void reserve(decltype(allocation_) capacity);

  // resizes the actual size to to total_size
  // returns an error if total_size > allocation_
  // should call `reserve` first to be safe
  void resize(decltype(size_) total_size);

  // extends the current allocation by `extend_size`
  // does not affect size_
  void extend_allocation(decltype(size_) extend_size);

  // extends the current size by extend_size
  // returns an error if size_+extend_size > allocation
  // should call `extend_allocation` first to be safe
  void extend_size(decltype(size_) extend_size);

  char& operator[](std::size_t idx) const;

  char* release_raw() { return buf_.release(); }

  void reset();
  virtual const char* data() const override;
  virtual std::size_t size() const override;
  virtual char* mutable_data() override;
  std::size_t capacity() const;
};
}  // namespace agd
