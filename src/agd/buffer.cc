#include "buffer.h"
#include <cstddef>
#include <cstring>
#include <memory>
#include "util.h"

namespace agd {

using namespace std;
using namespace errors;

Buffer::Buffer(size_t initial_size, size_t extend_extra)
    : size_(0), allocation_(initial_size), extend_extra_(extend_extra) {
  // FIXME in an ideal world, allocation_ should be checked to be positive
  buf_ = std::make_unique<char[]>(allocation_);
  num_allocs_++;
}

Status Buffer::WriteBuffer(const char* content, size_t content_size) {
  if (allocation_ < content_size) {
    allocation_ = content_size + extend_extra_;
    buf_.reset(new char[allocation_]());  // reset() -> old buf will be deleted
    num_allocs_++;
  }
  memcpy(buf_.get(), content, content_size);
  size_ = content_size;
  return Status::OK();
}

Status Buffer::AppendBuffer(const char* content, size_t content_size) {
  auto old_size = size_;
  extend_size(content_size);
  memcpy(&buf_.get()[old_size], content, content_size);
  return Status::OK();
}

const char* Buffer::data() const { return buf_.get(); }

char* Buffer::mutable_data() { return buf_.get(); }

size_t Buffer::size() const { return size_; }

void Buffer::reset() { size_ = 0; }

char& Buffer::operator[](size_t idx) const { return buf_[idx]; }

void Buffer::reserve(size_t capacity) {
  if (capacity > allocation_) {
    allocation_ = capacity + extend_extra_;
    decltype(buf_) a = std::make_unique<char[]>(allocation_);
    memcpy(a.get(), buf_.get(), size_);
    buf_.swap(a);
    num_allocs_++;
  }
}

void Buffer::resize(size_t total_size) {
  reserve(total_size);
  size_ = total_size;
}

void Buffer::extend_allocation(size_t extend_size) {
  reserve(allocation_ + extend_size);
}

void Buffer::extend_size(size_t extend_size) {
  return resize(size_ + extend_size);
}

size_t Buffer::capacity() const { return allocation_; }

std::size_t Buffer::num_allocs() const { return num_allocs_; }

}  // namespace agd