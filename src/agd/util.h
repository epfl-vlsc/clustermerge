#pragma once

#include <vector>
#include "errors.h"
#include "data.h"

namespace agd {

  template <typename T>
    inline void safe_reserve(std::vector<T> &v, const std::size_t ideal_length, const size_t additive_increase_bonus_size = 2 * 1024 * 1024) {
    if (v.capacity() < ideal_length) {
      v.reserve(ideal_length + additive_increase_bonus_size);
    }
  }
} // namespace agd {
