
#pragma once

#include "absl/strings/string_view.h"

class Sequence {
 public:
  Sequence(absl::string_view seq, std::string& coverage,
           const std::string& genome, uint32_t genome_size)
      : sequence_(seq),
        coverage_(coverage),
        genome_(genome),
        genome_size_(genome_size) {}

  Sequence(const Sequence& other) = default;
  // not sure if move is even necessary, since members
  // are basically all just POD-ish
  Sequence(Sequence&& other) = default;

  const std::string& Genome() const { return genome_; }
  absl::string_view Seq() const { return sequence_; }

 private:
  absl::string_view sequence_;
  std::string& coverage_;
  const std::string& genome_;
  uint32_t genome_size_;
};