
#pragma once

#include "absl/strings/string_view.h"

class Sequence {
 public:
  Sequence(absl::string_view seq, /*std::string& coverage,*/
           const std::string& genome, uint32_t genome_size, uint32_t genome_index, uint32_t id)
      : sequence_(seq),
        //coverage_(coverage),
        genome_(genome),
        genome_size_(genome_size),
        genome_index_(genome_index),
        id_(id) {}

  Sequence(const Sequence& other) = default;
  // not sure if move is even necessary, since members
  // are basically all just POD-ish
  Sequence(Sequence&& other) = default;

  const std::string& Genome() const { return genome_; }
  uint32_t GenomeSize() const { return genome_size_; }
  uint32_t GenomeIndex() const { return genome_index_; }
  absl::string_view Seq() const { return sequence_; }
  uint32_t ID() const { return id_; }

 private:
  absl::string_view sequence_;
  //std::string& coverage_;
  const std::string& genome_;
  uint32_t genome_size_;
  uint32_t genome_index_;
  uint32_t id_;
};