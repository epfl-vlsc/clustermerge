#pragma once

#include <tuple>
#include <sstream>
#include "src/common/sequence.h"

struct __attribute__((__packed__)) Match {
  int seq1_min;
  int seq1_max;
  int seq2_min;
  int seq2_max;
  double score;
  double distance;
  double variance;
  size_t cluster_size;
  inline bool operator==(const Match& rhs) {
    return seq1_min == rhs.seq1_min && seq1_max == rhs.seq1_max &&
           seq2_min == rhs.seq2_min && seq2_max == rhs.seq2_max &&
           score == rhs.score && distance == rhs.distance &&
           variance == rhs.variance;
  }
  inline bool operator!=(const Match& rhs) { return !(*this == rhs); }

  std::string ToString() {
    std::ostringstream s;
    s << "s1m: " << seq1_min << ", s1M: " << seq1_max << ", s2m: " << seq2_min
      << ", s2M: " << seq2_max << ", score: " << score << ", dist: " << distance
      << ", var: " << variance;
    return s.str();
  }
};

// dist needs to include which seqs were aligned
struct __attribute__((__packed__)) DistMatchResult {
  Match m;
  int abs_seq_1;
  int abs_seq_2;
};

class AllAllBase {

public:

  typedef std::tuple<const Sequence*, const Sequence*, size_t> WorkItem;
  
  virtual void EnqueueAlignment(const WorkItem& item) = 0;
};