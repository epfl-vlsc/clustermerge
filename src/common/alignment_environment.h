
#pragma once
#include "src/agd/errors.h"
extern "C" {
#include "swps3/EstimatePam.h"
#include "swps3/swps3.h"
}
#include <vector>

#include "json.hpp"

using json = nlohmann::json;

#define DIMSIZE MATRIX_DIM
#define MDIM DIMSIZE* DIMSIZE

struct AlignmentEnvironment {
  double gap_open;
  double gap_extend;
  double pam_distance;
  double threshold;
  double* matrix = nullptr;  // 26 by 26 matrix of scores
  int8_t gap_open_int8;
  int8_t gap_ext_int8;
  int8_t* matrix_int8 = nullptr;
  int16_t gap_open_int16;
  int16_t gap_ext_int16;
  int16_t* matrix_int16 = nullptr;
};

class AlignmentEnvironments {
  // pointers here own no data
 public:
  AlignmentEnvironments() {}
  // no allow copy, iz bad
  AlignmentEnvironments(const AlignmentEnvironments& envs) = delete;

  void EstimPam(char* seq1, char* seq2, int len, double result[3]) const;
  const AlignmentEnvironment& FindNearest(double pam) const;
  const AlignmentEnvironment& LogPamEnv() const;
  const AlignmentEnvironment& JustScoreEnv() const;
  const AlignmentEnvironment& LogPamJustScoreEnv() const;

  void InitFromJSON(const json& logpam_json, const json& matrices_json,
                    double threshold);

  // override the just score env to use blosum62
  void UseBlosum(const json& blosum_json, double threshold);

 private:
  void CreateDayMatrices(std::vector<double>& gap_open,
                         std::vector<double>& gap_ext,
                         std::vector<double>& pam_dist,
                         std::vector<double*>& matrices);
  std::vector<AlignmentEnvironment> envs_;
  DayMatrix* day_matrices_;
  // double* logpam1_matrix_;
  AlignmentEnvironment logpam_env_;
  AlignmentEnvironment just_score_env_;
  AlignmentEnvironment logpam_just_score_;
};
