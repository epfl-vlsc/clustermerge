
#include "alignment_environment.h"
#include <iostream>
#include <string>

extern "C" {
#include "swps3/DynProgr_scalar.h"
#include "swps3/DynProgr_sse_double.h"
#include "swps3/swps3.h"
}

// len is seq1 len
void AlignmentEnvironments::EstimPam(char* seq1, char* seq2, int len,
                                     double result[3]) const {
  EstimatePam(seq1, seq2, len, day_matrices_, (int)envs_.size(),
              logpam_env_.matrix, result);
}

using namespace std;

double ByteFactor(double* matrix, double threshold) {
  auto min = fabs(*std::min_element(matrix, matrix + MDIM));
  return 255.0f / (threshold + min);
}

double ShortFactor(double threshold) { return 65535.0 / threshold; }

int8_t ScaleByte(double value, double factor) {
  auto ret = ceil(value * factor);
  if (ret < -128) return -128;
  return int8_t(ret);
}

int16_t ScaleShort(double value, double factor) {
  auto ret = ceil(value * factor);
  if (ret < -32767) return -32767;
  return int16_t(ret);
}

int8_t* CreateScaled(double* matrix, double threshold, double gap,
                     double gap_ext, int8_t& gapi8, int8_t& gap_exti8) {
  auto factor = ByteFactor(matrix, threshold);
  int8_t* ret = new int8_t[MDIM];
  gapi8 = ScaleByte(gap, factor);
  gap_exti8 = ScaleByte(gap_ext, factor);
  for (size_t i = 0; i < MDIM; i++) {
    ret[i] = ScaleByte(matrix[i], factor);
  }
  return ret;
}

int16_t* CreateScaled(double* matrix, double threshold, double gap,
                      double gap_ext, int16_t& gapi16, int16_t& gap_exti16) {
  auto factor = ShortFactor(threshold);
  int16_t* ret = new int16_t[MDIM];
  gapi16 = ScaleShort(gap, factor);
  gap_exti16 = ScaleShort(gap_ext, factor);
  for (size_t i = 0; i < MDIM; i++) {
    ret[i] = ScaleShort(matrix[i], factor);
  }
  return ret;
}

const AlignmentEnvironment& AlignmentEnvironments::FindNearest(
    double pam) const {
  size_t i = 0;
  while (pam - envs_[i].pam_distance > 0.0f && i < envs_.size()) i++;
  if (i == envs_.size())
    return envs_[i - 1];
  else {
    if (fabs(envs_[i].pam_distance - pam) <
        fabs(envs_[i - 1].pam_distance - pam))
      return envs_[i];
    else
      return envs_[i - 1];
  }
}

const AlignmentEnvironment& AlignmentEnvironments::LogPamEnv() const {
  return logpam_env_;
}

const AlignmentEnvironment& AlignmentEnvironments::JustScoreEnv() const {
  return just_score_env_;
}

void AlignmentEnvironments::CreateDayMatrices(vector<double>& gap_open,
                                              vector<double>& gap_ext,
                                              vector<double>& pam_dist,
                                              vector<double*>& matrices) {
  day_matrices_ =
      createDayMatrices(&gap_open[0], &gap_ext[0], &pam_dist[0],
                        (long long*)&matrices[0], (int)gap_open.size() - 1);
}

void ReadJsonEnv(const json& json_env, AlignmentEnvironment* env) {
  env->gap_open = json_env["gap_open"];
  env->gap_extend = json_env["gap_ext"];
  env->pam_distance = json_env["pam_distance"];
  env->threshold = 85.0f;

  env->matrix = new double[MDIM];
  std::fill(env->matrix, env->matrix + MDIM, 0.0f);
  const auto& column_order = json_env["column_order"];
  const auto& compact_matrix = json_env["scores"];

  for (size_t i = 0; i < column_order.size(); i++) {
    if (compact_matrix[i].size() != column_order.size()) {
      cout << "column order doesnt match!! " << compact_matrix[i].dump()
           << " \n"
           << column_order << "\n";
      exit(0);
    }
    for (size_t j = 0; j < column_order.size(); j++) {
      // cout << "column order first char is " <<
      // column_order[i].get<string>().front() << "\n";
      env->matrix[(column_order[i].get<string>().front() - 'A') * DIMSIZE +
                  column_order[j].get<string>().front() - 'A'] =
          compact_matrix[i][j];
    }
  }

  env->matrix_int16 =
      CreateScaled(env->matrix, env->threshold, env->gap_open, env->gap_extend,
                   env->gap_open_int16, env->gap_ext_int16);
  env->matrix_int8 =
      CreateScaled(env->matrix, env->threshold, env->gap_open, env->gap_extend,
                   env->gap_open_int8, env->gap_ext_int8);
}

void AlignmentEnvironments::InitFromJSON(const json& logpam_json,
                                         const json& matrices_json,
                                         double threshold) {
  auto total_matrices = matrices_json["matrices"].size();
  envs_.resize(total_matrices);
  size_t i = 0;
  for (const auto& env : matrices_json["matrices"]) {
    ReadJsonEnv(env, &envs_[i]);
    i++;
  }

  ReadJsonEnv(logpam_json, &logpam_env_);

  just_score_env_.gap_open = -37.64 + 7.434 * log10(224);
  just_score_env_.gap_extend = -1.3961;
  just_score_env_.pam_distance = 224;
  just_score_env_.threshold = threshold * 0.75f;
  // cout << "just score threshold is " << just_score_env_.threshold << "\n";

  just_score_env_.matrix = new double[MDIM];
  CreateOrigDayMatrix(logpam_env_.matrix, 224, just_score_env_.matrix);
  just_score_env_.matrix_int8 =
      CreateScaled(just_score_env_.matrix, just_score_env_.threshold,
                   just_score_env_.gap_open, just_score_env_.gap_extend,
                   just_score_env_.gap_open_int8, just_score_env_.gap_ext_int8);
  just_score_env_.matrix_int16 = CreateScaled(
      just_score_env_.matrix, just_score_env_.threshold,
      just_score_env_.gap_open, just_score_env_.gap_extend,
      just_score_env_.gap_open_int16, just_score_env_.gap_ext_int16);

  /*cout << "just score int 16\n";
  for (size_t i = 0; i < DIMSIZE; i++) {
    for (size_t j = 0; j < DIMSIZE; j++) {
      printf("%d ", just_score_env_.matrix_int16[i*DIMSIZE + j]);
    }
    printf("\n");
  }*/

  for (size_t i = 0; i < envs_.size(); i++) {
    if (envs_[i].matrix == nullptr) {
      cout << "NULL MAT POINTER!!";
    }
  }
  vector<double> gaps, gap_extends, pam_dists;
  vector<double*> double_matrices;
  double_matrices.push_back(nullptr);
  gaps.push_back(0.0f);
  gap_extends.push_back(0.0f);
  pam_dists.push_back(0.0f);
  for (const auto& env : envs_) {
    gaps.push_back(env.gap_open);
    gap_extends.push_back(env.gap_extend);
    pam_dists.push_back(env.pam_distance);
    double_matrices.push_back(env.matrix);
  }

  CreateDayMatrices(gaps, gap_extends, pam_dists, double_matrices);
}

void AlignmentEnvironments::UseBlosum(const json& blosum_json,
                                      double threshold) {
  std::fill(just_score_env_.matrix, just_score_env_.matrix + MDIM, 0.0f);
  const auto& column_order = blosum_json["column_order"];
  const auto& compact_matrix = blosum_json["scores"];
  just_score_env_.threshold = threshold;
  just_score_env_.gap_open = (double)blosum_json["gap_open"];
  just_score_env_.gap_extend = (double)blosum_json["gap_ext"];
  // cout << "blosum override just score threshold is " <<
  // just_score_env_.threshold << "\n";

  for (size_t i = 0; i < column_order.size(); i++) {
    if (compact_matrix[i].size() != column_order.size()) {
      cout << "column order doesnt match!! " << compact_matrix[i].dump()
           << " \n"
           << column_order << "\n";
      exit(0);
    }
    for (size_t j = 0; j < column_order.size(); j++) {
      // cout << "column order first char is " <<
      // column_order[i].get<string>().front() << "\n";
      just_score_env_
          .matrix[(column_order[i].get<string>().front() - 'A') * DIMSIZE +
                  column_order[j].get<string>().front() - 'A'] =
          (double)compact_matrix[i][j];
    }
  }

  just_score_env_.matrix_int16 = CreateScaled(
      just_score_env_.matrix, just_score_env_.threshold,
      just_score_env_.gap_open, just_score_env_.gap_extend,
      just_score_env_.gap_open_int16, just_score_env_.gap_ext_int16);
  just_score_env_.matrix_int8 =
      CreateScaled(just_score_env_.matrix, just_score_env_.threshold,
                   just_score_env_.gap_open, just_score_env_.gap_extend,
                   just_score_env_.gap_open_int8, just_score_env_.gap_ext_int8);
}