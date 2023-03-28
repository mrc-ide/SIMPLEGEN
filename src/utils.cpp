
#include "utils.h"

using namespace std;

//------------------------------------------------
vector<Sampler> make_sampler_vec(vector<vector<double>> &prob_array, int sampler_draws) {
  int n = prob_array.size();
  vector<Sampler> ret = vector<Sampler>(n);
  for (int i = 0; i < n; ++i) {
    ret[i] = Sampler(prob_array[i], sampler_draws);
  }
  return ret;
}

//------------------------------------------------
cpp11::writable::doubles_matrix<> copy_mat(vector<vector<double>> &c_mat) {
  int n_row = int(c_mat.size());
  int n_col = int(c_mat[0].size());
  cpp11::writable::doubles_matrix<> ret(n_row, n_col);
  for (int i = 0; i < n_row; ++i) {
    for (int j = 0; j < n_col; ++j) {
      ret(i, j) = c_mat[i][j];
    }
  }
  return ret;
}
cpp11::writable::integers_matrix<> copy_mat(vector<vector<int>> &c_mat) {
  int n_row = int(c_mat.size());
  int n_col = int(c_mat[0].size());
  cpp11::writable::integers_matrix<> ret(n_row, n_col);
  for (int i = 0; i < n_row; ++i) {
    for (int j = 0; j < n_col; ++j) {
      ret(i, j) = c_mat[i][j];
    }
  }
  return ret;
}
