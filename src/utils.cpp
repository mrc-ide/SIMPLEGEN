
#include "utils.hpp"

using namespace std;

//------------------------------------------------
// create vector of sampler objects from array of probabilities
vector<Sampler> make_sampler_vec(vector<vector<double>> &prob_array, int sampler_draws) {
  int n = prob_array.size();
  vector<Sampler> ret = vector<Sampler>(n);
  for (int i = 0; i < n; ++i) {
    ret[i] = Sampler(prob_array[i], sampler_draws);
  }
  return ret;
}
