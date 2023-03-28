
#pragma once

#include "Sampler_v6.h"

//------------------------------------------------
// create vector of sampler objects from array of probabilities
std::vector<Sampler> make_sampler_vec(std::vector<std::vector<double>> &prob_array, int sampler_draws);

//------------------------------------------------
// dumb function for making manual copy of std matrix into cpp11 matrix
cpp11::writable::doubles_matrix<> copy_mat(std::vector<std::vector<double>> &c_mat);
cpp11::writable::integers_matrix<> copy_mat(std::vector<std::vector<int>> &c_mat);
