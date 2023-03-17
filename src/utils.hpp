
#pragma once

#include "Sampler_v6.hpp"

//------------------------------------------------
// create vector of sampler objects from array of probabilities
std::vector<Sampler> make_sampler_vec(std::vector<std::vector<double>> &prob_array, int sampler_draws);
