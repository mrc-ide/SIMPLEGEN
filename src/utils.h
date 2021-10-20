
#pragma once

#include "Sampler_v4.h"

//------------------------------------------------
// create vector of sampler objects from array of probabilities
std::vector<Sampler> make_sampler_vec(std::vector<std::vector<double>> prob_array, int sampler_draws);
