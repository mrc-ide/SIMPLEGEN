
#pragma once

#include <cpp11.hpp>
#include <dust/random/random.hpp>
#include <iostream>
#include "../misc_v17.h"
#include "../probability.h"
#include "../Sampler_v6.h"
#include "Parameters.h"
#include "Dispatcher.h"

//------------------------------------------------
cpp11::list sim_indiv(cpp11::list args, cpp11::list args_progress);
