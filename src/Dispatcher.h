
#pragma once

#include "misc_v6.h"
#include "Parameters.h"
#include "Sampler_v1.h"

#include <set>

//------------------------------------------------
// class for dispatching main simulation. Inherits parameters
class Dispatcher : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  Sampler sampler_duration_infection;
  
  // scheduler objects
  std::vector<std::set<int>> schedule_death;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void run_simulation();
  
};

