
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
  std::vector<std::vector<std::pair<int, int>>> schedule_Eh_to_Ih;
  std::vector<std::vector<std::pair<int, int>>> schedule_Ih_to_Sh;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_recovery;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void run_simulation();
  
};

