
#pragma once

#include "misc_v6.h"
#include "Parameters.h"
#include "Sampler_v1.h"
#include "Host.h"
#include "Mosquito.h"

#include <set>

//------------------------------------------------
// class for dispatching main simulation. Inherits parameters
class Dispatcher : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  Sampler sampler_duration_acute;
  
  // scheduler objects
  std::vector<std::set<int>> schedule_death;
  std::vector<std::vector<std::pair<int, int>>> schedule_Eh_to_Ih;
  std::vector<std::vector<std::pair<int, int>>> schedule_Ih_to_Sh;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_recovery;
  
  // counts of host types
  std::vector<int> H;
  std::vector<int> Sh;
  std::vector<int> Eh;
  std::vector<int> Ih;
  
  // population of human hosts
  std::vector<Host> host_pop;
  int next_host_ID;
  
  // store the integer index of hosts in each deme
  std::vector<std::vector<int>> host_index;
  std::vector<std::vector<int>> host_infective_index;
  
  // counts of mosquito types
  int M_total;
  std::vector<int> Sv;
  std::vector<int> Ev;
  std::vector<int> Iv;
  
  // objects for tracking mosquitoes that die in lag phase
  std::vector<std::vector<int>> Ev_death;
  
  // populations of mosquitoes at various stages
  std::vector<std::vector<std::vector<Mosquito>>> Ev_pop;
  std::vector<std::vector<Mosquito>> Iv_pop;
  
  // objects for storing daily values
  std::vector<std::vector<std::vector<double>>> daily_values;
  
  // misc
  std::vector<double> EIR;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void run_simulation();
  
};

