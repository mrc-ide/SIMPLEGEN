
#pragma once

#include "misc_v10.h"
#include "Parameters.h"
#include "Sampler_v4.h"
#include "Host.h"
#include "Mosquito.h"

#include <set>
#include <fstream>

//------------------------------------------------
// class for dispatching main simulation. Inherits parameters
class Dispatcher {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters * params;
  
  // filestream to transmission record
  std::ofstream transmission_record;
  
  // unique IDs for each inoculation
  int next_inoc_ID;
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  std::vector<Sampler> sampler_duration_acute;
  std::vector<Sampler> sampler_duration_chronic;
  std::vector<Sampler> sampler_time_treatment_acute;
  std::vector<Sampler> sampler_time_treatment_chronic;
  std::vector<Sampler> sampler_duration_prophylactic;
  
  // counts of host types
  std::vector<int> H;
  std::vector<int> Sh;
  std::vector<int> Eh;
  std::vector<int> Ah;
  std::vector<int> Ch;
  std::vector<int> Ph;
  
  // further counts of host types
  std::vector<double> Ah_detectable_microscopy;
  std::vector<double> Ch_detectable_microscopy;
  std::vector<double> Ah_detectable_PCR;
  std::vector<double> Ch_detectable_PCR;
  
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
  
  // objects for storing results
  std::vector<std::vector<double>> daily_numer;
  std::vector<std::vector<double>> daily_denom;
  std::vector<double> sweep_numer;
  std::vector<double> sweep_denom;
  
  // misc
  std::vector<double> daily_EIR;
  std::vector<double> prob_infectious_bite;

  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher() {};
  
  // methods
  void init(Parameters &params);
  void run_simulation(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void get_sample_details(int t, int deme, int n, Diagnostic diag);
  void get_age_distribution(int t_index);
  
};

