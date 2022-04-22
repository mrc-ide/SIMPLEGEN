
#pragma once

#include "misc_v14.h"
#include "sim_Parameters.h"
#include "Sampler_v5.h"
#include "sim_Host.h"
#include "sim_Mosquito.h"
#include "utils.h"

#include <set>
#include <fstream>

//------------------------------------------------
// class for dispatching main simulation. Inherits parameters
class sim_Dispatcher {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  sim_Parameters * params;
  
  // filestream to transmission record
  std::ofstream transmission_record;
  
  // used to define unique IDs for hosts, mosquitoes and infections
  int next_host_ID;
  int next_mosq_ID;
  int next_infection_ID;
  
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
  
  // vector of all human hosts
  std::vector<sim_Host> host_pop;
  
  // the integer index of hosts in each deme. host_infective_index is a subset
  // of host_index
  std::vector<std::vector<int>> host_index;
  std::vector<std::vector<int>> host_infective_index;
  
  // counts of mosquito types
  std::vector<int> Sv;
  std::vector<int> Ev;
  std::vector<int> Iv;
  
  // objects for tracking mosquitoes that die in lag phase
  std::vector<std::vector<int>> Ev_death;
  
  // arrays of mosquitoes at various stages
  std::vector<std::vector<std::vector<sim_Mosquito>>> Ev_pop;
  std::vector<std::vector<sim_Mosquito>> Iv_pop;
  
  // objects for storing results
  std::vector<std::vector<double>> daily_output;
  std::vector<double> sweep_output;
  Rcpp::List surveys_indlevel_output;
  Rcpp::List surveys_indlevel_output_infection_IDs;
  
  // misc
  std::vector<double> daily_EIR;
  std::vector<double> prob_infectious_bite;

  
  // PUBLIC FUNCTIONS
  
  // constructors
  sim_Dispatcher() {};
  
  // methods
  void init(sim_Parameters &params);
  void open_trans_record();
  void run_simulation(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void get_sample_details(int t, int deme, int n, Diagnostic diag);
  void get_age_distribution(int t_index);
  
};

