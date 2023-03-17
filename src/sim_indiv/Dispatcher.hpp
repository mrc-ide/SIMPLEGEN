
#pragma once

//#include <set>
//#include <fstream>
#include <cpp11.hpp>
#include "Parameters.hpp"
//#include "misc_v14.h"
//#include "Sampler_v5.h"
//#include "sim_Host.h"
//#include "sim_Mosquito.h"
//#include "utils.h"

//------------------------------------------------
// class for dispatching main simulation
class Dispatcher {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters* params;
  
  // filestream to transmission record
  std::ofstream transmission_record;
  
  // used to define unique IDs for hosts, mosquitoes and inoculations
  int next_host_ID;
  int next_mosq_ID;
  int next_inoculation_ID;
  
  // counts of host types in each deme
  std::vector<int> H;
  std::vector<int> Sh;
  std::vector<int> Eh;
  std::vector<int> Ah;
  std::vector<int> Ch;
  std::vector<int> Ph;
  /*
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
  */
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher() {};
  
  // methods
  void init(Parameters &params);
  void open_trans_record();
  void run_simulation(cpp11::list &args_progress);
  //void get_sample_details(int t, int deme, int n, Diagnostic diag);
  //void get_age_distribution(int t_index);
  
};

