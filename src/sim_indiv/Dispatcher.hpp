
#pragma once

#include <cpp11.hpp>
#include "Parameters.hpp"
#include "Host.hpp"
#include "Mosquito_pop.hpp"
#include "Host_pop.hpp"
#include "Sweep.hpp"
#include "../utils.hpp"

//------------------------------------------------
// class for dispatching main simulation
class Dispatcher {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters* params;
  
  // filestream to transmission record
  std::ofstream transmission_record;
  
  // population of human hosts over all demes
  Host_pop host_pop;
  
  // population of mosquitoes per deme
  std::vector<Mosquito_pop> mosq_pop;
  
  // objects for sampling from human and mosquito populations in either "daily"
  // or "sweep" formats
  Sweep daily;
  Sweep sweeps;
  
  // output
  std::vector<std::vector<double>> daily_output;
  cpp11::writable::list tmp_output;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher() {};
  
  // methods
  void init(Parameters &params_);
  void open_trans_record();
  void run_simulation(cpp11::list &args_progress);
};

