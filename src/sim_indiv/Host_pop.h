
#pragma once

#include "Parameters.h"
#include "Host.h"
#include "Mosquito.h"
//#include "Mosquito_pop.h"
#include "../misc_v17.h"

class Mosquito_pop; // forwards declaration

//------------------------------------------------
// class defining population of hosts over all demes
class Host_pop {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointers to parameters and transmission record
  Parameters* params;
  std::ofstream* transmission_record;
  
  // total number of hosts in each deme
  std::vector<int> H;
  
  // vector of all human hosts
  std::vector<Host> host_vec;
  
  // the integer index of hosts in each deme
  std::vector<std::vector<int>> host_index;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host_pop() {};
  
  // main events
  void init(Parameters &params_, std::ofstream &transmission_record_);
  void seed_infections();
  void update_hosts(int t);
  void draw_new_infections(Mosquito_pop &mosq_pop, int t, int k);
};
