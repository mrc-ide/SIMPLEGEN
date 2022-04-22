
#pragma once

#include <vector>

#include "sim_Host.h"

//------------------------------------------------
// class defining mosquito
class sim_Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  int mosquito_ID;
  int infection_ID;
  
  int source_time;
  int source_deme;
  int source_host_ID;
  std::vector<int> source_infection_ID_vec;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  sim_Mosquito() {};
  
  // main events
  void infection(int t, int &next_infection_ID, sim_Host &host);
  void write_buffer(std::ofstream &transmission_record);
  void print_status();
  
  // getters and setters
  void set_mosquito_ID(int &mosquito_ID);
  
};
