
#pragma once

#include "Parameters.h"
#include "Host.h"
#include "../misc_v17.h"

//------------------------------------------------
// class defining mosquito
class Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  
  int mosq_ID;
  int inoc_ID;
  
  bool infectious_on;
  int time_infectious;
  int time_death;
  
  int source_time;
  int source_deme;
  int source_host_ID;
  std::vector<int> source_inoc_ID;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito();
  
  // main events
  void init_from_host(Parameters &params, int t, int t_death, Host &host);
  void init_from_mosq(Mosquito &mosq);
  void Ev_to_Iv();
  void write_buffer(std::ofstream &transmission_record);
  void print_summary();
};
