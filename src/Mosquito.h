
#pragma once

#include <vector>

#include "Host.h"

//------------------------------------------------
// class defining mosquito
class Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  int mosquito_ID;
  int infection_ID;
  
  int source_time;
  int source_host_ID;
  std::vector<int> source_infection_ID_vec;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito() {};
  
  // main events
  void infection(int t, int &next_infection_ID, Host &host);
  void write_buffer(std::ofstream &transmission_record);
  
  // getters and setters
  void set_mosquito_ID(int &mosquito_ID);
  
};
