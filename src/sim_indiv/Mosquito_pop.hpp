
#pragma once

#include "Parameters.hpp"
#include "Host.hpp"
#include "Mosquito.hpp"
#include "../misc_v17.hpp"
//#include "Host_pop.hpp"

class Host_pop; // forward declaration

//------------------------------------------------
// class defining population of mosquitoes in a single deme
class Mosquito_pop {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters* params;
  
  // counts of mosquito types
  int M;
  int Sv;
  int Ev;
  int Iv;
  
  // ring buffer of deaths in latent stage
  int ringbuffer_index;
  std::vector<int> ringbuffer;
  
  // vector of mosquitoes
  int memory_increment;
  int mosq_vec_n;
  std::vector<Mosquito> mosq_vec;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito_pop() {};
  
  // main events
  void init(Parameters &params_, int M_);
  void apply_ringbuffer_death();
  void update_pop(int t);
  void draw_new_infections(Host_pop &host_pop, int t, int k);
  void print_summary();
};
