
#pragma once

#include "misc_v6.h"
#include "Parameters.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// class for dispatching main simulation. Inherits parameters
class Dispatcher : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scheduler objects
  std::vector<std::set<int>> schedule_death;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher() {};
  
  // methods
  void run_simulation();
  
};

