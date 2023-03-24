
#pragma once

#include "Parameters.hpp"
#include "../misc_v17.hpp"

//------------------------------------------------
// enumerate possible asexual and sexual infection states
enum State_asexual {Inactive_asexual, Liverstage_asexual, Acute_asexual, Chronic_asexual};
enum State_sexual {Inactive_sexual, Acute_sexual, Chronic_sexual};

//------------------------------------------------
// leightweght class defining a single inoculation inside a human host
class Host_inoc {
  
public:
  
  // PUBLIC OBJECTS
  
  int ID;
  bool active;
  
  State_asexual state_asexual;
  State_sexual state_sexual;
  bool transition_AC;
  int time_start_acute;
  int time_stop_acute;
  int time_start_chronic;
  int time_stop_chronic;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host_inoc();
  
  // main methods
  void reset();
  void print_summary();
  
};
