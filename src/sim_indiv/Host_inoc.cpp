
#include "Host_inoc.hpp"

using namespace std;


//------------------------------------------------
// constructor
Host_inoc::Host_inoc() {
  reset();
}

//------------------------------------------------
// initialise or reset inoculation
void Host_inoc::reset() {
  ID = -1;
  active = false;
  state_asexual = Inactive_asexual;
  state_sexual = Inactive_sexual;
  transition_AC = false;
  time_start_acute = -1;
  time_stop_acute = -1;
  time_start_chronic = -1;
  time_stop_chronic = -1;
}

//------------------------------------------------
// print summary
void Host_inoc::print_summary() {
  print("ID:", ID);
  print("active:", active);
  print("state_asexual:", state_asexual);
  print("state_sexual:", state_sexual);
  print("transition_AC:", transition_AC);
  print("time_start_acute:", time_start_acute);
  print("time_stop_acute:", time_stop_acute);
  print("time_start_chronic:", time_start_chronic);
  print("time_stop_chronic:", time_stop_chronic);
}
