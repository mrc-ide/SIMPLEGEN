
#pragma once

#include "Parameters.hpp"
#include "Host_inoc.hpp"
#include "../misc_v17.hpp"

//------------------------------------------------
// enumerate possible host states
enum State_host {Host_Sh, Host_Eh, Host_Ah, Host_Ch, Host_Ph};

//------------------------------------------------
// class defining human host
class Host {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointers to parameters and transmission record
  Parameters* params;
  std::ofstream* transmission_record;
  
  // identifiers
  int host_ID;    // unique ID of this host, modified upon death
  int home_deme;  // deme into which this host was born
  int deme;       // deme in which this host currently resides
  
  // host status now, previously, and time at which it changed
  State_host host_state;
  State_host host_state_previous;
  int time_host_state_change;
  
  // dates of birth and death
  int time_birth;
  int time_death;
  
  // treatment status
  int time_treatment;
  bool prophylaxis_on;
  int time_prophylaxis_stop;
  
  // cumulative count of how many times this host has been infected
  int cumul_infections;
  
  // host characteristics
  double treatment_seeking_prob;  // host-specific treatment seeking probability
  
  // inoculation slots
  bool any_inoc_active;
  std::vector<Host_inoc> inoc;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // initialisation
  void init(Parameters &params_, std::ofstream &transmission_record_, int deme_);
  void draw_starting_age();
  void draw_treatment_seeking();
  
  // major life events
  void migrate(int new_deme);
  void death(int t);
  
  // update
  void update(int t);
  void treatment(int t);
  void end_prophylaxis(int t);
  
  // infection
  void infect(int t, int mosquito_ID, int mosquito_inoc_ID);
  
  // getters
  double get_prob_infection(int t);
  double get_infectivity(int t);
  double get_detectability_microscopy(int t);
  double get_detectability_PCR(int t);
  int get_age_years(int t);
  bool get_any_inoc_active();
  State_host get_host_state();
  State_host get_host_state_previous();
  int get_time_host_state_change();
  int get_time_treatment();
  
  // misc
  void recalc_host_state(int t);
  void print_summary();
};
