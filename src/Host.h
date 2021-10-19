
#pragma once

#include "Parameters.h"
#include "Sampler_v4.h"

#include <vector>
#include <set>
#include <fstream>

//------------------------------------------------
// enumerate possible asexual and sexual infection states
enum State_asexual {Inactive_asexual, Liverstage_asexual, Acute_asexual, Chronic_asexual};
enum State_sexual {Inactive_sexual, Acute_sexual, Chronic_sexual};

//------------------------------------------------
// enumerate possible host states
enum State_host {Host_Sh, Host_Eh, Host_Ah, Host_Ch, Host_Ph};

//------------------------------------------------
// enumerate possible events
enum Event {Event_Eh_to_Ah, Event_Eh_to_Ch,
            Event_Ah_to_Ch, Event_Ah_to_Sh,
            Event_Ch_to_Sh,
            Event_begin_infective_acute, Event_begin_infective_chronic, Event_end_infective,
            Event_treatment};

//------------------------------------------------
// class defining human host
class Host {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters * params;
  
  // identifiers
  int index;      // where in the population (vector of hosts) this host resides. This index never changes
  int host_ID;    // unique ID, modified upon death
  int home_deme;  // deme into which this host was born
  int deme;       // deme in which this host currently resides
  
  // pointers to indices of hosts and infective hosts in each deme
  std::vector<std::vector<int>>* host_index_ptr;
  std::vector<std::vector<int>>* host_infective_index_ptr;
  
  // pointers to sampler objects. Used to make efficient random draws from
  // global distributions
  Sampler* sampler_age_stable_ptr;
  Sampler* sampler_age_death_ptr;
  std::vector<Sampler>* sampler_duration_acute_ptr;
  std::vector<Sampler>* sampler_duration_chronic_ptr;
  std::vector<Sampler>* sampler_time_treatment_acute_ptr;
  std::vector<Sampler>* sampler_time_treatment_chronic_ptr;
  std::vector<Sampler>* sampler_duration_prophylatic_ptr;
  
  // cumulative count of how many times this host has been bitten by infective
  // mosquito (cumul_infections) and how many times an infection has taken hold
  // (cumul_inoculations)
  int cumul_infections;
  int cumul_inoculations;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // prophylactic status
  bool prophylaxis_on;
  int t_prophylaxis_stop;
  
  // host status now, previously, and time at which it changed
  State_host host_state;
  State_host host_state_previous;
  int time_host_state_change;
  
  // host characteristics
  double treatment_seeking;
  
  // inoculation slots
  std::vector<int> inoc_ID_vec;
  std::vector<bool> inoc_active;
  std::vector<State_asexual> inoc_state_asexual;
  std::vector<State_sexual> inoc_state_sexual;
  std::vector<int> inoc_time_asexual;
  std::vector<int> inoc_time_sexual;
  
  // events
  std::vector<std::map<Event, int>> inoc_events;
  int t_next_inoc_event;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // initialisation
  void init(Parameters &params_,
            int index, int &host_ID, int deme,
            std::vector<std::vector<int>> &host_index,
            std::vector<std::vector<int>> &host_infective_index,
            Sampler &sampler_age_stable,
            Sampler &sampler_age_death,
            std::vector<Sampler> &sampler_duration_acute,
            std::vector<Sampler> &sampler_duration_chronic,
            std::vector<Sampler> &sampler_time_treatment_acute,
            std::vector<Sampler> &sampler_time_treatment_chronic,
            std::vector<Sampler> &sampler_duration_prophylactic);
  void draw_starting_age();
  void draw_treatment_seeking();
  
  // event schedulers
  void new_inoc_event(int t, Event this_event, int this_slot);
  void check_inoc_event(int t);
  
  // host-level events
  void death(int &host_ID, int t);
  void migrate(int new_deme);
  void denovo_infection(int t, int &next_inoc_ID, std::ofstream &transmission_record);
  void infection(int t, int &next_inoc_ID, int mosquito_ID, int mosquito_infection_ID, std::ofstream &transmission_record);
  void treatment(int t);
  void end_prophylaxis();
  
  // inoculation-level events
  void new_Eh(int this_slot, int t, int &next_inoc_ID);
  void Eh_to_Ah(int this_slot, int t);
  void Eh_to_Ch(int this_slot, int t);
  void Ah_to_Ch(int this_slot, int t);
  void Ah_to_Sh(int this_slot, int t);
  void Ch_to_Sh(int this_slot, int t);
  void begin_infective_acute(int this_slot, int t);
  void begin_infective_chronic(int this_slot, int t);
  void end_infective(int this_slot);
  void recalculate_host_state(int t);
  
  // getters, setters and random draws
  int get_free_inoc_slot();
  int get_n_active_inoc();
  std::vector<int> get_inoc_ID_vec();
  State_host get_host_state();
  int get_n_liverstage();
  int get_n_bloodstage_acute();
  int get_n_bloodstage_chronic();
  int get_n_bloodstage();
  int get_n_asexual();
  int get_n_infective();
  double get_prob_infection();
  double get_prob_acute();
  double get_prob_AC();
  int draw_duration_acute();
  int draw_duration_chronic();
  int draw_time_treatment_acute();
  int draw_time_treatment_chronic();
  int draw_duration_prophylaxis();
  double get_detectability_microscopy_acute(int t);
  double get_detectability_microscopy_chronic(int t);
  double get_detectability_PCR_acute(int t);
  double get_detectability_PCR_chronic(int t);
  double get_infectivity(int t);
  int get_age(int t);
  
  // outputs
  void update_output(Measure measure,
                     Model_state state,
                     Diagnostic diagnostic,
                     int age_min,
                     int age_max,
                     int t,
                     double &numer,
                     double &denom);
  
  // diagnostics and checks
  void print_inoc_events();
  void print_inoc_state();
  void check_host_infective_index(int x);
  
};
