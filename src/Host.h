
#pragma once

#include "Mosquito.h"
#include "Parameters.h"
#include "Sampler_v2.h"

#include <vector>
#include <set>
#include <fstream>

//------------------------------------------------
// enumerate possible asexual and sexual inoculation status
enum Status_asexual {Inactive_asexual, Liverstage_asexual, Acute_asexual, Chronic_asexual};
enum Status_sexual {Inactive_sexual, Acute_sexual, Chronic_sexual};

//------------------------------------------------
// enumerate possible host status
enum Status_host {Host_Sh, Host_Eh, Host_Ah, Host_Ch, Host_Ph};

//------------------------------------------------
// enumerate possible events
// NB, the order of this enum is important as inoc_events will be explored in
// this order. For example, it is OK for Event_Eh_to_Ah and Event_Ah_to_Sh to
// occur at the same time point as long as they are this way round, but not the
// other way round
enum Event {Event_Eh_to_Ah, Event_Eh_to_Ch,
            Event_Ah_to_Ch, Event_Ah_to_Sh,
            Event_Ch_to_Sh,
            Event_begin_infective_acute, Event_begin_infective_chronic, Event_end_infective,
            Event_treatment};

//------------------------------------------------
// class defining human host
class Host : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // identifiers
  int index;      // where in the population (vector of hosts) this host resides. The index never changes
  int host_ID;    // unique ID, incremented upon death
  int home_deme;  // deme into which this host is born
  int deme;       // deme in which this host currently resides
  
  // pointer to indices of infective hosts in each deme
  std::vector<std::vector<int>>* host_infective_index_ptr;
  
  // pointers to sampler objects, for efficiently drawing from global
  // distributions
  Sampler* sampler_age_stable_ptr;
  Sampler* sampler_age_death_ptr;
  std::vector<Sampler>* sampler_duration_acute_ptr;
  std::vector<Sampler>* sampler_duration_chronic_ptr;
  std::vector<Sampler>* sampler_time_treatment_acute_ptr;
  std::vector<Sampler>* sampler_time_treatment_chronic_ptr;
  Sampler* sampler_duration_prophylatic_ptr;
  
  // cumulative count of how many times this host has been bitten by infective
  // mosquito (infection_index) and how many times an infection has taken hold
  // (inoc_index)
  int infection_index;
  int inoc_index;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // prophylactic status
  bool prophylaxis_on;
  int t_prophylaxis_stop;
  
  // host characteristics
  double treatment_seeking;
  
  // inoculation slots
  std::vector<int> inoc_ID_vec;
  std::vector<bool> inoc_active;
  std::vector<Status_asexual> inoc_status_asexual;
  std::vector<Status_sexual> inoc_status_sexual;
  std::vector<int> inoc_time_infective;
  
  // events
  std::vector<std::map<Event, int>> inoc_events;
  int t_next_inoc_event;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // other methods
  void init(int index, int &host_ID, int deme,
            std::vector<std::vector<int>> &host_infective_index,
            Sampler &sampler_age_stable,
            Sampler &sampler_age_death,
            std::vector<Sampler> &sampler_duration_acute,
            std::vector<Sampler> &sampler_duration_chronic,
            std::vector<Sampler> &sampler_time_treatment_acute,
            std::vector<Sampler> &sampler_time_treatment_chronic,
            Sampler &sampler_duration_prophylactic);
  void draw_starting_age();
  void draw_treatment_seeking();
  
  void new_inoc_event(int t, Event this_event, int this_slot);
  void check_inoc_event(int t);
  
  void death(int &host_ID, int t);
  void denovo_infection(int t, int &next_inoc_ID, std::ofstream &transmission_record);
  void infection(int t, int &next_inoc_ID, Mosquito &mosq, std::ofstream &transmission_record);
  void treatment(int t);
  void end_prophylaxis();
  
  void Eh_to_Ah(int this_slot);
  void Eh_to_Ch(int this_slot);
  void Ah_to_Ch(int this_slot);
  void Ah_to_Sh(int this_slot);
  void Ch_to_Sh(int this_slot);
  void begin_infective_acute(int this_slot, int t);
  void begin_infective_chronic(int this_slot, int t);
  void end_infective(int this_slot);
  
  // getters and setters
  std::vector<int> get_inoc_ID_vec();
  int get_n_active_inoc();
  Status_host get_host_status();
  int get_n_liverstage();
  int get_n_bloodstage_acute();
  int get_n_bloodstage_chronic();
  int get_n_bloodstage();
  int get_n_asexual();
  int get_n_infective();
  double get_prob_infection();
  double get_prob_acute();
  double get_prob_AC();
  int get_duration_acute();
  int get_duration_chronic();
  int get_time_treatment_acute();
  int get_time_treatment_chronic();
  int get_duration_prophylaxis();
  double get_infectivity(int t);
  int get_free_inoc_slot();
  int get_age(int t);
  
  // print methods
  void print_inoc_events();
  void print_inoc_status();
};
