
#pragma once

#include "Parameters.h"
#include "Sampler_v1.h"

#include <vector>
#include <set>

//------------------------------------------------
// enumerate possible asexual and sexual inoculation status
enum Status_asexual {Inactive_asexual, Liverstage_asexual, Acute_asexual, Chronic_asexual};
enum Status_sexual {Inactive_sexual, Acute_sexual, Chronic_sexual};

//------------------------------------------------
// enumerate possible host status
enum Status_host {Host_Sh, Host_Eh, Host_Ah, Host_Ch};

//------------------------------------------------
// class defining human host
class Host : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // identifiers
  int index;      // where in the population (vector of hosts) this host resides
  int ID;         // unique ID, incremented upon death
  int home_deme;  // deme into which this host is born
  int deme;       // deme in which this host currently resides
  
  // pointers to deme counts, for modifying counts e.g. upon death
  std::vector<int>* Sh_ptr;
  std::vector<int>* Eh_ptr;
  std::vector<int>* Ah_ptr;
  std::vector<int>* Ch_ptr;
  
  // pointer to indices of infective hosts in each deme
  std::vector<std::vector<int>>* host_infective_index_ptr;
  
  // pointers to sampler objects, for efficiently drawing from global
  // distributions
  Sampler* sampler_age_stable_ptr;
  Sampler* sampler_age_death_ptr;
  std::vector<Sampler>* sampler_duration_acute_ptr;
  std::vector<Sampler>* sampler_duration_chronic_ptr;
  
  // pointers to scheduler objects, for adding events to schedulers
  std::vector<std::set<int>>* schedule_death_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Eh_to_Ah_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Eh_to_Ch_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Ah_to_Ch_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Ah_to_Sh_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Ch_to_Sh_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_acute_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_chronic_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_recovery_ptr;
  
  // cumulative count of how many times this host has been bitten by infective
  // mosquito (infection_index) and how many times an infection has taken hold
  // (inoc_index)
  int infection_index;
  int inoc_index;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // inoculation objects
  std::vector<bool> inoc_active;
  std::vector<Status_asexual> inoc_status_asexual;
  std::vector<Status_sexual> inoc_status_sexual;
  std::vector<int> inoc_time_infective;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // other methods
  void init(int index, int &ID, int deme,
            std::vector<int> &Sh, std::vector<int> &Eh, std::vector<int> &Ah, std::vector<int> &Ch,
            std::vector<std::vector<int>> &host_infective_index,
            std::vector<std::set<int>> &schedule_death,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Eh_to_Ah,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Eh_to_Ch,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Ah_to_Ch,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Ah_to_Sh,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Ch_to_Sh,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective_acute,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective_chronic,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective_recovery,
            Sampler &sampler_age_stable, Sampler &sampler_age_death,
            std::vector<Sampler> &sampler_duration_acute, std::vector<Sampler> &sampler_duration_chronic);
  
  void draw_starting_age();
  void death(int &ID, int birth_day);
  void denovo_infection(int t);
  void infection(int t);
  void Eh_to_Ah(int this_slot);
  void Eh_to_Ch(int this_slot);
  void Ah_to_Ch(int this_slot);
  void Ah_to_Sh(int this_slot);
  void Ch_to_Sh(int this_slot);
  void begin_infective_acute(int this_slot, int t);
  void begin_infective_chronic(int this_slot, int t);
  void end_infective(int this_slot);
  
  // getters and setters
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
  double get_infectivity(int t);
  int get_free_inoc_slot();
};
