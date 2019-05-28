
#pragma once

#include "Parameters.h"
#include "Sampler_v1.h"

#include <vector>
#include <set>

//------------------------------------------------
// enumerate possible asexual and sexual innoculation status
enum Status_asexual {Inactive_asexual, Liverstage_asexual, Bloodstage_asexual};
enum Status_sexual {Inactive_sexual, Active_sexual};

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
  std::vector<int>* Ih_ptr;
  
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
  std::vector<std::vector<std::pair<int, int>>>* schedule_Eh_to_Ih_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Ih_to_Sh_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_recovery_ptr;
  
  // indices relating to global distributions. For example, the probability of
  // this host becoming infected is equal to prob_infection[prob_infection_index]
  int prob_infection_index;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // innoculation objects
  std::vector<bool> innoc_active;
  std::vector<Status_asexual> innoc_status_asexual;
  std::vector<Status_sexual> innoc_status_sexual;
  std::vector<int> innoc_time_infective;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // other methods
  void init(int index, int &ID, int deme,
            std::vector<int> &Sh, std::vector<int> &Eh, std::vector<int> &Ih,
            std::vector<std::vector<int>> &host_infective_index,
            std::vector<std::set<int>> &schedule_death,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Eh_to_Ih,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Ih_to_Sh,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective_recovery,
            Sampler &sampler_age_stable, Sampler &sampler_age_death,
            std::vector<Sampler> &sampler_duration_acute, std::vector<Sampler> &sampler_duration_chronic);
  
  void draw_starting_age();
  void death(int &ID, int birth_day);
  void denovo_infection(int t);
  void infection(int t);
  void Eh_to_Ih(int this_slot);
  void Ih_to_Sh(int this_slot);
  void begin_infective(int this_slot, int t);
  void end_infective(int this_slot);
  void update_prob_infection();
  
  // getters and setters
  int get_n_innoculations();
  int get_n_liverstage();
  int get_n_bloodstage();
  int get_n_asexual();
  int get_n_infective();
  double get_prob_infection();
  double get_infectivity(int t);
  int get_free_innoc_slot();
};
