
#include "Host.h"
#include "misc_v6.h"
#include "probability_v2.h"

using namespace std;

//------------------------------------------------
// initialise host
void Host::init(int index, int &ID, int deme,
                vector<int> &Sh, vector<int> &Eh, vector<int> &Ih,
                vector<vector<int>> &host_infective_index,
                vector<set<int>> &schedule_death,
                vector<vector<pair<int, int>>> &schedule_Eh_to_Ih,
                vector<vector<pair<int, int>>> &schedule_Ih_to_Sh,
                vector<vector<pair<int, int>>> &schedule_infective,
                vector<vector<pair<int, int>>> &schedule_infective_recovery,
                Sampler &sampler_age_stable, Sampler &sampler_age_death,
                vector<Sampler> &sampler_duration_acute, vector<Sampler> &sampler_duration_chronic) {
  
  // identifiers
  this->index = index;
  this->ID = ID++;
  home_deme = deme;
  this->deme = deme;
  
  // pointers to deme counts, for modifying counts e.g. upon death
  Sh_ptr = &Sh;
  Eh_ptr = &Eh;
  Ih_ptr = &Ih;
  
  // pointer to indices of infective hosts in each deme
  host_infective_index_ptr = &host_infective_index;
  
  // pointers to sampler objects, for efficiently drawing from global
  // distributions
  sampler_age_stable_ptr = &sampler_age_stable;
  sampler_age_death_ptr = &sampler_age_death;
  sampler_duration_acute_ptr = &sampler_duration_acute;
  sampler_duration_chronic_ptr = &sampler_duration_chronic;
  
  // pointers to scheduler objects, for adding events to schedulers
  schedule_death_ptr = &schedule_death;
  schedule_Eh_to_Ih_ptr = &schedule_Eh_to_Ih;
  schedule_Ih_to_Sh_ptr = &schedule_Ih_to_Sh;
  schedule_infective_ptr = &schedule_infective;
  schedule_infective_recovery_ptr = &schedule_infective_recovery;
  
  // indices relating to global distributions. For example, the probability of
  // this host becoming infected is equal to prob_infection[prob_infection_index]
  prob_infection_index = 0;
  
  // draw age from stable demography distribution
  draw_starting_age();
  
  // initialise innoculation objects
  innoc_active = vector<bool>(max_innoculations, false);
  innoc_status_asexual = vector<Status_asexual>(max_innoculations, Inactive_asexual);
  innoc_status_sexual = vector<Status_sexual>(max_innoculations, Inactive_sexual);
  innoc_time_infective = vector<int>(max_innoculations);
  
}

//------------------------------------------------
// draw birth and death days given a known stable demography distribution, and
// schedule death for future. Only needed at start of simulation, after which
// hosts that die are instantly re-born and have death drawn from entire life
// table.
void Host::draw_starting_age() {
  
  // draw age from stable demography distribution
  int age_years = sampler_age_stable_ptr->draw();
  int extra_days = sample2(0, 364);
  int age_days = age_years*365 + extra_days;
  
  // draw duration of life from demography distribution looking forward from
  // current age. This is tricky, as we must account for the fact that if we
  // are already part way into an age group we have a reduced probability of
  // dying within that age group.
  //
  // math derivation:
  // assume constant rate of death r over 1-year period. Total probability of
  // dying this year is therefore p = 1 - exp(-r). We know p from life table,
  // from which we can derive r = -log(1-p). If we are already a proportion x
  // through this year, then the probability of dying in the remaining time is
  // Pr(die) = 1 - exp(-r(1-x)). Sustituting in r and simplifying we get
  // Pr(die) = 1 - (1-p)exp(1-x). Notice that if p=1 from the life table then we
  // are certain to die this year, irrespective of how far through the year we
  // have got already.
  int life_days = 0;
  double prop_year_remaining = 1.0 - extra_days/365.0; // (x in the above derivation)
  double prob_die_this_year = 1.0 - (1.0 - life_table[age_years])*exp(1.0 - prop_year_remaining);
  if (rbernoulli1(prob_die_this_year)) {
    life_days = age_years*365 + sample2(extra_days, 364);
  } else {
    // if we do not die in year age_years then loop through all remaining years,
    // performing Bernoulli draw from probability of dying in that year
    for (int i = (age_years+1); i < int(life_table.size()); ++i) {
      if (rbernoulli1(life_table[i])) {
        life_days = i*365 + sample2(0, 364);
        break;
      }
    }
  }
  
  // calculate final birth and death days from age at time 0 and life_days
  birth_day = -age_days;
  death_day = life_days - age_days;
  
  // add death_day to scheduler
  if (death_day < max_time) {
    (*schedule_death_ptr)[death_day].insert(index);
  }
  
}

//------------------------------------------------
// death
void Host::death(int &ID, int birth_day) {
  
  // drop from infective list if necessary
  if (get_n_infective() > 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
  // update deme counts to reflect death
  if (get_n_bloodstage() > 0) {
    (*Ih_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  } else if (get_n_liverstage() > 0) {
    (*Eh_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  }
  
  // new unique ID
  this->ID = ID++;
  
  // make current deme home deme
  home_deme = deme;
  
  // reset indices relating to global distributions
  prob_infection_index = 0;
  
  // set date of birth
  this->birth_day = birth_day;
  
  // draw life duration from demography distribution
  int life_years = sampler_age_death_ptr->draw();
  int life_days = life_years*365 + sample2(0, 364);
  if (life_days == 0) {  // cannot die on same day born
    life_days = 1;
  }
  death_day = birth_day + life_days;
  
  // add new death_day to scheduler
  if (death_day < max_time) {
    (*schedule_death_ptr)[death_day].insert(index);
  }
  
  // reset innoculation objects
  fill(innoc_active.begin(), innoc_active.end(), false);
  fill(innoc_status_asexual.begin(), innoc_status_asexual.end(), Inactive_asexual);
  fill(innoc_status_sexual.begin(), innoc_status_sexual.end(), Inactive_sexual);
  fill(innoc_time_infective.begin(), innoc_time_infective.end(), 0);
  
}

//------------------------------------------------
// de-novo infection
void Host::denovo_infection(int t) {
  
  infection(t);
  
  /*
  // generating starting genotype in a dummy mosquito
  Mosquito dummy_mosquito;
  dummy_mosquito.init(param_ptr);
  dummy_mosquito.denovo_infection();
  
  // carry out infection
  new_infection(dummy_mosquito, 0);
  */
}

//------------------------------------------------
// new infection
//void Host::new_infection(Mosquito &mosq, int t) {
void Host::infection(int t) {
  
  // update prob_infection_index irrespective of whether infection takes hold
  update_prob_infection();
  
  // return if already at max_innoculations
  if (get_n_innoculations() == max_innoculations) {
    return;
  }
  
  // update deme counts
  if (get_n_asexual() == 0) {
    (*Sh_ptr)[deme]--;
    (*Eh_ptr)[deme]++;
  }
  
  // get next free innoculation slot
  int this_slot = get_free_innoc_slot();
  
  // add new innoculation
  innoc_active[this_slot] = true;
  innoc_status_asexual[this_slot] = Liverstage_asexual;
  
  // draw duration of infection
  int duration_infection = (*sampler_duration_acute_ptr)[0].draw() + 1;
  
  // get times of future events
  int t1 = t + u;                    // begin bloodstage
  int t2 = t1 + duration_infection;  // end bloodstage
  int t3 = t1 + g;                   // begin infective
  int t4 = t3 + duration_infection;  // end infective
  
  // schedule move to Ih
  if (t1 < death_day && t1 < max_time) {
    (*schedule_Eh_to_Ih_ptr)[t1].emplace_back(index, this_slot);
  }
  
  // schedule bloodstage recovery
  if (t2 < death_day && t2 < max_time) {
    (*schedule_Ih_to_Sh_ptr)[t2].emplace_back(index, this_slot);
  }
  
  // schedule begin infective
  if (t3 < death_day && t3 < max_time) {
    (*schedule_infective_ptr)[t3].emplace_back(index, this_slot);
  }
  
  // schedule end infective
  if (t4 < death_day && t4 < max_time) {
    (*schedule_infective_recovery_ptr)[t4].emplace_back(index, this_slot);
  }
  
}

//------------------------------------------------
// move innoculation from Eh state to Ih
void Host::Eh_to_Ih(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage() == 0) {
    (*Eh_ptr)[deme]--;
    (*Ih_ptr)[deme]++;
  }
  
  // update status
  innoc_status_asexual[this_slot] = Bloodstage_asexual;
  
}

//------------------------------------------------
// move innoculation from Ih state to Sh
void Host::Ih_to_Sh(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage() == 1) {
    (*Ih_ptr)[deme]--;
    if (get_n_liverstage() == 0) {
      (*Sh_ptr)[deme]++;
    } else {
      (*Eh_ptr)[deme]++;
    }
  }
  
  // update status
  innoc_status_asexual[this_slot] = Inactive_asexual;
  
}

//------------------------------------------------
// begin infective period
void Host::begin_infective(int this_slot, int t) {
  
  // update host status
  innoc_status_sexual[this_slot] = Active_sexual;
  
  // store time
  innoc_time_infective[this_slot] = t;
  
  // if newly infective then add to infectives list
  if (get_n_infective() == 1) {
    (*host_infective_index_ptr)[deme].push_back(index);
  }
  
}

//------------------------------------------------
// end infective period
void Host::end_infective(int this_slot) {
  
  // update host status
  innoc_status_sexual[this_slot] = Inactive_sexual;
  innoc_active[this_slot] = false;
  
  // if no longer infective then drop from infectives list
  if (get_n_infective() == 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
}

//------------------------------------------------
// update probabilty of infection
void Host::update_prob_infection() {
  if (prob_infection_index < (n_prob_infection-1)) {
    prob_infection_index++;
  }
}

//------------------------------------------------
// get total number of innoculations
int Host::get_n_innoculations() {
  return sum_bool(innoc_active);
}

//------------------------------------------------
// get total number of liverstage innoculations
int Host::get_n_liverstage() {
  int ret = 0;
  for (int i = 0; i < int(innoc_status_asexual.size()); ++i) {
    if (innoc_status_asexual[i] == Liverstage_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of bloodstage innoculations
int Host::get_n_bloodstage() {
  int ret = 0;
  for (int i = 0; i < int(innoc_status_asexual.size()); ++i) {
    if (innoc_status_asexual[i] == Bloodstage_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of asexual stage innoculations
int Host::get_n_asexual() {
  return get_n_liverstage() + get_n_bloodstage();
}

//------------------------------------------------
// get total number of infective (sexual stage) innoculations
int Host::get_n_infective() {
  int ret = 0;
  for (int i = 0; i < int(innoc_status_sexual.size()); ++i) {
    if (innoc_status_sexual[i] == Active_sexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  return prob_infection[prob_infection_index];
}

//------------------------------------------------
// get current infectivity
double Host::get_infectivity(int t) {
  
  // TODO - infectivity from distribution
  /*
  // get probability of being infected by sexual-stage innoculation
  int ret = 0;
  for (int i = 0; i < int(innoc_status_sexual.size()); ++i) {
    if (innoc_status_sexual[i] == Active_sexual) {
      int time_diff = t - innoc_time_infective[i];
      //infectivity_acute
    }
  }
  */
  
  return max_infectivity;
}

//------------------------------------------------
// get next free innoculation slot
int Host::get_free_innoc_slot() {
  
  // loop until find inactive slot
  int ret = 0;
  for (int i = 0; i < max_innoculations; ++i) {
    if (!innoc_active[i]) {
      break;
    }
    ret++;
  }
  
  // error check for outside range
  if (ret == max_innoculations) {
    Rcpp::stop("could not find free innoculation slot");
  }
  
  return ret;
}
