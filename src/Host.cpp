
#include "Host.h"
#include "misc_v6.h"
#include "probability_v2.h"

using namespace std;

//------------------------------------------------
// initialise host
void Host::init(int index, int &ID, int deme,
                vector<int> &Sh, vector<int> &Eh, vector<int> &Ah, vector<int> &Ch,
                vector<vector<int>> &host_infective_index,
                vector<set<int>> &schedule_death,
                vector<vector<pair<int, int>>> &schedule_Eh_to_Ah,
                vector<vector<pair<int, int>>> &schedule_Eh_to_Ch,
                vector<vector<pair<int, int>>> &schedule_Ah_to_Ch,
                vector<vector<pair<int, int>>> &schedule_Ah_to_Sh,
                vector<vector<pair<int, int>>> &schedule_Ch_to_Sh,
                vector<vector<pair<int, int>>> &schedule_infective_acute,
                vector<vector<pair<int, int>>> &schedule_infective_chronic,
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
  Ah_ptr = &Ah;
  Ch_ptr = &Ch;
  
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
  schedule_Eh_to_Ah_ptr = &schedule_Eh_to_Ah;
  schedule_Eh_to_Ch_ptr = &schedule_Eh_to_Ch;
  schedule_Ah_to_Ch_ptr = &schedule_Ah_to_Ch;
  schedule_Ah_to_Sh_ptr = &schedule_Ah_to_Sh;
  schedule_Ch_to_Sh_ptr = &schedule_Ch_to_Sh;
  schedule_infective_acute_ptr = &schedule_infective_acute;
  schedule_infective_chronic_ptr = &schedule_infective_chronic;
  schedule_infective_recovery_ptr = &schedule_infective_recovery;
  
  // cumulative count of how many times this host has been bitten by infective
  // mosquito (infection_index) and how many times an infection has taken hold
  // (inoc_index)
  infection_index = 0;
  inoc_index = 0;
  
  // draw age from stable demography distribution
  draw_starting_age();
  
  // initialise inoculation objects
  inoc_active = vector<bool>(max_inoculations, false);
  inoc_status_asexual = vector<Status_asexual>(max_inoculations, Inactive_asexual);
  inoc_status_sexual = vector<Status_sexual>(max_inoculations, Inactive_sexual);
  inoc_time_infective = vector<int>(max_inoculations);
  
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
  if (death_day == 0) {  // cannot die on day 0
    death_day = 1;
  }
  
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
  if (get_n_bloodstage_acute() > 0) {
    (*Ah_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  } else if (get_n_bloodstage_chronic() > 0) {
    (*Ch_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  } else if (get_n_liverstage() > 0) {
    (*Eh_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  }
  
  // new unique ID
  this->ID = ID++;
  
  // make current deme home deme
  home_deme = deme;
  
  // reset cumulative counts
  infection_index = 0;
  inoc_index = 0;
  
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
  
  // reset inoculation objects
  fill(inoc_active.begin(), inoc_active.end(), false);
  fill(inoc_status_asexual.begin(), inoc_status_asexual.end(), Inactive_asexual);
  fill(inoc_status_sexual.begin(), inoc_status_sexual.end(), Inactive_sexual);
  fill(inoc_time_infective.begin(), inoc_time_infective.end(), 0);
  
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
  
  // update infection_index irrespective of whether exit due to reaching max
  // inoculations
  infection_index++;
  
  // return if already at max_inoculations
  if (get_n_active_inoc() == max_inoculations) {
    return;
  }
  
  // update deme counts
  if (get_n_asexual() == 0) {
    (*Sh_ptr)[deme]--;
    (*Eh_ptr)[deme]++;
  }
  
  // get next free inoculation slot
  int this_slot = get_free_inoc_slot();
  
  // add new inoculation
  inoc_active[this_slot] = true;
  inoc_status_asexual[this_slot] = Liverstage_asexual;
  
  // schedule future events
  int t1 = t + u;
  if (t1 < (death_day-1) && t1 < max_time) {
    
    bool acute = rbernoulli1(get_prob_acute());
    if (acute) {
      
      (*schedule_Eh_to_Ah_ptr)[t1].emplace_back(index, this_slot);
      int t2 = t1 + g;
      if (t2 < (death_day-1) && t2 < max_time) {
        (*schedule_infective_acute_ptr)[t2].emplace_back(index, this_slot);
      }
      
      int duration_acute = get_duration_acute();
      int t3 = t1 + duration_acute;
      if (t3 < (death_day-1) && t3 < max_time) {
        
        bool acute_to_chronic = rbernoulli1(get_prob_AC());
        if (acute_to_chronic) {
          
          (*schedule_Ah_to_Ch_ptr)[t3].emplace_back(index, this_slot);
          int t4 = t3 + g;
          if (t4 < (death_day-1) && t4 < max_time) {
            (*schedule_infective_chronic_ptr)[t4].emplace_back(index, this_slot);
          }
          
          int duration_chronic = get_duration_chronic();
          int t5 = t3 + duration_chronic;
          if (t5 < (death_day-1) && t5 < max_time) {
            
            (*schedule_Ch_to_Sh_ptr)[t5].emplace_back(index, this_slot);
            int t6 = t5 + g;
            if (t6 < death_day && t6 < max_time) {
              (*schedule_infective_recovery_ptr)[t6].emplace_back(index, this_slot);
            }
            
          }
          
        } else {
          
          (*schedule_Ah_to_Sh_ptr)[t3].emplace_back(index, this_slot);
          int t4 = t3 + g;
          if (t4 < (death_day-1) && t4 < max_time) {
            (*schedule_infective_recovery_ptr)[t4].emplace_back(index, this_slot);
          }
          
        }
      }
      
    } else {
      
      (*schedule_Eh_to_Ch_ptr)[t1].emplace_back(index, this_slot);
      int t2 = t1 + g;
      if (t2 < (death_day-1) && t2 < max_time) {
        (*schedule_infective_chronic_ptr)[t2].emplace_back(index, this_slot);
      }
      
      int duration_chronic = get_duration_chronic();
      int t3 = t1 + duration_chronic;
      if (t3 < (death_day-1) && t3 < max_time) {
        
        (*schedule_Ch_to_Sh_ptr)[t3].emplace_back(index, this_slot);
        int t4 = t3 + g;
        if (t4 < (death_day-1) && t4 < max_time) {
          (*schedule_infective_recovery_ptr)[t4].emplace_back(index, this_slot);
        }
        
      }
      
    }
  }
  
}

//------------------------------------------------
// move inoculation from Eh state to Ah
void Host::Eh_to_Ah(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage_acute() == 0) {
    (*Ah_ptr)[deme]++;
    if (get_n_bloodstage_chronic() > 0) {
      (*Ch_ptr)[deme]--;
    } else {
      (*Eh_ptr)[deme]--;
    }
  }
  
  // update status
  inoc_status_asexual[this_slot] = Acute_asexual;
  
}

//------------------------------------------------
// move inoculation from Eh state to Ch
void Host::Eh_to_Ch(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage_acute() == 0) {
    if (get_n_bloodstage_chronic() == 0) {
      (*Eh_ptr)[deme]--;
      (*Ch_ptr)[deme]++;
    }
  }
  
  // update status
  inoc_status_asexual[this_slot] = Chronic_asexual;
  
}

//------------------------------------------------
// move inoculation from Ah state to Ch
void Host::Ah_to_Ch(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage_acute() == 1) {
    (*Ah_ptr)[deme]--;
    (*Ch_ptr)[deme]++;
  }
  
  // update status
  inoc_status_asexual[this_slot] = Chronic_asexual;
  
}

//------------------------------------------------
// move inoculation from Ah state to Sh
void Host::Ah_to_Sh(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage_acute() == 1) {
    (*Ah_ptr)[deme]--;
    if (get_n_bloodstage_chronic() > 0) {
      (*Ch_ptr)[deme]++;
    } else if (get_n_liverstage() > 0) {
      (*Eh_ptr)[deme]++;
    } else {
      (*Sh_ptr)[deme]++;
    }
  }
  
  // update status
  inoc_status_asexual[this_slot] = Inactive_asexual;
  
}

//------------------------------------------------
// move inoculation from Ch state to Sh
void Host::Ch_to_Sh(int this_slot) {
  
  // update deme counts
  if (get_n_bloodstage_acute() == 0) {
    if (get_n_bloodstage_chronic() == 1) {
      (*Ch_ptr)[deme]--;
      if (get_n_liverstage() > 0) {
        (*Eh_ptr)[deme]++;
      } else {
        (*Sh_ptr)[deme]++;
      }
    }
  }
  
  // update status
  inoc_status_asexual[this_slot] = Inactive_asexual;
  
}

//------------------------------------------------
// begin acutely infective period
void Host::begin_infective_acute(int this_slot, int t) {
  
  // update host status
  inoc_status_sexual[this_slot] = Acute_sexual;
  
  // store time
  inoc_time_infective[this_slot] = t;
  
  // if newly infective then add to infectives list
  if (get_n_infective() == 1) {
    (*host_infective_index_ptr)[deme].push_back(index);
  }
  
}

//------------------------------------------------
// begin chronically infective period
void Host::begin_infective_chronic(int this_slot, int t) {
  
  // get current host status
  Status_sexual current_status = inoc_status_sexual[this_slot];
  
  // update host status
  inoc_status_sexual[this_slot] = Chronic_sexual;
  
  // store time
  inoc_time_infective[this_slot] = t;
  
  // if newly infective then add to infectives list
  if (current_status == Inactive_sexual && get_n_infective() == 1) {
    (*host_infective_index_ptr)[deme].push_back(index);
  }
  
}

//------------------------------------------------
// end infective period
void Host::end_infective(int this_slot) {
  
  // update host status
  inoc_status_sexual[this_slot] = Inactive_sexual;
  inoc_active[this_slot] = false;
  
  // if no longer infective then drop from infectives list
  if (get_n_infective() == 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
}

//------------------------------------------------
// get total number of active inoculations
int Host::get_n_active_inoc() {
  return sum_bool(inoc_active);
}

//------------------------------------------------
// get host status
Status_host Host::get_host_status() {
  Status_host ret = Host_Sh;
  for (int i = 0; i < max_inoculations; ++i) {
    if (!inoc_active[i]) {
      continue;
    }
    if (ret == Host_Sh && inoc_status_asexual[i] == Liverstage_asexual) {
      ret = Host_Eh;
    }
    if ((ret == Host_Sh || ret == Host_Eh) && inoc_status_asexual[i] == Chronic_asexual) {
      ret = Host_Ch;
    }
    if ((ret == Host_Sh || ret == Host_Eh || ret == Host_Ch) && inoc_status_asexual[i] == Acute_asexual) {
      ret = Host_Ah;
      break;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of liverstage inoculations
int Host::get_n_liverstage() {
  int ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_asexual[i] == Liverstage_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of acute bloodstage inoculations
int Host::get_n_bloodstage_acute() {
  int ret = 0;
  for (int i = 0; i < int(inoc_status_asexual.size()); ++i) {
    if (inoc_status_asexual[i] == Acute_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of chronic bloodstage inoculations
int Host::get_n_bloodstage_chronic() {
  int ret = 0;
  for (int i = 0; i < int(inoc_status_asexual.size()); ++i) {
    if (inoc_status_asexual[i] == Chronic_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of bloodstage inoculations
int Host::get_n_bloodstage() {
  return get_n_bloodstage_acute() + get_n_bloodstage_chronic();
}

//------------------------------------------------
// get total number of asexual stage inoculations
int Host::get_n_asexual() {
  return get_n_liverstage() + get_n_bloodstage();
}

//------------------------------------------------
// get total number of infective (sexual stage) inoculations
int Host::get_n_infective() {
  int ret = 0;
  for (int i = 0; i < int(inoc_status_sexual.size()); ++i) {
    if (inoc_status_sexual[i] == Acute_sexual || inoc_status_sexual[i] == Chronic_sexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  if (infection_index < (n_prob_infection-1)) {
    return prob_infection[infection_index];
  }
  return prob_infection[n_prob_infection-1];
}

//------------------------------------------------
// get current probability of going to acute infection
double Host::get_prob_acute() {
  if (inoc_index < (n_prob_acute-1)) {
    return prob_acute[inoc_index];
  }
  return prob_acute[n_prob_acute-1];
}

//------------------------------------------------
// get current probability of going to chronic from acute infection
double Host::get_prob_AC() {
  if (inoc_index < (n_prob_AC-1)) {
    return prob_AC[inoc_index];
  }
  return prob_AC[n_prob_AC-1];
}

//------------------------------------------------
// get duration of acute disease
int Host::get_duration_acute() {
  int ret;
  if (inoc_index < (n_duration_acute-1)) {
    ret = (*sampler_duration_acute_ptr)[inoc_index].draw() + 1;
  }
  ret = (*sampler_duration_acute_ptr)[n_duration_acute-1].draw() + 1;
  return ret;
}

//------------------------------------------------
// get duration of chronic disease
int Host::get_duration_chronic() {
  int ret;
  if (inoc_index < (n_duration_chronic-1)) {
    ret = (*sampler_duration_chronic_ptr)[inoc_index].draw() + 1;
  }
  ret = (*sampler_duration_chronic_ptr)[n_duration_chronic-1].draw() + 1;
  return ret;
}

//------------------------------------------------
// get current infectivity
double Host::get_infectivity(int t) {
  
  // TODO - infectivity from distribution
  /*
  // get probability of being infected by sexual-stage inoculation
  int ret = 0;
  for (int i = 0; i < int(inoc_status_sexual.size()); ++i) {
    if (inoc_status_sexual[i] == Active_sexual) {
      int time_diff = t - inoc_time_infective[i];
      //infectivity_acute
    }
  }
  */
  
  return max_infectivity;
}

//------------------------------------------------
// get next free inoculation slot
int Host::get_free_inoc_slot() {
  
  // loop until find inactive slot
  int ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (!inoc_active[i]) {
      break;
    }
    ret++;
  }
  
  // error check for outside range
  if (ret == max_inoculations) {
    Rcpp::stop("could not find free inoculation slot");
  }
  
  return ret;
}
