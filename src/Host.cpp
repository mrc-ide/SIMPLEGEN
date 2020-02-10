
#include "Host.h"
#include "misc_v9.h"
#include "probability_v10.h"

using namespace std;


// ################################################################################################
// STARTUP

//------------------------------------------------
// initialise host
void Host::init(int index, int &host_ID, int deme,
                vector<vector<int>> &host_index,
                vector<vector<int>> &host_infective_index,
                Sampler &sampler_age_stable,
                Sampler &sampler_age_death,
                vector<Sampler> &sampler_duration_acute,
                vector<Sampler> &sampler_duration_chronic,
                vector<Sampler> &sampler_time_treatment_acute,
                vector<Sampler> &sampler_time_treatment_chronic,
                Sampler &sampler_duration_prophylatic) {
  
  // identifiers
  this->index = index;
  this->host_ID = host_ID++;
  home_deme = deme;
  this->deme = deme;
  
  // pointer to indices of hosts and infective hosts in each deme
  host_index_ptr = &host_index;
  host_infective_index_ptr = &host_infective_index;
  
  // pointers to sampler objects, for efficiently drawing from global
  // distributions
  sampler_age_stable_ptr = &sampler_age_stable;
  sampler_age_death_ptr = &sampler_age_death;
  sampler_duration_acute_ptr = &sampler_duration_acute;
  sampler_duration_chronic_ptr = &sampler_duration_chronic;
  sampler_time_treatment_acute_ptr = &sampler_time_treatment_acute;
  sampler_time_treatment_chronic_ptr = &sampler_time_treatment_chronic;
  sampler_duration_prophylatic_ptr = &sampler_duration_prophylatic;
  
  // cumulative count of how many times this host has been bitten by infective
  // mosquito (infection_index) and how many times an infection has taken hold
  // (inoc_index)
  infection_index = 0;
  inoc_index = 0;
  
  // prophylactic status
  prophylaxis_on = false;
  t_prophylaxis_stop = 0;
  
  // initialise host-specific treatment-seeking probability
  draw_treatment_seeking();
  
  // initialise inoculation slots
  inoc_ID_vec = vector<int>(max_inoculations);
  inoc_active = vector<bool>(max_inoculations, false);
  inoc_status_asexual = vector<Status_asexual>(max_inoculations, Inactive_asexual);
  inoc_status_sexual = vector<Status_sexual>(max_inoculations, Inactive_sexual);
  inoc_time_asexual = vector<int>(max_inoculations);
  inoc_time_sexual = vector<int>(max_inoculations);
  
  // initialise events
  inoc_events = vector<map<Event, int>>(max_inoculations);
  t_next_inoc_event = max_time + 1;
  
  // draw age from stable demography distribution
  draw_starting_age();
  
}

//------------------------------------------------
// draw birth and death days given a known stable demography distribution, and
// schedule death for future. Only needed at start of simulation, after which
// hosts that die are instantly re-born and have death drawn from entire life
// table
void Host::draw_starting_age() {
  
  // draw age from stable demography distribution
  int age_years = sampler_age_stable_ptr->draw();
  int extra_days = sample2(0, 364);
  int age_days = age_years*365 + extra_days;
  
  // draw duration of life from demography distribution looking forward from
  // current age. This is tricky, as we must account for the fact that if we are
  // already part way into an age group then we have a reduced probability of
  // dying within that age group.
  //
  // math derivation:
  // assume constant rate of death r over 1-year period. Total probability of
  // dying this year is therefore p = 1 - exp(-r). We know p from life table,
  // from which we can derive r = -log(1-p). If we are already a proportion x
  // through this year, then the probability of dying in the remaining time is
  // Pr(die) = 1 - exp(-r(1-x)). Sustituting in r and simplifying we get Pr(die)
  // = 1 - (1-p)^(1-x). Notice that if p = 1 from the life table then we are
  // certain to die this year, irrespective of how far through the year we have
  // got already.
  int life_days = 0;
  double prop_year_remaining = 1.0 - extra_days/365.0; // (x in the above derivation)
  double prob_die_this_year = 1.0 - pow(1.0 - life_table[age_years], 1.0 - prop_year_remaining);
  if (rbernoulli1(prob_die_this_year)) {
    life_days = age_years*365 + sample2(extra_days, 364);
  } else {
    // if we do not die in year age_years then loop through all remaining years,
    // performing Bernoulli draw from probability of dying in that year
    for (unsigned int i = (age_years+1); i < life_table.size(); ++i) {
      if (rbernoulli1(life_table[i])) {
        life_days = i*365 + sample2(0, 364);
        break;
      }
    }
  }
  
  // calculate final birth and death days from age at time 0 and life_days
  birth_day = -age_days;
  death_day = life_days - age_days;
  
}

//------------------------------------------------
// draw host-specific treatment seeking parameter, either from Beta distribution
// or Dirac delta distribution
void Host::draw_treatment_seeking() {
  
  if (treatment_seeking_mean == 0 || treatment_seeking_mean == 1 || treatment_seeking_sd == 0) {
    treatment_seeking = treatment_seeking_mean;
  } else {
    double m = treatment_seeking_mean;
    double s = treatment_seeking_sd;
    double alpha = m*m*(1 - m)/(s*s) - m;
    double beta = m*(1 - m)*(1 - m)/(s*s) - (1 - m);
    treatment_seeking = rbeta1(alpha, beta);
  }
}

// ################################################################################################
// EVENT SCHEDULERS
// 
// rather than testing each generation to see if a host undergoes an event (e.g.
// acute infection transitioning to chronic), instead we draw the timings of all
// future events at the point of infection. These events are appended to
// inoc_events for a particular inoculation slot. The value t_next_inoc_event is
// used to describe the time of the very next event over all inoculation slots,
// meaning all we have to do in a given generation is test whether the host has
// an event scheduled for this day, and if not we move on to the next host.
// 
// When a scheduled event occurs we carry out the appropriate function, and
// simultaneously recalculate t_next_inoc_event to be the time of the next
// event. Similarly, when a new event is scheduled we append it to inoc_events,
// and simultaneously ask whether the new event occurs before the current
// t_next_inoc_event, in which case it replaces it.

//------------------------------------------------
// add new inoc event to list
void Host::new_inoc_event(int t, Event this_event, int this_slot) {
  t_next_inoc_event = (t < t_next_inoc_event) ? t : t_next_inoc_event;
  inoc_events[this_slot].insert(make_pair(this_event, t));
}

//------------------------------------------------
// check for any inoc-level events that need to be carried out at time t
void Host::check_inoc_event(int t) {
  
  // return if nothing scheduled at time t
  if (t_next_inoc_event != t) {
    return;
  }
  
  // deal with treatment at time t first
  bool treatment_break = false;
  for (int i = 0; i < max_inoculations; ++i) {
    for (const auto & x : inoc_events[i]) {
      if (x.first == Event_treatment && x.second == t) {
        treatment(t);
        treatment_break = true;
        break;
      }
    }
    if (treatment_break) {
      break;
    }
  }
  
  // find other events that are scheduled for time t
  t_next_inoc_event = max_time + 1;
  for (int i = 0; i < max_inoculations; ++i) {
    
    // check for event at this slot that is scheduled for time t
    for (auto it = inoc_events[i].begin(); it != inoc_events[i].end();) {
      if (it->second == t) {
        
        // switch depending on event
        switch(it->first) {
        case Event_Eh_to_Ah:
          Eh_to_Ah(i, t);
          break;
        case Event_Eh_to_Ch:
          Eh_to_Ch(i, t);
          break;
        case Event_Ah_to_Ch:
          Ah_to_Ch(i, t);
          break;
        case Event_Ah_to_Sh:
          Ah_to_Sh(i, t);
          break;
        case Event_Ch_to_Sh:
          Ch_to_Sh(i, t);
          break;
        case Event_begin_infective_acute:
          begin_infective_acute(i, t);
          break;
        case Event_begin_infective_chronic:
          begin_infective_chronic(i, t);
          break;
        case Event_end_infective:
          end_infective(i);
          break;
        default:
          break;
        }
        
        // drop this event from inoc_events[i]
        it = inoc_events[i].erase(it);
      } else {
        if (it->second < t_next_inoc_event) {
          t_next_inoc_event = it->second;
        }
        ++it;
      }
    }
    
  }  // end i loop
  
}


// ################################################################################################
// HOST-LEVEL EVENTS
// 
// Functions for carrying out the main work of different types of event. Some of
// these, for example infection, apply to a single inoculation slot and result
// in new events being scheduled. Others, such as treatment and host death,
// affect multiple inoculation slots, along with other host properties.
// 
// Some of these functions can be fiddly. For example, on treatment we move all
// inoc_status_asexual slots to Inactive_asexual, which is simple enough, but we
// also need to go through all future scheduled events for the progression of
// this inoculation, and update or remove them as needed.

//------------------------------------------------
// death
void Host::death(int &host_ID, int t) {
  
  // move host back to home deme, dealing with infective list etc.
  migrate(home_deme);
  
  // set current deme equal to home deme
  //deme = home_deme;
  
  // new unique ID
  this->host_ID = host_ID++;
  
  // reset cumulative counts
  infection_index = 0;
  inoc_index = 0;
  
  // reset prophylactic status
  prophylaxis_on = false;
  t_prophylaxis_stop = 0;
  
  // new host will be re-born today
  birth_day = t;
  
  // draw life duration from demography distribution
  int life_years = sampler_age_death_ptr->draw();
  int life_days = life_years*365 + sample2(0, 364);
  
  // cannot die on same day born
  if (life_days == 0) {
    life_days = 1;
  }
  death_day = birth_day + life_days;
  
  // re-draw host-specific treatment-seeking probability
  draw_treatment_seeking();
  
  // reset inoculation slots
  fill(inoc_ID_vec.begin(), inoc_ID_vec.end(), 0);
  fill(inoc_active.begin(), inoc_active.end(), false);
  fill(inoc_status_asexual.begin(), inoc_status_asexual.end(), Inactive_asexual);
  fill(inoc_status_sexual.begin(), inoc_status_sexual.end(), Inactive_sexual);
  fill(inoc_time_asexual.begin(), inoc_time_asexual.end(), 0);
  fill(inoc_time_sexual.begin(), inoc_time_sexual.end(), 0);
  
  // clear all scheduled inoc events
  for (int i = 0; i < max_inoculations; ++i) {
    inoc_events[i].clear();
  }
  t_next_inoc_event = max_time + 1;
  
}

//------------------------------------------------
// de-novo infection
void Host::denovo_infection(int t, int &next_inoc_ID, std::ofstream &transmission_record) {
  
  // generate dummy mosquito with an empty vector of inoculations
  vector<int> empty_vec;
  Mosquito dummy_mosquito(empty_vec);
  
  // carry out infection from dummy mosquito
  infection(t, next_inoc_ID, dummy_mosquito, transmission_record);
  
}

//------------------------------------------------
// infection
void Host::infection(int t, int &next_inoc_ID, Mosquito &mosq, std::ofstream &transmission_record) {
  
  //print(get_n_active_inoc());
  // return if already at max_inoculations
  if (get_n_active_inoc() == max_inoculations) {
    infection_index++;
    return;
  }
  
  // return if in prophlactic state
  if (prophylaxis_on) {
    infection_index++;
    return;
  }
  
  // get next free inoculation slot
  int this_slot = get_free_inoc_slot();
  
  // add new inoculation
  inoc_ID_vec[this_slot] = next_inoc_ID;
  inoc_active[this_slot] = true;
  inoc_status_asexual[this_slot] = Liverstage_asexual;
  
  // schedule disease progression. This is where we look through the tree of all
  // possible future trajectories, for example whether the inoculation
  // transitions to acute stage or directly to chronic stage etc, and schedule
  // these events to happen
  int time_consider_treatment = 0;
  bool seek_treatment = false;
  int t1 = t + u;
  bool acute = rbernoulli1(get_prob_acute());
  if (acute) {
    
    // schedule change of state
    new_inoc_event(t1, Event_Eh_to_Ah, this_slot);
    
    // schedule become acutely infective
    int t2 = t1 + g;
    new_inoc_event(t2, Event_begin_infective_acute, this_slot);
    
    // draw duration of acute phase
    int dur_acute = get_duration_acute();
    int t3 = t1 + dur_acute;
    
    // draw time to considering seeking treatment
    time_consider_treatment = get_time_treatment_acute();
    
    // determine whether host actively seeks treatment
    if (time_consider_treatment <= dur_acute) {
      seek_treatment = rbernoulli1(treatment_seeking);
    }
    
    // whether to transition to chronic stage prior to recovery
    bool acute_to_chronic = rbernoulli1(get_prob_AC());
    if (acute_to_chronic) {
      
      // schedule change of state
      new_inoc_event(t3, Event_Ah_to_Ch, this_slot);
      
      // schedule become chronically infective
      int t4 = t3 + g;
      new_inoc_event(t4, Event_begin_infective_chronic, this_slot);
      
      // draw duration of chronic phase
      int dur_chronic = get_duration_chronic();
      int t5 = t3 + dur_chronic;
      
      // schedule change of state
      new_inoc_event(t5, Event_Ch_to_Sh, this_slot);
      
      // schedule recover from infective
      int t6 = t5 + g;
      new_inoc_event(t6, Event_end_infective, this_slot);
      
    } else {  // transition directly from acute to recovery
      
      // schedule change of state
      new_inoc_event(t3, Event_Ah_to_Sh, this_slot);
      
      // schedule recover from infective
      int t4 = t3 + g;
      new_inoc_event(t4, Event_end_infective, this_slot);
      
    }
    
  } else {  // transition directly to chronic stage
    
    // schedule change of state
    new_inoc_event(t1, Event_Eh_to_Ch, this_slot);
    
    // schedule become chronically infective
    int t2 = t1 + g;
    new_inoc_event(t2, Event_begin_infective_chronic, this_slot);
    
    // draw duration of chronic phase
    int dur_chronic = get_duration_chronic();
    int t3 = t1 + dur_chronic;
    
    // draw time to considering seeking treatment
    time_consider_treatment = get_time_treatment_chronic();
    
    // determine whether host actively seeks treatment
    if (time_consider_treatment <= dur_chronic) {
      seek_treatment = rbernoulli1(treatment_seeking);
    }
    
    // schedule change of state
    new_inoc_event(t3, Event_Ch_to_Sh, this_slot);
    
    // schedule recover from infective
    int t4 = t3 + g;
    new_inoc_event(t4, Event_end_infective, this_slot);
    
  }
  
  // schedule treatment
  if (seek_treatment) {
    int t2 = t1 + time_consider_treatment;
    new_inoc_event(t2, Event_treatment, this_slot);
  }
  
  // update indices
  infection_index++;
  inoc_index++;
  //print(inoc_index);
  
  // print to transmission record
  if (save_transmission_record) {
    transmission_record << next_inoc_ID;
    for (const auto &x : mosq.inoc_ID) {
      transmission_record << " " << x;
    }
    transmission_record << ";";
  }
  
  // update inoc ID
  next_inoc_ID++;
  
}

//------------------------------------------------
// treatment
void Host::treatment(int t) {
  
  // draw duration of prophylaxis
  int dur_prophylaxis = get_duration_prophylaxis();
  int t2 = t + dur_prophylaxis;
  
  // activate prophylaxis
  if (dur_prophylaxis > 0) {
    prophylaxis_on = true;
    t_prophylaxis_stop = t2;
  }
  
  // loop through inoculations
  for (int i = 0; i < max_inoculations; ++i) {
    
    // anything that is currently in the acute or chronic stages is cured
    if (inoc_status_asexual[i] == Acute_asexual || inoc_status_asexual[i] == Chronic_asexual) {
      
      // if due to become infective at a future timepoint then store this
      // timepoint
      bool due_acute_infective = (inoc_events[i].count(Event_begin_infective_acute) != 0);
      bool due_chronic_infective = (inoc_events[i].count(Event_begin_infective_chronic) != 0);
      int t_infective_start;
      if (due_acute_infective) {
        t_infective_start = inoc_events[i][Event_begin_infective_acute];
      } else if (due_chronic_infective) {
        t_infective_start = inoc_events[i][Event_begin_infective_chronic];
      }
      
      // clear bloodstage infection
      inoc_status_asexual[i] = Inactive_asexual;
      
      // clear all scheduled events
      inoc_events[i].clear();
      
      // reinstate infectious start if needed
      if (due_acute_infective) {
        new_inoc_event(t_infective_start, Event_begin_infective_acute, i);
      } else if (due_chronic_infective) {
        new_inoc_event(t_infective_start, Event_begin_infective_chronic, i);
      }
      
      // schedule new infectious end
      new_inoc_event(t + g, Event_end_infective, i);
      
    }
    
    // anything that is scheduled to emerge from the liver during the
    // prophylactic period is cured immediately upon emergence
    if (inoc_status_asexual[i] == Liverstage_asexual) {
      
      if (inoc_events[i].count(Event_Eh_to_Ah) != 0 && inoc_events[i][Event_Eh_to_Ah] <= t2) {
        
        // store time at which due to emerge
        int t_emerge = inoc_events[i][Event_Eh_to_Ah];
        
        // clear all scheduled events
        inoc_events[i].clear();
        
        // reschedule progression directly to clearance
        new_inoc_event(t_emerge, Event_Eh_to_Ah, i);
        new_inoc_event(t_emerge, Event_Ah_to_Sh, i);
        new_inoc_event(t_emerge, Event_end_infective, i);
        
      } else if (inoc_events[i].count(Event_Eh_to_Ch) != 0 && inoc_events[i][Event_Eh_to_Ch] <= t2) {
        
        // store time at which due to emerge
        int t_emerge = inoc_events[i][Event_Eh_to_Ch];
        
        // clear all scheduled events
        inoc_events[i].clear();
        
        // reschedule progression directly to clearance
        new_inoc_event(t_emerge, Event_Eh_to_Ch, i);
        new_inoc_event(t_emerge, Event_Ch_to_Sh, i);
        new_inoc_event(t_emerge, Event_end_infective, i);
        
      }
    }
    
  }  // end i loop
  
  // recalculate time of next event
  t_next_inoc_event = max_time + 1;
  for (int i = 0; i < max_inoculations; ++i) {
    for (const auto & x : inoc_events[i]) {
      if (x.second < t) {
        print("time =", t);
        print_inoc_status();
        print_inoc_events();
        Rcpp::stop("in Host::treatment, events found that predate current time");
      }
      if (x.second < t_next_inoc_event) {
        t_next_inoc_event = x.second;
      }
    }
  }
  
}

//------------------------------------------------
// end prophylactic period
void Host::end_prophylaxis() {
  prophylaxis_on = false;
}

//------------------------------------------------
// migrate to new deme
void Host::migrate(int new_deme) {
  
  // return if no migration
  if (new_deme == deme) {
    return;
  }
  
  // move host index between demes
  erase_remove((*host_index_ptr)[deme], index);
  (*host_index_ptr)[new_deme].push_back(index);
  
  // if infective move host index between demes
  if (get_n_infective() != 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
    (*host_infective_index_ptr)[new_deme].push_back(index);
  }
  
  // update deme
  deme = new_deme;
  
}

// ################################################################################################
// INOCULATION-LEVEL EVENTS

//------------------------------------------------
// move inoculation from Eh state to Ah
void Host::Eh_to_Ah(int this_slot, int t) {
  
  // update status
  inoc_status_asexual[this_slot] = Acute_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
}

//------------------------------------------------
// move inoculation from Eh state to Ch
void Host::Eh_to_Ch(int this_slot, int t) {
  
  // update status
  inoc_status_asexual[this_slot] = Chronic_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
}

//------------------------------------------------
// move inoculation from Ah state to Ch
void Host::Ah_to_Ch(int this_slot, int t) {
  
  // update status
  inoc_status_asexual[this_slot] = Chronic_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
}

//------------------------------------------------
// move inoculation from Ah state to Sh
void Host::Ah_to_Sh(int this_slot, int t) {
  
  // update status
  inoc_status_asexual[this_slot] = Inactive_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
}

//------------------------------------------------
// move inoculation from Ch state to Sh
void Host::Ch_to_Sh(int this_slot, int t) {
  
  // update status
  inoc_status_asexual[this_slot] = Inactive_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
}

//------------------------------------------------
// begin acutely infective period
void Host::begin_infective_acute(int this_slot, int t) {
  
  // update host status
  inoc_status_sexual[this_slot] = Acute_sexual;
  
  // store time
  inoc_time_sexual[this_slot] = t;
  
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
  inoc_time_sexual[this_slot] = t;
  
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
  
  // reset times
  inoc_time_asexual[this_slot] = 0;
  inoc_time_sexual[this_slot] = 0;
  
  // if no longer infective then drop from infectives list
  if (get_n_infective() == 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
}


// ################################################################################################
// GETTERS AND SETTERS

//------------------------------------------------
// return vector of inoc IDs taken from infective inoculations
vector<int> Host::get_inoc_ID_vec() {
  vector<int> ret;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_sexual[i] != Inactive_sexual) {
      ret.push_back(inoc_ID_vec[i]);
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of active inoculations
int Host::get_n_active_inoc() {
  return sum_bool(inoc_active);
}

//------------------------------------------------
// get host status
Status_host Host::get_host_status() {
  if (prophylaxis_on) {
    return Host_Ph;
  }
  Status_host ret = Host_Sh;
  for (int i = 0; i < max_inoculations; ++i) {
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
  for (unsigned int i = 0; i < inoc_status_asexual.size(); ++i) {
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
  for (unsigned int i = 0; i < inoc_status_asexual.size(); ++i) {
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
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_sexual[i] == Acute_sexual || inoc_status_sexual[i] == Chronic_sexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  int index = (infection_index < n_prob_infection) ? infection_index : n_prob_infection - 1;
  return prob_infection[index];
}

//------------------------------------------------
// get current probability of going to acute infection
double Host::get_prob_acute() {
  int index = (inoc_index < n_prob_acute) ? inoc_index : n_prob_acute - 1;
  return prob_acute[index];
}

//------------------------------------------------
// get current probability of going to chronic from acute infection
double Host::get_prob_AC() {
  int index = (inoc_index < n_prob_AC) ? inoc_index : n_prob_AC - 1;
  return prob_AC[index];
}

//------------------------------------------------
// get duration of acute disease
int Host::get_duration_acute() {
  int index = (inoc_index < n_duration_acute) ? inoc_index : n_duration_acute - 1;
  return (*sampler_duration_acute_ptr)[index].draw();
}

//------------------------------------------------
// get duration of chronic disease
int Host::get_duration_chronic() {
  int index = (inoc_index < n_duration_chronic) ? inoc_index : n_duration_chronic - 1;
  return (*sampler_duration_chronic_ptr)[index].draw();
}

//------------------------------------------------
// get time until treatment of acute disease
int Host::get_time_treatment_acute() {
  int index = (inoc_index < n_time_treatment_acute) ? inoc_index : n_time_treatment_acute - 1;
  return (*sampler_time_treatment_acute_ptr)[index].draw();
}

//------------------------------------------------
// get time until treatment of chronic disease
int Host::get_time_treatment_chronic() {
  int index = (inoc_index < n_time_treatment_chronic) ? inoc_index : n_time_treatment_chronic - 1;
  return (*sampler_time_treatment_chronic_ptr)[index].draw();
}

//------------------------------------------------
// get duration of prophylaxis
int Host::get_duration_prophylaxis() {
  return (*sampler_duration_prophylatic_ptr).draw();
}

//------------------------------------------------
// get current detectability by microscopy over acute inoculations
double Host::get_detectability_microscopy_acute(int t) {
  
  // acute detectability is the max detectability over all acute inoculations
  double ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_asexual[i] == Acute_asexual) {
      
      // get time since became infective
      int time_diff = t - inoc_time_asexual[i];
      
      // get infectivity from appropriate distribution
      int tmp1 = (inoc_index < n_detectability_microscopy_acute) ? inoc_index : n_detectability_microscopy_acute - 1;
      int tmp2 = (time_diff < int(detectability_microscopy_acute[tmp1].size())) ? time_diff : int(detectability_microscopy_acute[tmp1].size()) - 1;
      double inoc_detectability = detectability_microscopy_acute[tmp1][tmp2];
      
      // ret is max of ret and inoc_detectability
      ret = (ret > inoc_detectability) ? ret : inoc_detectability;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get current detectability by microscopy over chronic inoculations
double Host::get_detectability_microscopy_chronic(int t) {
  
  // acute detectability is the max detectability over all chronic inoculations
  double ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_asexual[i] == Chronic_asexual) {
      
      // get time since became infective
      int time_diff = t - inoc_time_asexual[i];
      
      // get infectivity from appropriate distribution
      int tmp1 = (inoc_index < n_detectability_microscopy_chronic) ? inoc_index : n_detectability_microscopy_chronic - 1;
      int tmp2 = (time_diff < int(detectability_microscopy_chronic[tmp1].size())) ? time_diff : int(detectability_microscopy_chronic[tmp1].size()) - 1;
      double inoc_detectability = detectability_microscopy_chronic[tmp1][tmp2];
      
      // ret is max of ret and inoc_detectability
      ret = (ret > inoc_detectability) ? ret : inoc_detectability;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get current detectability by PCR over acute inoculations
double Host::get_detectability_PCR_acute(int t) {
  
  // acute detectability is the max detectability over all acute inoculations
  double ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_asexual[i] == Acute_asexual) {
      
      // get time since became infective
      int time_diff = t - inoc_time_asexual[i];
      
      // get infectivity from appropriate distribution
      int tmp1 = (inoc_index < n_detectability_PCR_acute) ? inoc_index : n_detectability_PCR_acute - 1;
      int tmp2 = (time_diff < int(detectability_PCR_acute[tmp1].size())) ? time_diff : int(detectability_PCR_acute[tmp1].size()) - 1;
      double inoc_detectability = detectability_PCR_acute[tmp1][tmp2];
      
      // ret is max of ret and inoc_detectability
      ret = (ret > inoc_detectability) ? ret : inoc_detectability;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get current detectability by PCR over chronic inoculations
double Host::get_detectability_PCR_chronic(int t) {
  
  // acute detectability is the max detectability over all chronic inoculations
  double ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_asexual[i] == Chronic_asexual) {
      
      // get time since became infective
      int time_diff = t - inoc_time_asexual[i];
      
      // get infectivity from appropriate distribution
      int tmp1 = (inoc_index < n_detectability_PCR_chronic) ? inoc_index : n_detectability_PCR_chronic - 1;
      int tmp2 = (time_diff < int(detectability_PCR_chronic[tmp1].size())) ? time_diff : int(detectability_PCR_chronic[tmp1].size()) - 1;
      double inoc_detectability = detectability_PCR_chronic[tmp1][tmp2];
      
      // ret is max of ret and inoc_detectability
      ret = (ret > inoc_detectability) ? ret : inoc_detectability;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get current infectivity
double Host::get_infectivity(int t) {
  
  // host infectivity is the max infectivity over all inoculations
  double ret = 0;
  for (int i = 0; i < max_inoculations; ++i) {
    if (inoc_status_sexual[i] == Acute_sexual || inoc_status_sexual[i] == Chronic_sexual) {
      
      // get time since became infective
      int time_diff = t - inoc_time_sexual[i];
      
      // split by acute vs. chronic infectivity
      double inoc_infectivity;
      if (inoc_status_sexual[i] == Acute_sexual) {
        
        // get infectivity from appropriate distribution
        int tmp1 = (inoc_index < n_infectivity_acute) ? inoc_index : n_infectivity_acute-1;
        int tmp2 = (time_diff < int(infectivity_acute[tmp1].size())) ? time_diff : int(infectivity_acute[tmp1].size())-1;
        inoc_infectivity = infectivity_acute[tmp1][tmp2];
        
      } else {
        
        // get infectivity from appropriate distribution
        int tmp1 = (inoc_index < n_infectivity_chronic) ? inoc_index : n_infectivity_chronic-1;
        int tmp2 = (time_diff < int(infectivity_chronic[tmp1].size())) ? time_diff : int(infectivity_chronic[tmp1].size())-1;
        inoc_infectivity = infectivity_chronic[tmp1][tmp2];
        
      }
      
      // ret is max of ret and inoc_infectivity
      ret = (ret > inoc_infectivity) ? ret : inoc_infectivity;
    }
  }
  
  return ret;
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

//------------------------------------------------
// get host age (years) at time t
int Host::get_age(int t) {
  return (t - birth_day)/365;
}

//------------------------------------------------
// print innoc_events
void Host::print_inoc_events() {
  
  for (int i = 0; i < max_inoculations; ++i) {
    for (const auto & x : inoc_events[i]) {
      Rcpp::Rcout << "[" << x.first << ", " << x.second << "] ";
    }
    Rcpp::Rcout << "\n";
  }
}

//------------------------------------------------
// print status of inoc slots
void Host::print_inoc_status() {
  Rcpp::Rcout << "host_ID: " << host_ID << "\n";
  Rcpp::Rcout << "inoc_ID: ";
  print_vector(inoc_ID_vec);
  Rcpp::Rcout << "active : ";
  print_vector(inoc_active);
  Rcpp::Rcout << "asex   : ";
  print_vector(inoc_status_asexual);
  Rcpp::Rcout << "sex    : ";
  print_vector(inoc_status_sexual);
  Rcpp::Rcout << "time   : ";
  print_vector(inoc_time_sexual);
}
