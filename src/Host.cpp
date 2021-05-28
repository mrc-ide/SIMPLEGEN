
#include "Host.h"
#include "misc_v10.h"
#include "probability_v11.h"

using namespace std;


// #############################################################################
// #                                                                           #
// #   INITIALISATION                                                          #
// #                                                                           #
// #############################################################################

//------------------------------------------------
// initialise host
void Host::init(Parameters &params_,
                int index, int &host_ID, int deme,
                vector<vector<int>> &host_index,
                vector<vector<int>> &host_infective_index,
                Sampler &sampler_age_stable,
                Sampler &sampler_age_death,
                vector<Sampler> &sampler_duration_acute,
                vector<Sampler> &sampler_duration_chronic,
                vector<Sampler> &sampler_time_treatment_acute,
                vector<Sampler> &sampler_time_treatment_chronic,
                vector<Sampler> &sampler_duration_prophylatic) {
  
  // store pointer to parameters
  params = &params_;
  
  // identifiers. See host.h for further description
  this->index = index;
  this->host_ID = host_ID++;
  home_deme = deme;
  this->deme = deme;
  
  // pointers to indices of hosts and infective hosts in each deme
  host_index_ptr = &host_index;
  host_infective_index_ptr = &host_infective_index;
  
  // pointers to sampler objects. Used to make efficient random draws from
  // global distributions
  sampler_age_stable_ptr = &sampler_age_stable;
  sampler_age_death_ptr = &sampler_age_death;
  sampler_duration_acute_ptr = &sampler_duration_acute;
  sampler_duration_chronic_ptr = &sampler_duration_chronic;
  sampler_time_treatment_acute_ptr = &sampler_time_treatment_acute;
  sampler_time_treatment_chronic_ptr = &sampler_time_treatment_chronic;
  sampler_duration_prophylatic_ptr = &sampler_duration_prophylatic;
  
  // cumulative count of how many times this host has been bitten by infective
  // mosquito (cumul_infections) and how many times an infection has taken hold
  // (cumul_inoculations)
  cumul_infections = 0;
  cumul_inoculations = 0;
  
  // draw birth and death days from stable demography distribution
  draw_starting_age();
  
  // prophylactic status
  prophylaxis_on = false;
  t_prophylaxis_stop = 0;
  
  // host state now, previously, and time at which it changed
  host_state = Host_Sh;
  host_state_previous = Host_Sh;
  time_host_state_change = 0;
  
  // draw host-specific treatment-seeking probability
  draw_treatment_seeking();
  
  // initialise inoculation slots
  inoc_ID_vec = vector<int>(params->max_inoculations);
  inoc_active = vector<bool>(params->max_inoculations, false);
  inoc_state_asexual = vector<State_asexual>(params->max_inoculations, Inactive_asexual);
  inoc_state_sexual = vector<State_sexual>(params->max_inoculations, Inactive_sexual);
  inoc_time_asexual = vector<int>(params->max_inoculations);
  inoc_time_sexual = vector<int>(params->max_inoculations);
  
  // initialise events
  inoc_events = vector<map<Event, int>>(params->max_inoculations);
  t_next_inoc_event = params->max_time + 1;
  
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
  // = 1 - (1-p)^(1-x). Notice that as x tends to 1 the chance of dying this
  // year tends to 0, UNLESS p = 1 in which case death is certain this year
  // irrespective of x.
  int life_days = 0;
  double prop_year_remaining = 1.0 - extra_days/365.0; // (x in the above derivation)
  double prob_die_this_year = 1.0 - pow(1.0 - params->life_table[age_years], 1.0 - prop_year_remaining);
  if (rbernoulli1(prob_die_this_year)) {
    life_days = age_years*365 + sample2(extra_days, 364);
  } else {
    // if we do not die in the current year of age then loop through all
    // remaining years, performing Bernoulli draw from probability of dying in
    // that year
    for (unsigned int i = (age_years+1); i < params->life_table.size(); ++i) {
      if (rbernoulli1(params->life_table[i])) {
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
  
  if ((params->treatment_seeking_mean == 0) || (params->treatment_seeking_mean == 1) || (params->treatment_seeking_sd == 0)) {
    treatment_seeking = params->treatment_seeking_mean;
  } else {
    double m = params->treatment_seeking_mean;
    double s = params->treatment_seeking_sd;
    double alpha = m*m*(1 - m)/(s*s) - m;
    double beta = m*(1 - m)*(1 - m)/(s*s) - (1 - m);
    treatment_seeking = rbeta1(alpha, beta);
  }
}

// #############################################################################
// #                                                                           #
// #   EVENT SCHEDULERS                                                        #
// #                                                                           #
// #############################################################################

// rather than testing each generation to see if a host undergoes an event (e.g.
// acute infection transitioning to chronic), instead we draw the timings of all
// future events at the point of infection. These events are appended to the
// inoc_events vector. The value t_next_inoc_event is used to describe the time
// of the next event in the life of the host, meaning all we have to do in a
// given generation is test whether the host has an event scheduled for day
// t_next_inoc_event, and if not we move on to the next host.
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
  for (int i = 0; i < params->max_inoculations; ++i) {
    for (const auto & x : inoc_events[i]) {
      if ((x.first == Event_treatment) && (x.second == t)) {
        treatment(t);
        goto treatment_break;
      }
    }
  }
  treatment_break:
  
  // find other events that are scheduled for time t. Simultaneously recalculate
  // the time of next inoc event
  t_next_inoc_event = params->max_time + 1;
  for (int i = 0; i < params->max_inoculations; ++i) {
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


// #############################################################################
// #                                                                           #
// #   HOST-LEVEL EVENTS                                                       #
// #                                                                           #
// #############################################################################

// Functions for carrying out the main work of different types of event. Some of
// these, for example infection, apply to a single inoculation slot and result
// in new events being scheduled. Others, such as treatment and host death,
// affect multiple inoculation slots along with other host properties.
// 
// Some of these functions can be fiddly. For example, on treatment we move all
// inoc_state_asexual slots to Inactive_asexual, which is simple enough, but we
// also need to go through all future scheduled events for the progression of
// this inoculation and update or remove them as needed. For example,
// progression to chronic infection may be curtailed, and the time of asexual
// stages terminating will also likely change.

//------------------------------------------------
// death
void Host::death(int &host_ID, int t) {
  
  // move host back to home deme
  migrate(home_deme);
  
  // drop from infectives list if present
  erase_remove((*host_infective_index_ptr)[deme], index);
  
  // new unique ID
  this->host_ID = host_ID++;
  
  // reset cumulative counts
  cumul_infections = 0;
  cumul_inoculations = 0;
  
  // reset prophylactic status
  prophylaxis_on = false;
  t_prophylaxis_stop = 0;
  
  // reset host state
  host_state = Host_Sh;
  host_state_previous = Host_Sh;
  time_host_state_change = t;
  
  // reset inoculation slots
  fill(inoc_ID_vec.begin(), inoc_ID_vec.end(), 0);
  fill(inoc_active.begin(), inoc_active.end(), false);
  fill(inoc_state_asexual.begin(), inoc_state_asexual.end(), Inactive_asexual);
  fill(inoc_state_sexual.begin(), inoc_state_sexual.end(), Inactive_sexual);
  fill(inoc_time_asexual.begin(), inoc_time_asexual.end(), 0);
  fill(inoc_time_sexual.begin(), inoc_time_sexual.end(), 0);
  
  // clear all scheduled inoc events
  for (int i = 0; i < params->max_inoculations; ++i) {
    inoc_events[i].clear();
  }
  t_next_inoc_event = params->max_time + 1;
  
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
  
}

//------------------------------------------------
// migrate to new deme
void Host::migrate(int new_deme) {
  
  // return if migration is to current deme
  if (new_deme == deme) {
    return;
  }
  
  // exit if host index cannot be found in this deme
  if (!is_in_vector(index, (*host_index_ptr)[deme])) {
    Rcpp::stop("error in migrate() function: index not found in host_index vector");
  }
  
  // move host index between demes
  erase_remove((*host_index_ptr)[deme], index);
  (*host_index_ptr)[new_deme].push_back(index);
  
  // if infective then move host index between demes
  if (get_n_infective() != 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
    (*host_infective_index_ptr)[new_deme].push_back(index);
  }
  
  // update deme
  deme = new_deme;
  
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
  
  // update cumulative infections
  cumul_infections++;
  
  // return if already at max_inoculations
  if (get_n_active_inoc() == params->max_inoculations) {
    return;
  }
  
  // return if in prophlactic state
  if (prophylaxis_on) {
    return;
  }
  
  // update cumulative inoculations
  cumul_inoculations++;
  
  // print to transmission record
  if (params->save_transmission_record) {
    transmission_record << next_inoc_ID;
    for (const auto &x : mosq.inoc_ID) {
      transmission_record << " " << x;
    }
    transmission_record << ";";
  }
  
  // get next free inoculation slot
  int this_slot = get_free_inoc_slot();
  
  // add new inoculation and increment ID
  new_Eh(this_slot, t, next_inoc_ID);
  
  // schedule disease progression. This is where we look through the tree of all
  // possible future trajectories, for example whether the inoculation
  // transitions to acute stage or directly to chronic stage etc, and schedule
  // these events to happen
  int time_consider_treatment = 0;
  bool seek_treatment = false;
  int t1 = t + params->u;
  bool acute = rbernoulli1(get_prob_acute());
  if (acute) {
    
    // schedule change of state
    new_inoc_event(t1, Event_Eh_to_Ah, this_slot);
    
    // schedule become acutely infective
    int t2 = t1 + params->g;
    new_inoc_event(t2, Event_begin_infective_acute, this_slot);
    
    // draw duration of acute phase
    int dur_acute = draw_duration_acute();
    int t3 = t1 + dur_acute;
    
    // draw time to considering seeking treatment
    time_consider_treatment = draw_time_treatment_acute();
    
    // if consider seeking treatment before infection clears naturally,
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
      int t4 = t3 + params->g;
      new_inoc_event(t4, Event_begin_infective_chronic, this_slot);
      
      // draw duration of chronic phase
      int dur_chronic = draw_duration_chronic();
      int t5 = t3 + dur_chronic;
      
      // schedule change of state
      new_inoc_event(t5, Event_Ch_to_Sh, this_slot);
      
      // schedule recover from infective
      int t6 = t5 + params->g;
      new_inoc_event(t6, Event_end_infective, this_slot);
      
    } else {  // transition directly from acute to recovery
      
      // schedule change of state
      new_inoc_event(t3, Event_Ah_to_Sh, this_slot);
      
      // schedule recover from infective
      int t4 = t3 + params->g;
      new_inoc_event(t4, Event_end_infective, this_slot);
      
    }
    
  } else {  // transition directly to chronic stage
    
    // schedule change of state
    new_inoc_event(t1, Event_Eh_to_Ch, this_slot);
    
    // schedule become chronically infective
    int t2 = t1 + params->g;
    new_inoc_event(t2, Event_begin_infective_chronic, this_slot);
    
    // draw duration of chronic phase
    int dur_chronic = draw_duration_chronic();
    int t3 = t1 + dur_chronic;
    
    // draw time to considering seeking treatment
    time_consider_treatment = draw_time_treatment_chronic();
    
    // if consider seeking treatment before infection clears naturally,
    // determine whether host actively seeks treatment
    if (time_consider_treatment <= dur_chronic) {
      seek_treatment = rbernoulli1(treatment_seeking);
    }
    
    // schedule change of state
    new_inoc_event(t3, Event_Ch_to_Sh, this_slot);
    
    // schedule recover from infective
    int t4 = t3 + params->g;
    new_inoc_event(t4, Event_end_infective, this_slot);
    
  }
  
  // schedule treatment
  if (seek_treatment) {
    int t2 = t1 + time_consider_treatment;
    new_inoc_event(t2, Event_treatment, this_slot);
  }
  
}

//------------------------------------------------
// treatment
void Host::treatment(int t) {
  
  // draw duration of prophylaxis
  int dur_prophylaxis = draw_duration_prophylaxis();
  int t2 = t + dur_prophylaxis;
  
  // activate prophylaxis
  if (dur_prophylaxis > 0) {
    prophylaxis_on = true;
    t_prophylaxis_stop = t2;
  }
  
  // recalculate host state
  recalculate_host_state(t);
  
  // loop through inoculations
  for (int i = 0; i < params->max_inoculations; ++i) {
    
    // anything that is currently in the acute or chronic stages is cured
    if ((inoc_state_asexual[i] == Acute_asexual) || (inoc_state_asexual[i] == Chronic_asexual)) {
      
      // if due to become infective at a future timepoint then store this
      // timepoint
      bool due_acute_infective = (inoc_events[i].count(Event_begin_infective_acute) != 0);
      bool due_chronic_infective = (inoc_events[i].count(Event_begin_infective_chronic) != 0);
      int t_infective_start = 0;
      if (due_acute_infective) {
        t_infective_start = inoc_events[i][Event_begin_infective_acute];
      } else if (due_chronic_infective) {
        t_infective_start = inoc_events[i][Event_begin_infective_chronic];
      }
      
      // clear bloodstage infection
      inoc_state_asexual[i] = Inactive_asexual;
      
      // clear all scheduled events
      inoc_events[i].clear();
      
      // reinstate infectious start if needed
      if (due_acute_infective) {
        new_inoc_event(t_infective_start, Event_begin_infective_acute, i);
      } else if (due_chronic_infective) {
        new_inoc_event(t_infective_start, Event_begin_infective_chronic, i);
      }
      
      // schedule new infectious end
      new_inoc_event(t + params->g, Event_end_infective, i);
      
    }
    
    // anything that is scheduled to emerge from the liver during the
    // prophylactic period is cured immediately upon emergence
    if (inoc_state_asexual[i] == Liverstage_asexual) {
      
      // inoculations that will emerge to acute stage
      if ((inoc_events[i].count(Event_Eh_to_Ah) != 0) && (inoc_events[i][Event_Eh_to_Ah] <= t2)) {
        
        // store time at which due to emerge
        int t_emerge = inoc_events[i][Event_Eh_to_Ah];
        
        // clear all scheduled events
        inoc_events[i].clear();
        
        // reschedule progression directly to clearance
        new_inoc_event(t_emerge, Event_Eh_to_Ah, i);
        new_inoc_event(t_emerge, Event_Ah_to_Sh, i);
        new_inoc_event(t_emerge, Event_end_infective, i);
        
        // inoculations that will emerge to chronic stage
      } else if ((inoc_events[i].count(Event_Eh_to_Ch) != 0) && (inoc_events[i][Event_Eh_to_Ch] <= t2)) {
        
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
    
  }  // end i loop through inoculations
  
  // recalculate time of next event
  t_next_inoc_event = params->max_time + 1;
  for (int i = 0; i < params->max_inoculations; ++i) {
    for (const auto & x : inoc_events[i]) {
      
      // catch to ensure that events cannot predate current time
      if (x.second < t) {
        print("time =", t);
        print_inoc_state();
        print_inoc_events();
        Rcpp::stop("error in Host::treatment(), events found that predate current time");
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


// #############################################################################
// #                                                                           #
// #   INOCULATION-LEVEL EVENTS                                                #
// #                                                                           #
// #############################################################################


//------------------------------------------------
// new Eh inoculation
void Host::new_Eh(int this_slot, int t, int &next_inoc_ID) {
  
  // add inoculation ID and increment
  inoc_ID_vec[this_slot] = next_inoc_ID++;
  
  // mark inoculation as active
  inoc_active[this_slot] = true;
  
  // update inoculation state
  inoc_state_asexual[this_slot] = Liverstage_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
  // recalculate host state
  recalculate_host_state(t);
  
}

//------------------------------------------------
// move inoculation from Eh to Ah
void Host::Eh_to_Ah(int this_slot, int t) {
  
  // update inoculation state
  inoc_state_asexual[this_slot] = Acute_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
  // recalculate host state
  recalculate_host_state(t);
  
}

//------------------------------------------------
// move inoculation from Eh to Ch
void Host::Eh_to_Ch(int this_slot, int t) {
  
  // update inoculation state
  inoc_state_asexual[this_slot] = Chronic_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
  // recalculate host state
  recalculate_host_state(t);
  
}

//------------------------------------------------
// move inoculation from Ah to Ch
void Host::Ah_to_Ch(int this_slot, int t) {
  
  // update inoculation state
  inoc_state_asexual[this_slot] = Chronic_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
  // recalculate host state
  recalculate_host_state(t);
  
}

//------------------------------------------------
// move inoculation from Ah to Sh
void Host::Ah_to_Sh(int this_slot, int t) {
  
  // update inoculation state
  inoc_state_asexual[this_slot] = Inactive_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
  // recalculate host state
  recalculate_host_state(t);
  
}

//------------------------------------------------
// move inoculation from Ch to Sh
void Host::Ch_to_Sh(int this_slot, int t) {
  
  // update inoculation state
  inoc_state_asexual[this_slot] = Inactive_asexual;
  
  // store time
  inoc_time_asexual[this_slot] = t;
  
  // recalculate host state
  recalculate_host_state(t);
  
}

//------------------------------------------------
// begin acutely infective period
void Host::begin_infective_acute(int this_slot, int t) {
  
  // update inoculation state
  inoc_state_sexual[this_slot] = Acute_sexual;
  
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
  
  // get current inoculation state
  State_sexual current_state = inoc_state_sexual[this_slot];
  
  // update inoculation state
  inoc_state_sexual[this_slot] = Chronic_sexual;
  
  // store time
  inoc_time_sexual[this_slot] = t;
  
  // if newly infective then add to infectives list
  if (current_state == Inactive_sexual && get_n_infective() == 1) {
    (*host_infective_index_ptr)[deme].push_back(index);
  }
  
}

//------------------------------------------------
// end infective period
void Host::end_infective(int this_slot) {
  
  // update host state
  inoc_state_sexual[this_slot] = Inactive_sexual;
  inoc_active[this_slot] = false;
  
  // reset times
  inoc_time_asexual[this_slot] = 0;
  inoc_time_sexual[this_slot] = 0;
  
  // if no longer infective then drop from infectives list
  if (get_n_infective() == 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
}

//------------------------------------------------
// recalculate host state from inoculations
void Host::recalculate_host_state(int t) {
  
  // store host state before recalculating
  State_host state_before = host_state;
  
  // recalculate host state based on inoculation states
  if (prophylaxis_on) {
    host_state = Host_Ph;
  } else {
    host_state = Host_Sh;
    for (int i = 0; i < params->max_inoculations; ++i) {
      if ((host_state == Host_Sh) && (inoc_state_asexual[i] == Liverstage_asexual)) {
        host_state = Host_Eh;
      }
      if ((host_state == Host_Sh || host_state == Host_Eh) && (inoc_state_asexual[i] == Chronic_asexual)) {
        host_state = Host_Ch;
      }
      if ((host_state == Host_Sh || host_state == Host_Eh || host_state == Host_Ch) && (inoc_state_asexual[i] == Acute_asexual)) {
        host_state = Host_Ah;
        break;
      }
    }
  }
  
  // if host state has changed, make a note of the time
  if (host_state != state_before) {
    host_state_previous = state_before;
    time_host_state_change = t;
  }
  
}

// #############################################################################
// #                                                                           #
// #   GETTERS, SETTERS AND RANDOM DRAWS                                       #
// #                                                                           #
// #############################################################################

//------------------------------------------------
// get next free inoculation slot
int Host::get_free_inoc_slot() {
  
  // loop until find inactive slot
  int ret = 0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (!inoc_active[i]) {
      break;
    }
    ret++;
  }
  
  // error check for outside range
  if (ret == params->max_inoculations) {
    Rcpp::stop("could not find free inoculation slot");
  }
  
  return ret;
}

//------------------------------------------------
// get total number of active inoculations
int Host::get_n_active_inoc() {
  return sum_bool(inoc_active);
}

//------------------------------------------------
// return vector of inoc IDs taken from infective inoculations
vector<int> Host::get_inoc_ID_vec() {
  vector<int> ret;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_sexual[i] != Inactive_sexual) {
      ret.push_back(inoc_ID_vec[i]);
    }
  }
  return ret;
}

//------------------------------------------------
// get host state
State_host Host::get_host_state() {
  return host_state;
}

//------------------------------------------------
// get total number of liverstage inoculations
int Host::get_n_liverstage() {
  int ret = 0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_asexual[i] == Liverstage_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of acute bloodstage inoculations
int Host::get_n_bloodstage_acute() {
  int ret = 0;
  for (unsigned int i = 0; i < inoc_state_asexual.size(); ++i) {
    if (inoc_state_asexual[i] == Acute_asexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get total number of chronic bloodstage inoculations
int Host::get_n_bloodstage_chronic() {
  int ret = 0;
  for (unsigned int i = 0; i < inoc_state_asexual.size(); ++i) {
    if (inoc_state_asexual[i] == Chronic_asexual) {
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
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_sexual[i] == Acute_sexual || inoc_state_sexual[i] == Chronic_sexual) {
      ret ++;
    }
  }
  return ret;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  
  // situations in which zero chance of infection
  if ((get_n_active_inoc() == params->max_inoculations) || prophylaxis_on) {
    return 0.0;
  }
  
  // get probability from flexible distribution
  int tmp = (cumul_infections < params->n_prob_infection) ? cumul_infections : params->n_prob_infection - 1;
  return params->prob_infection[tmp];
}

//------------------------------------------------
// get current probability of going to acute infection
double Host::get_prob_acute() {
  int tmp = (cumul_inoculations < params->n_prob_acute) ? cumul_inoculations : params->n_prob_acute - 1;
  return params->prob_acute[tmp];
}

//------------------------------------------------
// get current probability of going to chronic from acute infection
double Host::get_prob_AC() {
  int tmp = (cumul_inoculations < params->n_prob_AC) ? cumul_inoculations : params->n_prob_AC - 1;
  return params->prob_AC[tmp];
}

//------------------------------------------------
// draw duration of acute disease
int Host::draw_duration_acute() {
  int tmp = (cumul_inoculations < params->n_duration_acute) ? cumul_inoculations : params->n_duration_acute - 1;
  return (*sampler_duration_acute_ptr)[tmp].draw();
}

//------------------------------------------------
// draw duration of chronic disease
int Host::draw_duration_chronic() {
  int tmp = (cumul_inoculations < params->n_duration_chronic) ? cumul_inoculations : params->n_duration_chronic - 1;
  return (*sampler_duration_chronic_ptr)[tmp].draw();
}

//------------------------------------------------
// draw time until treatment of acute disease
int Host::draw_time_treatment_acute() {
  int tmp = (cumul_inoculations < params->n_time_treatment_acute) ? cumul_inoculations : params->n_time_treatment_acute - 1;
  return (*sampler_time_treatment_acute_ptr)[tmp].draw();
}

//------------------------------------------------
// draw time until treatment of chronic disease
int Host::draw_time_treatment_chronic() {
  int tmp = (cumul_inoculations < params->n_time_treatment_chronic) ? cumul_inoculations : params->n_time_treatment_chronic - 1;
  return (*sampler_time_treatment_chronic_ptr)[tmp].draw();
}

//------------------------------------------------
// draw duration of prophylaxis
int Host::draw_duration_prophylaxis() {
  int tmp = (cumul_inoculations < params->n_duration_prophylactic) ? cumul_inoculations : params->n_duration_prophylactic - 1;
  return (*sampler_duration_prophylatic_ptr)[tmp].draw();
}

//------------------------------------------------
// get current detectability by microscopy over acute inoculations
double Host::get_detectability_microscopy_acute(int t) {
  
  // acute detectability is the max detectability over all acute inoculations
  double ret = 0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_asexual[i] == Acute_asexual) {
      
      // get time since entered acute state
      int time_diff = t - inoc_time_asexual[i];
      
      // get detectability from appropriate distribution
      int tmp1 = (cumul_inoculations < params->n_detectability_microscopy_acute) ? cumul_inoculations : params->n_detectability_microscopy_acute - 1;
      int tmp2 = (time_diff < int(params->detectability_microscopy_acute[tmp1].size())) ? time_diff : int(params->detectability_microscopy_acute[tmp1].size()) - 1;
      double inoc_detectability = params->detectability_microscopy_acute[tmp1][tmp2];
      
      // ret is max of ret and inoc_detectability
      ret = (ret > inoc_detectability) ? ret : inoc_detectability;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get current detectability by microscopy over chronic inoculations
double Host::get_detectability_microscopy_chronic(int t) {
  
  // chronic detectability is the max detectability over all chronic inoculations
  double ret = 0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_asexual[i] == Chronic_asexual) {
      
      // get time since entered chronic state
      int time_diff = t - inoc_time_asexual[i];
      
      // get detectability from appropriate distribution
      int tmp1 = (cumul_inoculations < params->n_detectability_microscopy_chronic) ? cumul_inoculations : params->n_detectability_microscopy_chronic - 1;
      int tmp2 = (time_diff < int(params->detectability_microscopy_chronic[tmp1].size())) ? time_diff : int(params->detectability_microscopy_chronic[tmp1].size()) - 1;
      double inoc_detectability = params->detectability_microscopy_chronic[tmp1][tmp2];
      
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
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_asexual[i] == Acute_asexual) {
      
      // get time since entered acute state
      int time_diff = t - inoc_time_asexual[i];
      
      // get detectability from appropriate distribution
      int tmp1 = (cumul_inoculations < params->n_detectability_PCR_acute) ? cumul_inoculations : params->n_detectability_PCR_acute - 1;
      int tmp2 = (time_diff < int(params->detectability_PCR_acute[tmp1].size())) ? time_diff : int(params->detectability_PCR_acute[tmp1].size()) - 1;
      double inoc_detectability = params->detectability_PCR_acute[tmp1][tmp2];
      
      // ret is max of ret and inoc_detectability
      ret = (ret > inoc_detectability) ? ret : inoc_detectability;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get current detectability by PCR over chronic inoculations
double Host::get_detectability_PCR_chronic(int t) {
  
  // chronic detectability is the max detectability over all chronic inoculations
  double ret = 0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc_state_asexual[i] == Chronic_asexual) {
      
      // get time since entered chronic state
      int time_diff = t - inoc_time_asexual[i];
      
      // get detectability from appropriate distribution
      int tmp1 = (cumul_inoculations < params->n_detectability_PCR_chronic) ? cumul_inoculations : params->n_detectability_PCR_chronic - 1;
      int tmp2 = (time_diff < int(params->detectability_PCR_chronic[tmp1].size())) ? time_diff : int(params->detectability_PCR_chronic[tmp1].size()) - 1;
      double inoc_detectability = params->detectability_PCR_chronic[tmp1][tmp2];
      
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
  for (int i = 0; i < params->max_inoculations; ++i) {
    if ((inoc_state_sexual[i] == Acute_sexual) || (inoc_state_sexual[i] == Chronic_sexual)) {
      
      // get time since became infective
      int time_diff = t - inoc_time_sexual[i];
      
      // split by acute vs. chronic infectivity
      double inoc_infectivity;
      if (inoc_state_sexual[i] == Acute_sexual) {
        
        // get infectivity from appropriate distribution
        int tmp1 = (cumul_inoculations < params->n_infectivity_acute) ? cumul_inoculations : params->n_infectivity_acute-1;
        int tmp2 = (time_diff < int(params->infectivity_acute[tmp1].size())) ? time_diff : int(params->infectivity_acute[tmp1].size())-1;
        inoc_infectivity = params->infectivity_acute[tmp1][tmp2];
        
      } else {
        
        // get infectivity from appropriate distribution
        int tmp1 = (cumul_inoculations < params->n_infectivity_chronic) ? cumul_inoculations : params->n_infectivity_chronic-1;
        int tmp2 = (time_diff < int(params->infectivity_chronic[tmp1].size())) ? time_diff : int(params->infectivity_chronic[tmp1].size())-1;
        inoc_infectivity = params->infectivity_chronic[tmp1][tmp2];
        
      }
      
      // ret is max of ret and inoc_infectivity
      ret = (ret > inoc_infectivity) ? ret : inoc_infectivity;
    }
  }
  
  return ret;
}

//------------------------------------------------
// get host age in years at time t
int Host::get_age(int t) {
  return (t - birth_day) / 365;
}


// #############################################################################
// #                                                                           #
// #   OUTPUTS                                                                 #
// #                                                                           #
// #############################################################################

//------------------------------------------------
// update numerator and denominator of output elements
void Host::update_output(Measure measure, Model_state state, Diagnostic diagnostic, int t, double &numer, double &denom) {
  
  // counts or prevalence
  if ((measure == Measure_count) || (measure == Measure_prevalence)) {
    
    // update numerator
    if ( (state == Model_A) && (get_host_state() == Host_Ah) ) {
      if (diagnostic == Diagnostic_true) {
        numer += 1.0;
      } else if (diagnostic == Diagnostic_microscopy) {
        numer += get_detectability_microscopy_acute(t);
      } else if (diagnostic == Diagnostic_PCR) {
        numer += get_detectability_PCR_acute(t);
      }
    } else if ( (state == Model_C) && (get_host_state() == Host_Ch) ) {
      if (diagnostic == Diagnostic_true) {
        numer += 1.0;
      } else if (diagnostic == Diagnostic_microscopy) {
        numer += get_detectability_microscopy_chronic(t);
      } else if (diagnostic == Diagnostic_PCR) {
        numer += get_detectability_PCR_chronic(t);
      }
    } else if ( ((state == Model_S) && (get_host_state() == Host_Sh)) ||
         ((state == Model_E) && (get_host_state() == Host_Eh)) ||
         ((state == Model_P) && (get_host_state() == Host_Ph)) ||
         (state == Model_H) ) {
      numer += 1.0;
    }
    
    // update denominator
    if (measure == Measure_prevalence) {
      denom += 1.0;
    }
    
  // incidence
  } else if (measure == Measure_incidence) {
    
    if (state == Model_E) {
      if ((host_state == Host_Eh) && (host_state_previous != Host_Eh) && (time_host_state_change == t)) {
        numer += 1.0;
      }
      if (host_state != Host_Eh) {
        denom += 1.0;
      }
    } else if (state == Model_A) {
      if ((host_state == Host_Ah) && (host_state_previous != Host_Ah) && (time_host_state_change == t)) {
        numer += 1.0;
      }
      if (host_state != Host_Ah) {
        denom += 1.0;
      }
    }  else if (state == Model_C) {
      if ((host_state == Host_Ch) && (host_state_previous != Host_Ch) && (time_host_state_change == t)) {
        numer += 1.0;
      }
      if (host_state != Host_Ch) {
        denom += 1.0;
      }
    } else if (state == Model_P) {
      if ((host_state == Host_Ph) && (host_state_previous != Host_Ph) && (time_host_state_change == t)) {
        numer += 1.0;
      }
      if (host_state != Host_Ph) {
        denom += 1.0;
      }
    }
    
  }
  
}


// #############################################################################
// #                                                                           #
// #   DIAGNOSTICS AND CHECKS                                                  #
// #                                                                           #
// #############################################################################

//------------------------------------------------
// print innoc_events
void Host::print_inoc_events() {
  
  for (int i = 0; i < params->max_inoculations; ++i) {
    for (const auto & x : inoc_events[i]) {
      Rcpp::Rcout << "[" << x.first << ", " << x.second << "] ";
    }
    Rcpp::Rcout << "\n";
  }
}

//------------------------------------------------
// print state of inoc slots
void Host::print_inoc_state() {
  Rcpp::Rcout << "host_ID: " << host_ID << "\n";
  Rcpp::Rcout << "inoc_ID: ";
  print_vector(inoc_ID_vec);
  Rcpp::Rcout << "active : ";
  print_vector(inoc_active);
  Rcpp::Rcout << "asex   : ";
  print_vector(inoc_state_asexual);
  Rcpp::Rcout << "sex    : ";
  print_vector(inoc_state_sexual);
  Rcpp::Rcout << "time   : ";
  print_vector(inoc_time_sexual);
}

//------------------------------------------------
// if host carries any infective inoculations, check that host index is present
// in host_infective_index
void Host::check_host_infective_index(int x) {
  
  if ((get_n_infective() != 0) && (!is_in_vector(index, (*host_infective_index_ptr)[deme]))) {
    print(x);
    Rcpp::stop("error in check_host_infective_index() function");
  }
  
}
