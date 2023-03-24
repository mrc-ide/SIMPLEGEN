
#include "Host.hpp"

using namespace std;

// link to globals
extern int next_host_ID;
extern int next_inoc_ID;

//------------------------------------------------
// initialise host
void Host::init(Parameters &params_, std::ofstream &transmission_record_, int deme_) {
  
  // store pointers
  params = &params_;
  transmission_record = &transmission_record_;
  
  // identifiers
  host_ID = next_host_ID++;    // unique ID of this host, modified upon death
  home_deme = deme_;           // deme into which this host was born
  deme = home_deme;            // deme in which this host currently resides
  
  // host state now, previously, and time at which it changed
  host_state = Host_Sh;
  host_state_previous = Host_Sh;
  time_host_state_change = 0;
  
  // draw birth and death days from stable demography distribution
  draw_starting_age();
  
  // treatment status
  time_treatment = -1;
  prophylaxis_on = false;
  time_prophylaxis_stop = -1;
  
  // cumulative count of how many times this host has been infected
  cumul_infections = 0;
  
  // draw host-specific treatment seeking probability
  draw_treatment_seeking();
  
  // initialise inoculation slots
  any_inoc_active = false;
  inoc = vector<Host_inoc>(params->max_inoculations);
}

//------------------------------------------------
// draw birth and death days given a known stable demography distribution. Only
// needed at start of simulation, after which hosts that die are instantly
// re-born and have death time drawn from entire life table
void Host::draw_starting_age() {
  
  // draw age from stable demography distribution
  int age_years = params->draw_age_stable();
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
  double prop_year_remaining = 1.0 - extra_days / 365.0; // (x in the above derivation)
  double prob_die_this_year = 1.0 - pow(1.0 - params->life_table[age_years], 1.0 - prop_year_remaining);
  if (rbernoulli1(prob_die_this_year)) {
    life_days = age_years*365 + sample2(extra_days, 364);
  } else {
    // if we do not die in the current year of age then loop through all
    // remaining years, performing Bernoulli draw from probability of dying in
    // that year
    for (size_t i = (age_years + 1); i < params->life_table.size(); ++i) {
      if (rbernoulli1(params->life_table[i])) {
        life_days = i*365 + sample2(0, 364);
        break;
      }
    }
  }
  
  // calculate final time of birth and death from age at time 0 and life_days
  time_birth = -age_days;
  time_death = life_days - age_days;
  
  // cannot die at time 0 otherwise will never be prompted to die again
  time_death = (time_death == 0) ? 1 : time_death;
}

//------------------------------------------------
// draw host-specific treatment seeking parameter, either from Beta distribution
// or Dirac delta distribution (fixed value)
void Host::draw_treatment_seeking() {
  
  if ((params->treatment_seeking_mean == 0) || (params->treatment_seeking_mean == 1) || (params->treatment_seeking_sd == 0)) {
    treatment_seeking_prob = params->treatment_seeking_mean;
  } else {
    double m = params->treatment_seeking_mean;
    double s = params->treatment_seeking_sd;
    double alpha = m*m*(1 - m)/(s*s) - m;
    double beta = m*(1 - m)*(1 - m)/(s*s) - (1 - m);
    treatment_seeking_prob = rbeta1(alpha, beta);
  }
}

//------------------------------------------------
// check for any events at this point in time. Note that some events e.g.
// migration or death are triggered at the Host_pop level as they require a
// change to the number of hosts per deme, and so these are not present here
void Host::update(int t) {
  
  // keep track of whether any updates are applied
  bool any_updates = false;
  
  // TODO - treatment
  //if (time_treatment == t) {
  //  treatment(t);
  //}
  
  // end of prophylactic period
  if (time_prophylaxis_stop == t) {
    any_updates = true;
    end_prophylaxis(t);
  }
  
  // inoculation-level events
  if (any_inoc_active) {
    for (int i = 0; i < params->max_inoculations; ++i) {
      if (inoc[i].active) {
        
        // start acute bloodstage. Clear everything if in prophylactic state,
        // i.e. never become infectious
        if (inoc[i].time_start_acute == t) {
          any_updates = true;
          if (prophylaxis_on) {
            inoc[i].reset();
          } else {
            inoc[i].state_asexual = Acute_asexual;
          }
        }
        
        // stop acute bloodstage. Transition to chronic or clear bloodstage
        if (inoc[i].time_stop_acute == t) {
          any_updates = true;
          //print(t, "stop acute bloodstage");
          if (inoc[i].transition_AC) {
            inoc[i].state_asexual = Chronic_asexual;
          } else {
            inoc[i].state_asexual = Inactive_asexual;
          }
        }
        
        // start chronic bloodstage. Clear everything if in prophylactic state,
        // i.e. never become infectious
        if (inoc[i].time_start_chronic == t) {
          any_updates = true;
          //print(t, "start chronic bloodstage");
          if (prophylaxis_on) {
            inoc[i].reset();
          } else {
            inoc[i].state_asexual = Chronic_asexual;
          }
        }
        
        // stop chronic bloodstage
        if (inoc[i].time_stop_chronic == t) {
          any_updates = true;
          //print(t, "stop chronic bloodstage");
          inoc[i].state_asexual = Inactive_asexual;
        }
        
        // start acute sexual stage
        if ((inoc[i].time_start_acute != -1) && ((inoc[i].time_start_acute + params->g) == t)) {
          any_updates = true;
          //print(t, "start acute sexual");
          inoc[i].state_sexual = Acute_sexual;
        }
        
        // stop acute sexual stage. Transition to chronic or clear
        if ((inoc[i].time_stop_acute != -1) && ((inoc[i].time_stop_acute + params->g) == t)) {
          any_updates = true;
          //print(t, "stop acute sexual");
          if (inoc[i].transition_AC) {
            inoc[i].state_sexual = Chronic_sexual;
          } else {
            inoc[i].reset();
          }
        }
        
        // start chronic sexual stage
        if ((inoc[i].time_start_chronic != -1) && ((inoc[i].time_start_chronic + params->g) == t)) {
          any_updates = true;
          //print(t, "start chronic sexual");
          inoc[i].state_sexual = Chronic_sexual;
        }
        
        // stop chronic sexual stage
        if ((inoc[i].time_stop_chronic != -1) && ((inoc[i].time_stop_chronic + params->g) == t)) {
          any_updates = true;
          //print(t, "stop chronic sexual");
          inoc[i].reset();
        }
        
      }
    }
  }
  
  // some actions if there have been any updates at this time
  if (any_updates) {
    
    // update any_inoc_active, as some may have become inactive due to complete
    // recovery from asexual and sexual stages
    any_inoc_active = get_any_inoc_active();
    
    // recalculate host state
    recalc_host_state(t);
  }
}

//------------------------------------------------
// migration
void Host::migrate(int new_deme) {
  deme = new_deme;
}

//------------------------------------------------
// death and simultaneous rebirth
void Host::death(int t) {
  
  // new unique ID
  host_ID = next_host_ID++;
  
  // reset host state
  host_state = Host_Sh;
  host_state_previous = Host_Sh;
  time_host_state_change = 0;
  
  // reborn into home deme
  migrate(home_deme);
  
  // reborn today. Draw life duration from demography distribution. If die on
  // same day as born then immediately resample
  time_birth = t;
  int life_years = 0;
  int life_days = 0;
  while (life_days == 0) {
    life_years = params->draw_age_death();
    life_days = life_years*365 + sample2(0, 364);
  }
  time_death = time_birth + life_days;
  
  // reset treatment status
  time_treatment = 0;
  prophylaxis_on = false;
  time_prophylaxis_stop = 0;
  
  // reset cumulative counts
  cumul_infections = 0;
  
  // re-draw host-specific treatment seeking probability
  draw_treatment_seeking();
  
  // reset inoculation slots
  any_inoc_active = false;
  for (int i = 0; i < params->max_inoculations; ++i) {
    inoc[i].reset();
  }
}

//------------------------------------------------
// treatment clears all bloodstage infections and curtails sexual stages
void Host::treatment(int t) {
  
  // activate prophylactic period and draw duration
  prophylaxis_on = true;
  time_prophylaxis_stop = t + params->draw_duration_prophylactic(cumul_infections);
  
  // update all inoculations as needed
  for (int i = 0; i < params->max_inoculations; ++i) {
    
    // clear active bloodstage inoculations, which will curtail corresponding
    // sexual stages. NB, this leaves liverstage inoculations alone, as
    // liverstage inoculations that emerge during the prophylactic period will
    // be dealt with when they emerge
    if (inoc[i].state_asexual == Acute_asexual) {
      inoc[i].state_asexual = Inactive_asexual;
      inoc[i].time_stop_acute = t;
      inoc[i].time_start_chronic = -1;
      inoc[i].time_stop_chronic = -1;
    } else if (inoc[i].state_asexual == Chronic_asexual) {
      inoc[i].state_asexual = Inactive_asexual;
      inoc[i].time_stop_chronic = t;
    }
  }
  
  // recalculate host state
  recalc_host_state(t);
}

//------------------------------------------------
// terminate prophylactic period
void Host::end_prophylaxis(int t) {
  prophylaxis_on = false;
  time_prophylaxis_stop = -1;
  
  // recalculate host state
  recalc_host_state(t);
}

//------------------------------------------------
// new infection from mosquito
void Host::infect(int t, int mosquito_ID, int mosquito_inoc_ID) {
  
  // get next free inoc slot, or return if no free slots
  int slot = -1;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (!inoc[i].active) {
      slot = i;
      break;
    }
  }
  if (slot == -1) {
    return;
  }
  
  // print to transmission record
  if (params->save_transmission_record) {
    *transmission_record << t << ",H," << host_ID << "," << mosquito_ID << "," << next_inoc_ID << "," << mosquito_inoc_ID << "," << deme + 1 << "\n";
  }
  
  // enter liverstage and schedule all future events relating to this infection
  any_inoc_active = true;
  inoc[slot].ID = next_inoc_ID++;
  inoc[slot].active = true;
  inoc[slot].state_asexual = Liverstage_asexual;
  double prob_acute = params->get_prob_acute(cumul_infections);
  if (rbernoulli1(prob_acute)) {
    int dur_acute = params->draw_duration_acute(cumul_infections);
    inoc[slot].time_start_acute = t + params->u;
    inoc[slot].time_stop_acute = t + params->u + dur_acute;
    double prob_AC = params->get_prob_AC(cumul_infections);
    inoc[slot].transition_AC = rbernoulli1(prob_AC);
    if (inoc[slot].transition_AC) {
      int dur_chronic = params->draw_duration_chronic(cumul_infections);
      inoc[slot].time_start_chronic = t + params->u + dur_acute;
      inoc[slot].time_stop_chronic = t + params->u + dur_acute + dur_chronic;
    }
  } else {
    int dur_chronic = params->draw_duration_chronic(cumul_infections);
    inoc[slot].time_start_chronic = t + params->u;
    inoc[slot].time_stop_chronic = t + params->u + dur_chronic;
  }
  
  // increment cumulative infection count
  cumul_infections++;
  
  // recalculate host state
  recalc_host_state(t);
}

//------------------------------------------------
// get probability of becoming infected given infectious bite
double Host::get_prob_infection(int t) {
  return params->get_prob_infection(cumul_infections);
}

//------------------------------------------------
// get onward infectivity at time t
double Host::get_infectivity(int t) {
  if (!any_inoc_active) {
    return 0.0;
  }
  // infectivity is the max infectivity over all inoculations
  double ret = 0.0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc[i].state_sexual == Acute_sexual) {
      double infectivity_acute = params->get_infectivity_acute(t, inoc[i].time_start_acute + params->g, cumul_infections);
      ret = max(ret, infectivity_acute);
    } else if (inoc[i].state_sexual == Chronic_sexual) {
      double infectivity_chronic = params->get_infectivity_chronic(t, inoc[i].time_start_chronic + params->g, cumul_infections);
      ret = max(ret, infectivity_chronic);
    }
  }
  return ret;
}

//------------------------------------------------
// get detectability by microscopy at time t
double Host::get_detectability_microscopy(int t) {
  if (!any_inoc_active) {
    return 0.0;
  }
  // detectability is the max detectability over all inoculations
  double ret = 0.0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc[i].state_asexual == Acute_asexual) {
      double detectability_acute = params->get_detectability_microscopy_acute(t, inoc[i].time_start_acute, cumul_infections);
      ret = max(ret, detectability_acute);
    } else if (inoc[i].state_asexual == Chronic_asexual) {
      double detectability_chronic = params->get_detectability_microscopy_chronic(t, inoc[i].time_start_chronic, cumul_infections);
      ret = max(ret, detectability_chronic);
    }
  }
  return ret;
}

//------------------------------------------------
// get detectability by PCR at time t
double Host::get_detectability_PCR(int t) {
  if (!any_inoc_active) {
    return 0.0;
  }
  // detectability is the max detectability over all inoculations
  double ret = 0.0;
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc[i].state_asexual == Acute_asexual) {
      double detectability_acute = params->get_detectability_PCR_acute(t, inoc[i].time_start_acute, cumul_infections);
      ret = max(ret, detectability_acute);
    } else if (inoc[i].state_asexual == Chronic_asexual) {
      double detectability_chronic = params->get_detectability_PCR_chronic(t, inoc[i].time_start_chronic, cumul_infections);
      ret = max(ret, detectability_chronic);
    }
  }
  return ret;
}

//------------------------------------------------
// return host age in years
int Host::get_age_years(int t) {
  return (t - time_birth) / 365;
}

//------------------------------------------------
// return if there are any active inoculations
bool Host::get_any_inoc_active() {
  for (int i = 0; i < params->max_inoculations; ++i) {
    if (inoc[i].active) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------
// return current host state
State_host Host::get_host_state() {
  return host_state;
}

//------------------------------------------------
// return previous host state
State_host Host::get_host_state_previous() {
  return host_state_previous;
}

//------------------------------------------------
// return time at which host state changed most recently
int Host::get_time_host_state_change() {
  return time_host_state_change;
}

//------------------------------------------------
// return time of treatment
int Host::get_time_treatment() {
  return time_treatment;
}

//------------------------------------------------
// recalculate host state based on inoculation states etc.
void Host::recalc_host_state(int t) {
  
  // go through scale: prophylactic, acute, chronic, exposed, susceptible
  State_host new_state = Host_Sh;
  if (prophylaxis_on) {
    new_state = Host_Ph;
  } else {
    for (int i = 0; i < params->max_inoculations; ++i) {
      if (!inoc[i].active) {
        continue;
      }
      if ((new_state == Host_Sh) && (inoc[i].state_asexual == Liverstage_asexual)) {
        new_state = Host_Eh;
      }
      if ((new_state == Host_Sh || new_state == Host_Eh) && (inoc[i].state_asexual == Chronic_asexual)) {
        new_state = Host_Ch;
      }
      if ((new_state == Host_Sh || new_state == Host_Eh || new_state == Host_Ch) && (inoc[i].state_asexual == Acute_asexual)) {
        new_state = Host_Ah;
        break;
      }
    }
  }
  
  // if host state has changed, make a note of the time
  if (new_state != host_state) {
    host_state_previous = host_state;
    host_state = new_state;
    time_host_state_change = t;
  }
}

//------------------------------------------------
// print summary
void Host::print_summary() {
  
  print("---------------------------");
  
  print("ID:", host_ID);
  print("home_deme:", home_deme);
  print("deme:", deme);
  
  print("time_birth:", time_birth);
  print("time_death:", time_death);
  
  print("time_treatment:", time_treatment);
  print("prophylaxis_on:", prophylaxis_on);
  print("time_prophylaxis_stop:", time_prophylaxis_stop);
  
  print("cumul_infections:", cumul_infections);
  print("treatment_seeking_prob:", treatment_seeking_prob);
  
  // inoculation level values
  print("");
  
  print("inoc.active:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].active << " ";
  }
  print("");
  
  print("inoc.ID:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].ID << " ";
  }
  print("");
  
  print("inoc.state_asexual:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].state_asexual << " ";
  }
  print("");
  
  print("inoc.state_sexual:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].state_sexual << " ";
  }
  print("");
  
  print("inoc.time_start_acute:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].time_start_acute << " ";
  }
  print("");
  
  print("inoc.time_stop_acute:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].time_stop_acute << " ";
  }
  print("");
  
  print("inoc.time_start_chronic:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].time_start_chronic << " ";
  }
  print("");
  
  print("inoc.time_stop_chronic:");
  for (int i = 0; i < params->max_inoculations; ++i) {
    std::cout << inoc[i].time_stop_chronic << " ";
  }
  print("");
  
  print("");
}
