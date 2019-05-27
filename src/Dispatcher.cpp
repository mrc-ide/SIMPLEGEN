
#include "Dispatcher.h"
#include "probability_v2.h"

using namespace std;

//------------------------------------------------
// constructor
Dispatcher::Dispatcher() {
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(age_stable, 1000);
  sampler_age_death = Sampler(age_death, 1000);
  sampler_duration_acute = Sampler(duration_acute[0], 1000);
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  schedule_death = vector<set<int>>(max_time);
  schedule_Eh_to_Ih = vector<vector<pair<int, int>>>(max_time);
  schedule_Ih_to_Sh = vector<vector<pair<int, int>>>(max_time);
  schedule_infective = vector<vector<pair<int, int>>>(max_time);
  schedule_infective_recovery = vector<vector<pair<int, int>>>(max_time);
  
  // counts of host types
  H = H_init;
  Sh = H;
  Eh = vector<int>(n_demes);
  Ih = vector<int>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single population we can
  // simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(sum(H));
  next_host_ID = 0;
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only
  host_index = vector<vector<int>>(n_demes);
  host_infective_index = vector<vector<int>>(n_demes);
  int tmp1 = 0;
  for (int k = 0; k < n_demes; ++k) {
    host_index[k] = seq_int(tmp1, tmp1+H[k]-1);
    reshuffle(host_index[k]);
    tmp1 += H[k];
  }
  
  // initialise hosts
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      int this_host = host_index[k][i];
      host_pop[this_host].init(this_host, next_host_ID, k,
                               Sh, Eh, Ih,
                               host_infective_index,
                               schedule_death, schedule_Eh_to_Ih, schedule_Ih_to_Sh,
                               schedule_infective, schedule_infective_recovery,
                               sampler_age_stable, sampler_age_death, sampler_duration_acute);
    }
  }
  
  // seed initial infections
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < seed_infections[k]; ++i) {
      host_pop[host_index[k][i]].denovo_infection(0);
    }
  }
  
  // counts of mosquito types
  M_total = sum(M);
  Sv = M;
  Ev = vector<int>(n_demes);
  Iv = vector<int>(n_demes);
  
  // objects for tracking mosquitoes that die in lag phase
  Ev_death = vector<vector<int>>(n_demes, vector<int>(v));
  
  // populations of mosquitoes at various stages
  Ev_pop = vector<vector<vector<Mosquito>>>(n_demes, vector<vector<Mosquito>>(v));
  Iv_pop = vector<vector<Mosquito>>(n_demes);
  
  // objects for storing daily values:
  // 0 = Sh, 1 = Eh, 2 = Ih, 3 = Sv, 4 = Ev, 5 = Iv, 6 = EIR
  daily_values = vector<vector<vector<double>>>(n_demes, vector<vector<double>>(max_time, vector<double>(7)));
  
  // misc
  EIR = vector<double>(n_demes);
  
}

//------------------------------------------------
// run main simulation
void Dispatcher::run_simulation() {
  
  // start message
  if (!silent) {
    print("Running simulation");
  }
  
  // initialise indices
  int v_ringbuffer = 0; // ring buffer that loops back to 0 when it exceeds (v-1)
  
  // loop through daily time steps
  for (int t = 0; t < max_time; ++t) {
    
    // update ring buffer index
    v_ringbuffer = (v_ringbuffer == v-1) ? 0 : v_ringbuffer + 1;
    
    // skip over first iteration to ensure user-defined values appear first in
    // results
    if (t != 0) {
      
      //-------- SCHEDULED HUMAN EVENTS --------
      
      // scheduled deaths
      for (auto it = schedule_death[t].begin(); it != schedule_death[t].end(); ++it) {
        int this_host = *it;
        host_pop[this_host].death(next_host_ID, t);
      }
      
      // scheduled Eh to Ih
      for (auto it = schedule_Eh_to_Ih[t].begin(); it != schedule_Eh_to_Ih[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Eh_to_Ih(this_slot);
      }
      
      // scheduled Ih to Sh
      for (auto it = schedule_Ih_to_Sh[t].begin(); it != schedule_Ih_to_Sh[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Ih_to_Sh(this_slot);
      }
      
      // scheduled become infective
      for (auto it = schedule_infective[t].begin(); it != schedule_infective[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].begin_infective(this_slot, t);
      }
      
      // scheduled infective recovery
      for (auto it = schedule_infective_recovery[t].begin(); it != schedule_infective_recovery[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].end_infective(this_slot);
      }
      
      
      // loop through demes
      for (int k = 0; k < n_demes; ++k) {
        
        //-------- NEW HUMAN EVENTS --------
        
        // get number of new infectious bites on humans
        EIR[k] = a*Iv[k]/double(H[k]);
        double prob_infectious_bite = 1 - exp(-EIR[k]);                // probability of new infectious bite on host
        int n_infectious_bites = rbinom1(H[k], prob_infectious_bite);  // total number of new infectious bites
        
        // apply new infectious bites
        for (int i = 0; i < n_infectious_bites; ++i) {
          
          // choose host at random
          int rnd1 = sample2(0, H[k]-1);
          int this_host = host_index[k][rnd1];
          
          // determine whether infectious bite is successful
          if (rbernoulli1(host_pop[this_host].get_prob_infection())) {
            
            // infect host
            host_pop[this_host].infection(t);
          }
          
        }  // end loop over infectious bites
        
        
        //-------- SCHEDULED MOSQUITO EVENTS --------
        
        // deaths in Ev
        Sv[k] += Ev_death[k][v_ringbuffer];
        Ev[k] -= Ev_death[k][v_ringbuffer];
        Ev_death[k][v_ringbuffer] = 0;
        
        // move Ev into Iv
        int delta_Ev = int(Ev_pop[k][v_ringbuffer].size());
        if (delta_Ev > 0) {
          Ev[k] -= delta_Ev;
          Iv[k] += delta_Ev;
          push_back_multiple(Iv_pop[k], Ev_pop[k][v_ringbuffer]);
          Ev_pop[k][v_ringbuffer].clear();
        }
        
        
        //-------- NEW MOSQUITO EVENTS --------
        
        // rate of mosquito biting infective host
        double rate_bite_infective = a*host_infective_index[k].size()/double(H[k]); 
        
        // draw number of mosquitoes that bite infective host or die (competing
        // hazards)
        double prob_bite_infective_or_death = 1 - exp(-(rate_bite_infective + mu));
        int n_bite_infective_or_death = rbinom1(Sv[k], prob_bite_infective_or_death);
        
        // draw number of mosquitoes that bite infective host, rather than dying.
        double relative_prob_bite_infective = rate_bite_infective/(rate_bite_infective + mu);
        int n_bite_infective = rbinom1(n_bite_infective_or_death, relative_prob_bite_infective);
        
        // if the infectivity of hosts towards mosquitoes was constant then we
        // could also filter out all infections that do not take hold in
        // mosquitoes at this early stage. However, in our model the infectivity
        // can vary with time and host characteristics. Therefore, instead
        // calculate the *maximum* possible infectivity, and filter based on
        // this number. We then need a second level of filtering to account for
        // the actual values.
        //
        // For example, if the infectivities of acute vs. chronicly infected
        // hosts are {0.1, 0.05}, and constant over time for the sake of
        // simplicity, then filter based on the value 0.1, i.e. draw the number
        // of infected mosquitoes N from binomial with probability 0.1. Then
        // loop through all N infections and draw from true probability of
        // infection, which is Bernoulli with probability {1.0, 0.5}. This
        // method should throw out a large number of infections that do not take
        // hold at an early stage, thereby reducing the number we need to loop
        // over.
        int n_query_infection = rbinom1(n_bite_infective, max_infectivity);
        
        // loop through query infective bites
        for (int i = 0; i < n_query_infection; ++i) {
          
          // choose host at random from infectives
          int rnd1 = sample2(0, host_infective_index[k].size()-1);
          int this_host = host_infective_index[k][rnd1];
          
          // get infectivity and draw whether infection takes hold in mosquito
          double host_infectivity = host_pop[this_host].get_infectivity(t);
          if (rbernoulli1(host_infectivity/max_infectivity)) {
            
            // update deme counts
            Sv[k]--;
            Ev[k]++;
            
            // the majority of new mosquito infections will die in lag phase.
            // Schedule these deaths to move back into Sv in future steps.
            // Otherwise add to Ev_pop
            int mosq_time_death = rgeom1(prob_mosq_death) + 1;
            if (mosq_time_death <= v) {
              
              // schedule death for future time
              Ev_death[k][(v_ringbuffer + mosq_time_death) % v]++;
              
            } else {
              
              // add to Ev_pop, to enter Iv_pop at future time
              int this_host_ID = host_pop[this_host].ID;
              Ev_pop[k][v_ringbuffer].emplace_back(this_host_ID, t);
              
            }
            
          }
        } // end loop through query infective bites
        
        // deaths in Iv
        int death_Iv = rbinom1(Iv[k], prob_mosq_death);
        Sv[k] += death_Iv;
        Iv[k] -= death_Iv;
        for (int i = 0; i < death_Iv; ++i) {
          int rnd1 = sample2(0, Iv[k]-1);
          quick_erase(Iv_pop[k], rnd1);
        }
        
      }  // end loop over demes
      
    }  // end skip over first time step
    
    
    //-------- STORE RESULTS --------
    
    // store daily values
    for (int k = 0; k < n_demes; ++k) {
      daily_values[k][t] = {double(Sh[k]), double(Eh[k]), double(Ih[k]),
                              double(Sv[k]), double(Ev[k]), double(Iv[k]),
                              EIR[k]};
    }
    
  }  // end loop through daily time steps
  
}

