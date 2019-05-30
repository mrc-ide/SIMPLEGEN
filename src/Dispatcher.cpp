
#include "Dispatcher.h"
#include "probability_v2.h"

using namespace std;

//------------------------------------------------
// constructor
Dispatcher::Dispatcher() {
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(age_stable, 1000);
  sampler_age_death = Sampler(age_death, 1000);
  sampler_duration_acute = vector<Sampler>(n_duration_acute);
  for (int i = 0; i < n_duration_acute; ++i) {
    sampler_duration_acute[i] = Sampler(duration_acute[i], 1000);
  }
  sampler_duration_chronic = vector<Sampler>(n_duration_chronic);
  for (int i = 0; i < n_duration_chronic; ++i) {
    sampler_duration_chronic[i] = Sampler(duration_chronic[i], 1000);
  }
  sampler_time_treatment_acute = vector<Sampler>(n_time_treatment_acute);
  for (int i = 0; i < n_time_treatment_acute; ++i) {
    sampler_time_treatment_acute[i] = Sampler(time_treatment_acute[i], 1000);
  }
  sampler_time_treatment_chronic = vector<Sampler>(n_time_treatment_chronic);
  for (int i = 0; i < n_time_treatment_chronic; ++i) {
    sampler_time_treatment_chronic[i] = Sampler(time_treatment_chronic[i], 1000);
  }
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  schedule_death = vector<set<int>>(max_time);
  schedule_Eh_to_Ah = vector<vector<pair<int, int>>>(max_time);
  schedule_Eh_to_Ch = vector<vector<pair<int, int>>>(max_time);
  schedule_Ah_to_Ch = vector<vector<pair<int, int>>>(max_time);
  schedule_Ah_to_Sh = vector<vector<pair<int, int>>>(max_time);
  schedule_Ch_to_Sh = vector<vector<pair<int, int>>>(max_time);
  schedule_Ah_to_Ph = vector<set<int>>(max_time);
  schedule_Ch_to_Ph = vector<set<int>>(max_time);
  schedule_Ph_to_Sh = vector<set<int>>(max_time);
  schedule_infective_acute = vector<vector<pair<int, int>>>(max_time);
  schedule_infective_chronic = vector<vector<pair<int, int>>>(max_time);
  schedule_infective_recovery = vector<vector<pair<int, int>>>(max_time);
  
  // counts of host types
  H = H_init;
  Sh = H;
  Eh = vector<int>(n_demes);
  Ah = vector<int>(n_demes);
  Ch = vector<int>(n_demes);
  Ph = vector<int>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single population we can
  // simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(sum(H));
  next_host_ID = 0;
  next_inoc_ID = 0;
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only. Infections are seeded in the latent
  // stage, therefore there are no infective hosts initially.
  host_index = vector<vector<int>>(n_demes);
  host_infective_index = vector<vector<int>>(n_demes);
  int tmp1 = 0;
  for (int k = 0; k < n_demes; ++k) {
    host_index[k] = seq_int(tmp1, tmp1+H[k]-1);
    tmp1 += H[k];
  }
  
  // initialise the host population
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      int this_host = host_index[k][i];
      host_pop[this_host].init(this_host, next_host_ID, k,
                               Sh, Eh, Ah, Ch,
                               host_infective_index,
                               schedule_death,
                               schedule_Eh_to_Ah, schedule_Eh_to_Ch, schedule_Ah_to_Ch,
                               schedule_Ah_to_Sh, schedule_Ch_to_Sh,
                               schedule_Ah_to_Ph, schedule_Ch_to_Ph, schedule_Ph_to_Sh,
                               schedule_infective_acute, schedule_infective_chronic, schedule_infective_recovery,
                               sampler_age_stable, sampler_age_death,
                               sampler_duration_acute, sampler_duration_chronic,
                               sampler_time_treatment_acute, sampler_time_treatment_chronic);
    }
  }
  
  // seed initial infections
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < seed_infections[k]; ++i) {
      host_pop[host_index[k][i]].denovo_infection(0, next_inoc_ID);
    }
  }
  
  // counts of mosquito types
  M_total = sum(M);
  Sv = M;
  Ev = vector<int>(n_demes);
  Iv = vector<int>(n_demes);
  
  // objects for tracking mosquitoes that die in lag phase (described below)
  Ev_death = vector<vector<int>>(n_demes, vector<int>(v));
  
  // populations of mosquitoes at various stages
  Ev_pop = vector<vector<vector<Mosquito>>>(n_demes, vector<vector<Mosquito>>(v));
  Iv_pop = vector<vector<Mosquito>>(n_demes);
  
  // objects for storing daily values:
  // 0 = Sh, 1 = Eh, 2 = Ah, 3 = Ch, 4 = Sv, 5 = Ev, 6 = Iv, 7 = EIR
  daily_values = vector<vector<vector<double>>>(n_demes, vector<vector<double>>(max_time, vector<double>(8)));
  
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
  
  // initialise ring buffer that loops back to 0 when it exceeds (v-1). Used for
  // keeping track of mosquito deaths in the latent phase.
  int v_ringbuffer = 0;
  
  // loop through daily time steps
  for (int t = 0; t < max_time; ++t) {
    
    // update ring buffer index
    v_ringbuffer = (v_ringbuffer == v-1) ? 0 : v_ringbuffer + 1;
    
    // skip over first iteration to ensure user-defined values appear first as
    // first result
    if (t != 0) {
      
      //-------- SCHEDULED HUMAN EVENTS --------
      
      // scheduled deaths
      for (auto it = schedule_death[t].begin(); it != schedule_death[t].end(); ++it) {
        int this_host = *it;
        host_pop[this_host].death(next_host_ID, t);
      }
      
      // scheduled Eh to Ah
      for (auto it = schedule_Eh_to_Ah[t].begin(); it != schedule_Eh_to_Ah[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Eh_to_Ah(this_slot);
      }
      
      // scheduled Eh to Ch
      for (auto it = schedule_Eh_to_Ch[t].begin(); it != schedule_Eh_to_Ch[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Eh_to_Ch(this_slot);
      }
      
      // scheduled Ah to Ch
      for (auto it = schedule_Ah_to_Ch[t].begin(); it != schedule_Ah_to_Ch[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Ah_to_Ch(this_slot);
      }
      
      // scheduled Ah to Sh
      for (auto it = schedule_Ah_to_Sh[t].begin(); it != schedule_Ah_to_Sh[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Ah_to_Sh(this_slot);
      }
      
      // scheduled Ch to Sh
      for (auto it = schedule_Ch_to_Sh[t].begin(); it != schedule_Ch_to_Sh[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Ch_to_Sh(this_slot);
      }
      
      // scheduled become acutely infective
      for (auto it = schedule_infective_acute[t].begin(); it != schedule_infective_acute[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].begin_infective_acute(this_slot, t);
      }
      
      // scheduled become chronically infective
      for (auto it = schedule_infective_chronic[t].begin(); it != schedule_infective_chronic[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].begin_infective_chronic(this_slot, t);
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
        
        // one method of drawing infections in humans would be to draw the total
        // number of infectious bites from Binomial(H[k], prob_infectious_bite),
        // then loop through all of these and see which infections take hold by
        // drawing from Bernoulli with probability given by the host-specific
        // prob_infection. However, this is wasteful as a large number of
        // infectious bites are rejected. On the other hand, if the
        // prob_infection was constant then we could draw from Binomial(H[k],
        // prob_infection*prob_infectious_bite), after which every bite would
        // lead to infection, however, we cannot do this as prob_infection is
        // host-specific and changes over inoculations. Therefore, as a
        // middleground, draw from Binomial(H[k],
        // max_prob_infection*prob_infectious_bite), where max_prob_infection is
        // the largest value that prob_infection could possibly take. Loop
        // through these query infectious bites and draw from a Bernoulli
        // distribution relative to this value.
        //
        // For example, if prob_infection = {0.1, 0.05} then filter based on the
        // value 0.1, i.e. draw the number of query infected hosts from
        // Binomial(H[k], 0.1). Then loop these query hosts and draw from the
        // relative probability of infection, which is Bernoulli with
        // probability {0.1, 0.05}/0.1 = {1.0, 0.5}.
        
        int host_query_infection = rbinom1(H[k], max_prob_infection*prob_infectious_bite);
        for (int i = 0; i < host_query_infection; ++i) {
          
          // choose host at random
          int rnd1 = sample2(0, H[k]-1);
          int this_host = host_index[k][rnd1];
          
          // determine whether infectious bite is successful
          if (rbernoulli1(host_pop[this_host].get_prob_infection()/max_prob_infection)) {
            
            // infect host
            host_pop[this_host].infection(t, next_inoc_ID);
          }
          
        }  // end loop over query infectious bites
        
        
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
        
        // use the same method of drawing query infections as used when
        // infecting human hosts (see above)
        int mosq_query_infection = rbinom1(n_bite_infective, max_infectivity);
        for (int i = 0; i < mosq_query_infection; ++i) {
          
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
              int this_host_ID = host_pop[this_host].host_ID;
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
      daily_values[k][t] = {double(Sh[k]), double(Eh[k]), double(Ah[k]), double(Ch[k]),
                              double(Sv[k]), double(Ev[k]), double(Iv[k]),
                              EIR[k]};
    }
    
  }  // end loop through daily time steps
  
}

