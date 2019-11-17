
#include "Dispatcher.h"
#include "probability_v9.h"

#include <fstream>

using namespace std;

//------------------------------------------------
// initialise
void Dispatcher::init() {
  
  // open filestream to write transmission record to file
  if (save_transmission_record) {
    
    // open filestream
    if (!silent) {
      print("Opening filestream to transmission record");
    }
    transmission_record.open(transmission_record_location);
    
    // check that open
    if (!transmission_record.is_open()) {
      Rcpp::stop("unable to create transmission record at specified location. Check the path exists, and that you have write access");
    }
  }
  
  // initialise unique IDs for each inoculation
  next_inoc_ID = 0;
  
  // objects for sampling from probability distributions
  int sampler_draws = 1000;
  sampler_age_stable = Sampler(age_stable, sampler_draws);
  sampler_age_death = Sampler(age_death, sampler_draws);
  
  sampler_duration_acute = vector<Sampler>(n_duration_acute);
  for (int i = 0; i < n_duration_acute; ++i) {
    sampler_duration_acute[i] = Sampler(duration_acute[i], sampler_draws);
  }
  sampler_duration_chronic = vector<Sampler>(n_duration_chronic);
  for (int i = 0; i < n_duration_chronic; ++i) {
    sampler_duration_chronic[i] = Sampler(duration_chronic[i], sampler_draws);
  }
  sampler_time_treatment_acute = vector<Sampler>(n_time_treatment_acute);
  for (int i = 0; i < n_time_treatment_acute; ++i) {
    sampler_time_treatment_acute[i] = Sampler(time_treatment_acute[i], sampler_draws);
  }
  sampler_time_treatment_chronic = vector<Sampler>(n_time_treatment_chronic);
  for (int i = 0; i < n_time_treatment_chronic; ++i) {
    sampler_time_treatment_chronic[i] = Sampler(time_treatment_chronic[i], sampler_draws);
  }
  sampler_duration_prophylactic = Sampler(duration_prophylactic, sampler_draws);
  
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
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only. Infections are seeded in the latent
  // stage, therefore there are no infective hosts initially.
  host_index = vector<vector<int>>(n_demes);
  host_infective_index = vector<vector<int>>(n_demes);
  int tmp1 = 0;
  for (int k = 0; k < n_demes; ++k) {
    host_index[k] = seq_int(tmp1, tmp1 + H[k] - 1);
    tmp1 += H[k];
  }
  
  // initialise the host population
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      int this_host = host_index[k][i];
      host_pop[this_host].init(this_host, next_host_ID, k,
                               host_infective_index,
                               sampler_age_stable, sampler_age_death,
                               sampler_duration_acute, sampler_duration_chronic,
                               sampler_time_treatment_acute, sampler_time_treatment_chronic,
                               sampler_duration_prophylactic);
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
  
  // objects for storing results
  // daily values: 0 = Sh, 1 = Eh, 2 = Ah, 3 = Ch, 4 = Ph, 5 = Sv, 6 = Ev, 7 = Iv, 8 = EIR
  daily_values = vector<vector<vector<double>>>(n_demes, vector<vector<double>>(max_time, vector<double>(9)));
  
  // misc
  EIR = vector<double>(n_demes);
  
}

//------------------------------------------------
// run main simulation
void Dispatcher::run_simulation(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // start message
  if (!silent) {
    print("Running simulation");
  }
  
  // extract R utility functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // initialise ring buffer that loops back to 0 when it exceeds (v-1). Used for
  // keeping track of mosquito deaths in the latent phase.
  int v_ringbuffer = 0;
  
  // loop through daily time steps
  int index_obtain_samples = 0;
  for (int t = 0; t < max_time; ++t) {
    
    // update progress bar
    if (!silent) {
      update_progress(args_progress, "pb_sim", t+1, max_time, true);
    }
    
    // update ring buffer index
    v_ringbuffer = (v_ringbuffer == v-1) ? 0 : v_ringbuffer + 1;
    
    // seed infections in first generation
    if (t == 0) {
      for (int k = 0; k < n_demes; ++k) {
        for (int i = 0; i < seed_infections[k]; ++i) {
          host_pop[host_index[k][i]].denovo_infection(t, next_inoc_ID, transmission_record);
        }
      }
    }
    
    
    //-------- STORE RESULTS --------
    // NB. results are stored at this early stage so that user-defined initial
    // conditions (e.g. the number of seeding infections) are stored as the
    // first result.
    
    // update counts of each host status in each deme
    update_host_counts();
    
    // store daily values
    for (int k = 0; k < n_demes; ++k) {
      daily_values[k][t] = {double(Sh[k]), double(Eh[k]), double(Ah[k]), double(Ch[k]), double(Ph[k]),
                            double(Sv[k]), double(Ev[k]), double(Iv[k]),
                            EIR[k]};
    }
    
    // sample inoc IDs
    if (obtain_samples) {
      while (ss_time[index_obtain_samples] == t) {
        
        // get sampling parameters
        int this_deme = ss_deme[index_obtain_samples];
        int this_n = ss_n[index_obtain_samples];
        
        // obtain samples
        for (int i = 0; i < this_n; ++i) {
          get_sample_details(t, this_deme);
        }
        
        // increment index
        index_obtain_samples++;
      }
    }
    
    
    
    //-------- MIGRATION --------
    // TODO - loop through hosts - check for migration
    
    
    // loop through demes
    for (int k = 0; k < n_demes; ++k) {
      
      
      //-------- NEW HUMAN EVENTS --------
      
      // get number of new infectious bites on humans
      EIR[k] = a*Iv[k]/double(H[k]);
      double prob_infectious_bite = 1 - exp(-EIR[k]);  // probability of new infectious bite on host
      
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
          
          // choose mosquito at random
          int rnd1 = sample2(0, Iv[k]-1);
          
          // infect host
          host_pop[this_host].infection(t, next_inoc_ID, Iv_pop[k][rnd1], transmission_record);
          
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
            
            // sample inoc IDs from host
            vector<int> inoc_ID_vec = host_pop[this_host].get_inoc_ID_vec();
            
            // add to Ev_pop, scheduled to enter Iv_pop at future time
            Ev_pop[k][v_ringbuffer].emplace_back(inoc_ID_vec);
            
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
      
      
      //-------- SCHEDULED HUMAN EVENTS --------
      
      // loop through all hosts
      for (int i = 0; i < int(host_pop.size()); ++i) {
        
        // check for host death
        if (host_pop[i].death_day == t) {
          host_pop[i].death(next_host_ID, t);
        }
        
        // check for host change of prophylactic state
        if (host_pop[i].t_prophylaxis_stop == t) {
          host_pop[i].end_prophylaxis();
        }
        
        // apply any scheduled inoc-level events
        host_pop[i].check_inoc_event(t);
        
      }
      
    }  // end loop over demes
    
    // line break at end of this time step
    transmission_record << "\n";
    
  }  // end loop through daily time steps
  
  // close filestream to transmission record
  if (save_transmission_record) {
    if (!silent) {
      print("Closing filestream to transmission record");
    }
    transmission_record.close();
  }
  
}

//------------------------------------------------
// update host counts
void Dispatcher::update_host_counts() {
  
  // reset counts
  fill(Sh.begin(), Sh.end(), 0);
  fill(Eh.begin(), Eh.end(), 0);
  fill(Ah.begin(), Ah.end(), 0);
  fill(Ch.begin(), Ch.end(), 0);
  fill(Ph.begin(), Ph.end(), 0);
  
  // loop through all hosts, update counts in given deme
  for (int i = 0; i < int(host_pop.size()); ++i) {
    int this_deme = host_pop[i].deme;
    switch(host_pop[i].get_host_status()) {
    case Host_Sh:
      Sh[this_deme]++;
      break;
    case Host_Eh:
      Eh[this_deme]++;
      break;
    case Host_Ah:
      Ah[this_deme]++;
      break;
    case Host_Ch:
      Ch[this_deme]++;
      break;
    case Host_Ph:
      Ph[this_deme]++;
      break;
    default:
      Rcpp::stop("invalid host status in update_host_counts()");
    }
  }
  
}

//------------------------------------------------
// draw sample from deme
void Dispatcher::get_sample_details(int t, int deme) {
  
  // store target time and deme
  vector<int> this_details;
  this_details.push_back(t);
  this_details.push_back(deme);
  
  // choose a random host and push back host ID
  int rnd1 = sample2(0, H[deme] - 1);
  int this_index = host_index[deme][rnd1];
  int this_host_ID = host_pop[this_index].host_ID;
  this_details.push_back(this_host_ID);
  
  // find if positive for malaria parasites and push back test results
  bool test_positive = false;
  Status_host this_host_status = host_pop[this_index].get_host_status();
  if (this_host_status == Host_Ah || this_host_status == Host_Ch) {
    test_positive = true;
  }
  this_details.push_back(test_positive);
  
  // if positive, push back inoc IDs
  if (test_positive) {
    for (int i = 0; i < max_inoculations; ++i) {
      Status_asexual this_asexual = host_pop[this_index].inoc_status_asexual[i];
      if (this_asexual == Acute_asexual || this_asexual == Chronic_asexual) {
        this_details.push_back(host_pop[this_index].inoc_ID_vec[i]);
      }
    }
  }
  
  // push to sample_details
  sample_details.push_back(this_details);
}

