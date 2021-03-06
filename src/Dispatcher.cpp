
#include "Dispatcher.h"
#include "probability_v11.h"

#include <fstream>

using namespace std;

//------------------------------------------------
// initialise
void Dispatcher::init(Parameters &params_) {
  
  // store pointer to parameters
  params = &params_;
  
  // open filestream to write transmission record to file
  if (params->save_transmission_record) {
    
    // open filestream
    if (!params->silent) {
      print("Opening filestream to transmission record");
    }
    transmission_record.open(params->transmission_record_location);
    
    // check that open
    if (!transmission_record.is_open()) {
      Rcpp::stop("unable to create transmission record at specified location. Check the path exists, and that you have write access");
    }
  }
  
  // initialise unique IDs for each inoculation
  next_inoc_ID = 1;
  
  // objects for sampling from probability distributions
  int sampler_draws = 1000;
  sampler_age_stable = Sampler(params->age_stable, sampler_draws);
  sampler_age_death = Sampler(params->age_death, sampler_draws);
  
  sampler_duration_acute = vector<Sampler>(params->n_duration_acute);
  for (int i = 0; i < params->n_duration_acute; ++i) {
    sampler_duration_acute[i] = Sampler(params->duration_acute[i], sampler_draws);
  }
  sampler_duration_chronic = vector<Sampler>(params->n_duration_chronic);
  for (int i = 0; i < params->n_duration_chronic; ++i) {
    sampler_duration_chronic[i] = Sampler(params->duration_chronic[i], sampler_draws);
  }
  sampler_time_treatment_acute = vector<Sampler>(params->n_time_treatment_acute);
  for (int i = 0; i < params->n_time_treatment_acute; ++i) {
    sampler_time_treatment_acute[i] = Sampler(params->time_treatment_acute[i], sampler_draws);
  }
  sampler_time_treatment_chronic = vector<Sampler>(params->n_time_treatment_chronic);
  for (int i = 0; i < params->n_time_treatment_chronic; ++i) {
    sampler_time_treatment_chronic[i] = Sampler(params->time_treatment_chronic[i], sampler_draws);
  }
  sampler_duration_prophylactic = vector<Sampler>(params->n_duration_prophylactic);
  for (int i = 0; i < params->n_duration_prophylactic; ++i) {
    sampler_duration_prophylactic[i] = Sampler(params->duration_prophylactic[i], sampler_draws);
  }
  
  // counts of host types
  H = params->H_init;
  Sh = H;
  Eh = vector<int>(params->n_demes);
  Ah = vector<int>(params->n_demes);
  Ch = vector<int>(params->n_demes);
  Ph = vector<int>(params->n_demes);
  
  // further counts of host types
  Ah_detectable_microscopy = vector<double>(params->n_demes);
  Ch_detectable_microscopy = vector<double>(params->n_demes);
  Ah_detectable_PCR = vector<double>(params->n_demes);
  Ch_detectable_PCR = vector<double>(params->n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single population we can
  // simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(sum(H));
  next_host_ID = 0;
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only. Infections are seeded in the latent
  // stage, therefore there are no infective hosts initially.
  host_index = vector<vector<int>>(params->n_demes);
  host_infective_index = vector<vector<int>>(params->n_demes);
  int tmp1 = 0;
  for (int k = 0; k < params->n_demes; ++k) {
    host_index[k] = seq_int(tmp1, tmp1 + H[k] - 1);
    tmp1 += H[k];
  }
  
  // initialise the host population
  for (int k = 0; k < params->n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      int this_index = host_index[k][i];
      host_pop[this_index].init(*params,
                                this_index, next_host_ID, k,
                                host_index,
                                host_infective_index,
                                sampler_age_stable, sampler_age_death,
                                sampler_duration_acute, sampler_duration_chronic,
                                sampler_time_treatment_acute, sampler_time_treatment_chronic,
                                sampler_duration_prophylactic);
    }
  }
  
  // counts of mosquito types
  M_total = sum(params->M);
  Sv = params->M;
  Ev = vector<int>(params->n_demes);
  Iv = vector<int>(params->n_demes);
  
  // objects for tracking mosquitoes that die in lag phase (process described below)
  Ev_death = vector<vector<int>>(params->n_demes, vector<int>(params->v));
  
  // populations of mosquitoes at various stages
  Ev_pop = vector<vector<vector<Mosquito>>>(params->n_demes, vector<vector<Mosquito>>(params->v));
  Iv_pop = vector<vector<Mosquito>>(params->n_demes);
  
  // objects for storing results
  daily_numer = vector<vector<double>>(params->max_time, vector<double>(params->n_daily_outputs));
  daily_denom = vector<vector<double>>(params->max_time, vector<double>(params->n_daily_outputs));
  sweep_numer = vector<double>(params->n_sweep_outputs);
  sweep_denom = vector<double>(params->n_sweep_outputs);
  
  // misc
  daily_EIR = vector<double>(params->n_demes);
  prob_infectious_bite = vector<double>(params->n_demes);
  
}

//------------------------------------------------
// run main simulation
void Dispatcher::run_simulation(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // start message
  if (!params->silent) {
    print("Running simulation");
  }
  
  // extract R utility functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // initialise ring buffer that loops back to 0 when it exceeds (v-1). Used for
  // keeping track of mosquito deaths in the latent phase.
  int v_ringbuffer = 0;
  
  // initialise vectors for storing daily output of each day
  vector<double> daily_numer_today(params->n_daily_outputs);
  vector<double> daily_denom_today(params->n_daily_outputs);
  
  // initialise indices that increment throughout simulation
  int sweep_index = 0;
  
  // loop through daily time steps
  for (int t = 0; t < params->max_time; ++t) {
    
    
    // update progress bar
    if (!params->silent) {
      int remainder = t % int(ceil(double(params->max_time)/100));
      if ((remainder == 0 && !params->pb_markdown) || ((t+1) == params->max_time)) {
        update_progress(args_progress, "pb_sim", t+1, params->max_time, true);
        if ((t+1) == params->max_time) {
          print("");
        }
      }
    }
    
    // update ring buffer index
    if (v_ringbuffer == (params->v - 1)) {
      v_ringbuffer = 0;
    } else {
      v_ringbuffer++;
    }
    
    // seed infections in first generation
    if (t == 0) {
      for (int k = 0; k < params->n_demes; ++k) {
        for (int i = 0; i < params->seed_infections[k]; ++i) {
          int this_host = host_index[k][i];
          host_pop[this_host].denovo_infection(t, next_inoc_ID, transmission_record);
        }
      }
    }
    
    
    //-------- MIGRATION --------
    
    
    // loop through all hosts, draw migration
    for (int i = 0; i < sum(H); ++i) {
      int this_deme = host_pop[i].deme;
      int new_deme = sample1(params->mig_mat[this_deme], 1.0);
      
      host_pop[i].migrate(new_deme);
    }
    
    
    //-------- MAIN LOOP THROUGH DEMES --------
    
    for (int k = 0; k < params->n_demes; ++k) {
      
      // recalculate number of hosts in this deme
      H[k] = host_index[k].size();
      if (H[k] == 0) {
        Rcpp::stop("empty deme");
      }
      
      
      //-------- NEW HUMAN EVENTS --------
      
      // get number of new infectious bites on humans
      daily_EIR[k] = params->a * Iv[k] / double(H[k]);
      
      // probability that each host receives an infectious bite
      prob_infectious_bite[k] = 1 - exp(-daily_EIR[k]);
      
      // one method of drawing infections in humans would be to draw the total
      // number of infectious bites from a Binomial(H[k], prob_infectious_bite)
      // distribution, then loop through all of these and see which infections
      // take hold by drawing from a Bernoulli with the probability given by the
      // host-specific prob_infection. However, this is wasteful as a large
      // number of infectious bites are rejected. On the other hand, if the
      // prob_infection was constant then we could draw from Binomial(H[k],
      // prob_infection*prob_infectious_bite), after which every bite would lead
      // to infection, however, we cannot do this as prob_infection is
      // host-specific and changes over inoculations. Therefore, as a
      // middleground, draw from Binomial(H[k],
      // max_prob_infection*prob_infectious_bite), where max_prob_infection is
      // the largest value that prob_infection could possibly take. Then loop
      // through these query infectious bites and draw from a Bernoulli
      // distribution relative to this value.
      //
      // For example, if prob_infection = {0.1, 0.05} then filter based on the
      // value 0.1, i.e. draw the number of query infected hosts from
      // Binomial(H[k], 0.1*prob_infectious_bite). Then loop through these query
      // hosts and draw from the relative probability of infection, which is
      // Bernoulli with probability {0.1, 0.05}/0.1 = {1.0, 0.5}.
      
      int host_query_infection = rbinom1(H[k], params->max_prob_infection * prob_infectious_bite[k]);
      for (int i = 0; i < host_query_infection; ++i) {
        
        // choose host at random
        int rnd1 = sample2(0, H[k] - 1);
        int this_host = host_index[k][rnd1];
        
        // determine whether infectious bite is successful
        if (rbernoulli1(host_pop[this_host].get_prob_infection() / params->max_prob_infection)) {
          
          // choose mosquito at random
          int rnd1 = sample2(0, Iv[k] - 1);
          
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
      double rate_bite_infective = params->a * host_infective_index[k].size() / double(H[k]); 
      
      // draw number of mosquitoes that bite infective host or die (competing
      // hazards)
      double prob_bite_infective_or_death = 1 - exp(-(rate_bite_infective + params->mu));
      int n_bite_infective_or_death = rbinom1(Sv[k], prob_bite_infective_or_death);
      
      // draw number of mosquitoes that bite infective host, rather than dying.
      double relative_prob_bite_infective = rate_bite_infective / (rate_bite_infective + params->mu);
      int n_bite_infective = rbinom1(n_bite_infective_or_death, relative_prob_bite_infective);
      
      // use the same method of drawing query infections as used when infecting
      // human hosts (see above), this time looping through mosquito infections
      int mosq_query_infection = rbinom1(n_bite_infective, params->max_infectivity);
      for (int i = 0; i < mosq_query_infection; ++i) {
        
        // choose host at random from infectives
        int rnd1 = sample2(0, host_infective_index[k].size() - 1);
        int this_host = host_infective_index[k][rnd1];
        
        // get infectivity and draw whether infection takes hold in mosquito
        double host_infectivity = host_pop[this_host].get_infectivity(t);
        if (rbernoulli1(host_infectivity / params->max_infectivity)) {
          
          // update deme counts
          Sv[k]--;
          Ev[k]++;
          
          // the majority of new mosquito infections will die in lag phase.
          // Schedule these deaths to move back into Sv in future steps.
          // Otherwise add to Ev_pop
          int mosq_time_death = rgeom1(params->prob_mosq_death) + 1;
          if (mosq_time_death <= params->v) {
            
            // schedule death for future time
            Ev_death[k][(v_ringbuffer + mosq_time_death) % params->v]++;
            
          } else {
            
            // sample inoc IDs from host
            vector<int> inoc_ID_vec = host_pop[this_host].get_inoc_ID_vec();
            
            // add to Ev_pop, scheduled to enter Iv_pop at future time
            Ev_pop[k][v_ringbuffer].emplace_back(inoc_ID_vec);
            
          }
          
        }
      } // end loop through query infective bites
      
      // deaths in Iv
      int death_Iv = rbinom1(Iv[k], params->prob_mosq_death);
      Sv[k] += death_Iv;
      Iv[k] -= death_Iv;
      for (int i = 0; i < death_Iv; ++i) {
        int rnd1 = sample2(0, Iv[k] - 1);
        quick_erase(Iv_pop[k], rnd1);
      }
      
    }  // end loop over demes
    
    
    //-------- SCHEDULED HUMAN EVENTS --------
    
    // loop through all hosts
    for (unsigned int i = 0; i < host_pop.size(); ++i) {
      
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
    
    
    //-------- STORE RESULTS --------
    
    // daily output
    if (params->any_daily_outputs) {
      
      // clear daily values vector
      fill(daily_numer_today.begin(), daily_numer_today.end(), 0.0);
      fill(daily_denom_today.begin(), daily_denom_today.end(), 0.0);
      
      // loop through hosts in each deme
      for (int k = 0; k < params->n_demes; ++k) {
        if (!params->daily_flag_deme[k]) {
          continue;
        }
        for (int i = 0; i < host_index[k].size(); ++i) {
          int this_host = host_index[k][i];
          int this_age = host_pop[this_host].get_age(t);
          
          // look up whether daily output required
          if (params->daily_map.count({k, this_age})) {
            vector<int> this_out = params->daily_map[{k, this_age}];
            
            // add to daily values
            for (unsigned int j = 0; j < this_out.size(); ++j) {
              host_pop[this_host].update_output(params->daily_measure[this_out[j]],
                                                params->daily_state[this_out[j]],
                                                params->daily_diagnostic[this_out[j]],
                                                params->daily_age_min[this_out[j]],
                                                params->daily_age_max[this_out[j]],
                                                params->daily_inoculations[this_out[j]],
                                                t,
                                                daily_numer_today[this_out[j]],
                                                daily_denom_today[this_out[j]]);
            }
            
          }
          
        }
      }
      
      
      // calculate EIR
      for (int i = 0; i < params->n_daily_outputs; i++) {
        if (params->daily_measure[i] == Measure_EIR) {
          int this_deme = params->daily_deme[i];
          if (this_deme == -1) {
            daily_numer_today[i] = 365 * mean(daily_EIR);
          } else {
            daily_numer_today[i] = 365 * daily_EIR[this_deme];
          }
        }
      }
      
      // store values
      daily_numer[t] = daily_numer_today;
      daily_denom[t] = daily_denom_today;
      
    }  // end daily outputs
    
    // sweep outputs
    if (params->any_sweep_outputs) {
      if (params->sweep_time_ordered[sweep_index] == t+1) {
        sweep_index++;
        
        // check all sweep outputs
        for (int i = 0; i < params->n_sweep_outputs; ++i) {
          if (params->sweep_time[i] == t+1) {
            
            double sweep_numer_today = 0.0;
            double sweep_denom_today = 0.0;
            
            int this_deme = params->sweep_deme[i];
            if (this_deme == -1) {
              for (unsigned int j = 0; j < host_pop.size(); ++j) {
                int this_age = host_pop[j].get_age(t);
                if ((this_age >= params->sweep_age_min[i]) && (this_age <= params->sweep_age_max[i])) {
                  host_pop[j].update_output(params->sweep_measure[i],
                                            params->sweep_state[i],
                                            params->sweep_diagnostic[i],
                                            params->sweep_age_min[i],
                                            params->sweep_age_max[i],
                                            params->sweep_inoculations[i],
                                            t, 
                                            sweep_numer_today, 
                                            sweep_denom_today);
                }
              }
            } else {
              for (unsigned int j = 0; j < host_index[this_deme].size(); ++j) {
                int this_host = host_index[this_deme][j];
                int this_age = host_pop[this_host].get_age(t);
                if ((this_age >= params->sweep_age_min[i]) && (this_age <= params->sweep_age_max[i])) {
                  host_pop[this_host].update_output(params->sweep_measure[i],
                                                    params->sweep_state[i],
                                                    params->sweep_diagnostic[i],
                                                    params->sweep_age_min[i],
                                                    params->sweep_age_max[i],
                                                    params->sweep_inoculations[i],
                                                    t, 
                                                    sweep_numer_today, 
                                                    sweep_denom_today);
                }
              }
            }
            
            // store values
            sweep_numer[i] = sweep_numer_today;
            sweep_denom[i] = sweep_denom_today;
            
          }
        }
        
      }
    }  // end sweep outputs
    
    /*
    // sample inoc IDs
    if (params->obtain_samples) {
      while (params->ss_time[index_survey] == t+1) {
        
        // get sampling parameters
        int this_deme = params->ss_deme[index_survey];
        int this_n = params->ss_n[index_survey];
        Diagnostic this_diagnosis = params->ss_diagnosis[index_survey];
        
        // obtain samples
        get_sample_details(t, this_deme, this_n, this_diagnosis);
        
        // increment index
        index_survey++;
      }
    }
    */
    
    // line break in transmission record at end of this time step
    if (params->save_transmission_record) {
      transmission_record << "\n";
    }
    
  }  // end main loop through daily time steps
  
  // close filestream to transmission record
  if (params->save_transmission_record) {
    if (!params->silent) {
      print("Closing filestream to transmission record");
    }
    transmission_record.close();
  }
  
}

//------------------------------------------------
// draw sample from deme
void Dispatcher::get_sample_details(int t, int deme, int n, Diagnostic diag) {
  /*
  // n cannot exceed human population size at this point in time
  if (n > H[deme]) {
    n = H[deme];
  }
  
  // draw vector by sampling without replacement
  vector<int> samp = sample4(n, 0, H[deme] - 1);
  
  // loop through all samples
  for (int i = 0; i < n; ++i) {
    
    // get host ID of this sample
    int this_index = host_index[deme][samp[i]];
    int this_host_ID = host_pop[this_index].host_ID;
    
    // find if positive for malaria parasites
    bool test_positive =  false;
    Status_host this_host_status = host_pop[this_index].get_host_status();
    if (this_host_status == Host_Ah) {
      
      // get positive prob for acute microscopy vs PCR
      double prob_positive = 0.0;
      if (diag == microscopy) {
        prob_positive = host_pop[this_index].get_detectability_microscopy_acute(t);
      } else if (diag == PCR) {
        prob_positive = host_pop[this_index].get_detectability_PCR_acute(t);
      }
      test_positive = rbernoulli1(prob_positive);
      
    } else if (this_host_status == Host_Ch) {
      
      // get positive prob for chronic microscopy vs PCR
      double prob_positive = 0.0;
      if (diag == microscopy) {
        prob_positive = host_pop[this_index].get_detectability_microscopy_chronic(t);
      } else if (diag == PCR) {
        prob_positive = host_pop[this_index].get_detectability_PCR_chronic(t);
      }
      test_positive = rbernoulli1(prob_positive);
      
    }
    
    // save basic details
    vector<int> this_details = {t+1, deme+1, this_host_ID, test_positive};
    
    // if positive then push back inoc IDs
    if (test_positive) {
      for (int j = 0; j < params->max_inoculations; ++j) {
        Status_asexual this_asexual = host_pop[this_index].inoc_status_asexual[j];
        if (this_asexual == Acute_asexual || this_asexual == Chronic_asexual) {
          this_details.push_back(host_pop[this_index].inoc_ID_vec[j]);
        }
      }
    }
    
    // push to sample_details
    //sample_details.push_back(this_details);
    
  }  // end i loop
  */
}
