
#include "Dispatcher.h"
#include "probability_v14.h"

#include <fstream>

using namespace std;

//------------------------------------------------
// initialise
void Dispatcher::init(Parameters &params_) {
  
  // store pointer to parameters
  params = &params_;
  
  // open filestream to write transmission record to file
  open_trans_record();
  
  // initialise unique IDs for hosts, mosquitoes and inoculations
  next_host_ID = 1;
  next_mosq_ID = 1;
  next_infection_ID = 1;
  
  // create objects for sampling from probability distributions
  int sampler_draws = 1000;
  sampler_age_stable = Sampler(params->age_stable, sampler_draws);
  sampler_age_death = Sampler(params->age_death, sampler_draws);
  sampler_duration_acute = make_sampler_vec(params->duration_acute, sampler_draws);
  sampler_duration_chronic = make_sampler_vec(params->duration_chronic, sampler_draws);
  sampler_time_treatment_acute = make_sampler_vec(params->time_treatment_acute, sampler_draws);
  sampler_time_treatment_chronic = make_sampler_vec(params->time_treatment_chronic, sampler_draws);
  sampler_duration_prophylactic = make_sampler_vec(params->duration_prophylactic, sampler_draws);
  
  // counts of host types
  H = params->H_init;
  Sh = H;
  Eh = vector<int>(params->n_demes);
  Ah = vector<int>(params->n_demes);
  Ch = vector<int>(params->n_demes);
  Ph = vector<int>(params->n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single population we can
  // simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(sum(H));
  
  // for each deme, store the integer index of all hosts in that deme. The
  // actual host object can be found by refering to this index within host_pop.
  // Do the same for infective hosts only. Infections are seeded in the latent
  // stage, therefore there are no infective hosts initially.
  host_index = vector<vector<int>>(params->n_demes);
  host_infective_index = vector<vector<int>>(params->n_demes);
  int tmp_int1 = 0;
  for (int k = 0; k < params->n_demes; ++k) {
    host_index[k] = seq_int(tmp_int1, tmp_int1 + H[k] - 1);
    tmp_int1 += H[k];
  }
  
  // initialise the host population
  for (int k = 0; k < params->n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      int this_index = host_index[k][i];
      host_pop[this_index].init(*params,
                                this_index, next_host_ID, k,
                                host_index, host_infective_index,
                                sampler_age_stable, sampler_age_death,
                                sampler_duration_acute, sampler_duration_chronic,
                                sampler_time_treatment_acute, sampler_time_treatment_chronic,
                                sampler_duration_prophylactic);
    }
  }
  
  // counts of mosquito types
  Sv = params->M;
  Ev = vector<int>(params->n_demes);
  Iv = vector<int>(params->n_demes);
  
  // objects for tracking mosquitoes that die in lag phase (process described below)
  Ev_death = vector<vector<int>>(params->n_demes, vector<int>(params->v));
  
  // populations of mosquitoes at various stages
  Ev_pop = vector<vector<vector<Mosquito>>>(params->n_demes, vector<vector<Mosquito>>(params->v));
  Iv_pop = vector<vector<Mosquito>>(params->n_demes);
  
  // objects for storing results
  daily_output = vector<vector<double>>(params->max_time, vector<double>(params->n_daily_outputs));
  sweep_output = vector<double>(params->n_sweep_outputs);
  
  // misc
  daily_EIR = vector<double>(params->n_demes);
  prob_infectious_bite = vector<double>(params->n_demes);
  
}

//------------------------------------------------
// open transmission record
void Dispatcher::open_trans_record() {
  
  // return if not saving transmission record
  if (!params->save_transmission_record) {
    return;
  }
  
  // open filestream
  if (!params->silent) {
    print("Opening filestream to transmission record");
  }
  transmission_record.open(params->transmission_record_location);
  
  // check that open
  if (!transmission_record.is_open()) {
    Rcpp::stop("unable to create transmission record at specified location. Check the path exists, and that you have write access");
  }
  
  // write header line of transmission record
  transmission_record << "time,event,human_ID,mosquito_ID,child_infection_ID,parent_infection_ID,deme\n";
  
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
  vector<double> daily_output_today(params->n_daily_outputs);
  
  // initialise indices that increment throughout simulation
  int sweep_index = 0;
  int surveys_expanded_index = 0;
  
  // initialise temporary objects used in survey output
  vector<set<int>> study_sampled_IDs(params->n_survey_outputs);
  vector<int> surveys_sample_size_remaining(params->n_survey_outputs);
  
  // loop through daily time steps
  for (int t = 0; t < params->max_time; ++t) {
    
    // update progress bar
    if (!params->silent) {
      int remainder = t % (int)ceil((double)params->max_time / 100);
      if ((remainder == 0 && !params->pb_markdown) || ((t + 1) == params->max_time)) {
        update_progress(args_progress, "pb_sim", t + 1, params->max_time, true);
        if ((t + 1) == params->max_time) {
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
          host_pop[this_host].denovo_infection(t, next_infection_ID, transmission_record);
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
      daily_EIR[k] = params->a * Iv[k] / (double)H[k];
      
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
      // host-specific and changes over infections. Therefore, as a
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
          
          // this mosquito must have been infected at an earlier time, and so it
          // may have info stored in its buffer. Write this buffer to
          // transmission record now so that chronologically it appears the
          // mosquito was infected before passing on infection.
          Iv_pop[k][rnd1].write_buffer(transmission_record);
          
          if (Iv_pop[k][rnd1].mosquito_ID == 122) {
            bar();
            Iv_pop[k][rnd1].print_status();
          }
          
          // infect host, and write this event to transmission record
          host_pop[this_host].infection(t, next_infection_ID,
                                        Iv_pop[k][rnd1].mosquito_ID,
                                        Iv_pop[k][rnd1].infection_ID,
                                        transmission_record);
          
        }
        
      }  // end loop over query infectious bites
      
      
      //-------- SCHEDULED MOSQUITO EVENTS --------
      
      // deaths in Ev are reborn in Sv
      Sv[k] += Ev_death[k][v_ringbuffer];
      Ev[k] -= Ev_death[k][v_ringbuffer];
      Ev_death[k][v_ringbuffer] = 0;
      
      // move Ev into Iv
      int delta_Ev = (int)Ev_pop[k][v_ringbuffer].size();
      if (delta_Ev > 0) {
        Ev[k] -= delta_Ev;
        Iv[k] += delta_Ev;
        push_back_multiple(Iv_pop[k], Ev_pop[k][v_ringbuffer]);
        Ev_pop[k][v_ringbuffer].clear();
      }
      
      
      //-------- NEW MOSQUITO EVENTS --------
      
      // rate of mosquito biting infective host
      int n_infective = (int)host_infective_index[k].size();
      double rate_bite_infective = params->a * n_infective / (double)H[k]; 
      
      // draw number of mosquitoes that bite infective host or die (these are
      // competing hazards)
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
          // Otherwise add to Ev_pop.
          int mosq_time_death = rgeom1(params->prob_mosq_death) + 1;
          if (mosq_time_death <= params->v) {
            
            // schedule death for future time
            Ev_death[k][(v_ringbuffer + mosq_time_death) % params->v]++;
            
          } else {
            
            // create new mosquito object and infect from host
            Mosquito m;
            m.set_mosquito_ID(next_mosq_ID);
            m.infection(t, next_infection_ID, host_pop[this_host]);
            
            if (m.mosquito_ID == 122) {
              foo();
              m.print_status();
              host_pop[this_host].print_inoc_state();
            }
            
            // add to Ev_pop, scheduled to enter Iv_pop at future time
            Ev_pop[k][v_ringbuffer].emplace_back(m);
            
          }
          
        }
      } // end loop through query infective bites
      
      // deaths in Iv are reborn in Sv
      int death_Iv = rbinom1(Iv[k], params->prob_mosq_death);
      Sv[k] += death_Iv;
      Iv[k] -= death_Iv;
      for (int i = 0; i < death_Iv; ++i) {
        int rnd1 = sample2(0, Iv[k] - 1);
        quick_erase(Iv_pop[k], rnd1);
      }
      
    }  // end loop over demes
    
    
    //-------- SCHEDULED HUMAN EVENTS --------
    
    // TODO - remove
    for (unsigned int i = 0; i < host_pop.size(); ++i) {
      Host h = host_pop[i];
      for (int j = 0; j < params->max_inoculations; ++j) {
        if (!h.inoc_active[j] && (h.inoc_state_sexual[j] != Inactive_sexual)) {
          print(i, j);
          foobar();
          Rcpp::stop("");
        }
      }
    }
    
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
      
      // apply any scheduled infection-level events
      host_pop[i].check_infection_event(t);
    }
    
    
    //-------- STORE RESULTS --------
    
    // daily output
    if (params->any_daily_outputs) {
      
      // clear daily values vector
      fill(daily_output_today.begin(), daily_output_today.end(), 0.0);
      
      // loop through demes
      for (int k = 0; k < params->n_demes; ++k) {
        
        // skip this deme if it not required in any output
        if (!params->daily_flag_deme[k]) {
          continue;
        }
        
        // loop through hosts in this deme
        for (int i = 0; i < host_index[k].size(); ++i) {
          int this_host = host_index[k][i];
          int this_age = host_pop[this_host].get_age(t);
          
          // skip this individual if it not required in any output
          if (params->daily_map.count({k, this_age}) == 0) {
            continue;
          }
          
          // get which rows of daily output this individual applies to
          vector<int> which_rows = params->daily_map[{k, this_age}];
          
          // get daily output for all required rows
          for (unsigned int j = 0; j < which_rows.size(); ++j) {
            int w = which_rows[j];
            daily_output_today[w] += host_pop[this_host].get_output(params->daily_measure[w],
                                                                    params->daily_state[w],
                                                                    params->daily_diagnostic[w],
                                                                    params->daily_age_min[w],
                                                                    params->daily_age_max[w],
                                                                    t);
          }
          
        }
      }
      
      // calculate EIR
      for (int i = 0; i < params->n_daily_outputs; i++) {
        if (params->daily_measure[i] == Measure_EIR) {
          int this_deme = params->daily_deme[i];
          if (this_deme == -1) {
            daily_output_today[i] = 365 * mean(daily_EIR);
          } else {
            daily_output_today[i] = 365 * daily_EIR[this_deme];
          }
        }
      }
      
      // divide through by denominator and change units
      for (int i = 0; i < params->n_daily_outputs; i++) {
        if (params->daily_measure[i] == Measure_prevalence) {
          int this_deme = params->daily_deme[i];
          if (this_deme == -1) {
            daily_output_today[i] *= 100.0 / double(sum(H));
          } else {
            daily_output_today[i] *= 100.0 / double(H[this_deme]);
          }
        } else if (params->daily_measure[i] == Measure_incidence) {
          int this_deme = params->daily_deme[i];
          if (this_deme == -1) {
            daily_output_today[i] *= 365.0 / double(sum(H));
          } else {
            daily_output_today[i] *= 365.0 / double(H[this_deme]);
          }
        }
      }
      
      // store values
      daily_output[t] = daily_output_today;
      
    }  // end daily outputs
    
    // sweep outputs
    if (params->any_sweep_outputs) {
      if (params->sweep_time_ordered[sweep_index] == (t + 1)) {
        sweep_index++;
        
        // check all sweep outputs
        for (int i = 0; i < params->n_sweep_outputs; ++i) {
          if (params->sweep_time[i] == (t + 1)) {
            
            // object for storing todays results
            double sweep_output_today = 0.0;
            
            int this_deme = params->sweep_deme[i];
            if (this_deme == -1) {
              for (unsigned int j = 0; j < host_pop.size(); ++j) {
                int this_age = host_pop[j].get_age(t);
                if ((this_age >= params->sweep_age_min[i]) && (this_age <= params->sweep_age_max[i])) {
                  sweep_output_today += host_pop[j].get_output(params->sweep_measure[i],
                                                               params->sweep_state[i],
                                                               params->sweep_diagnostic[i],
                                                               params->sweep_age_min[i],
                                                               params->sweep_age_max[i],
                                                               t);
                }
              }
            } else {
              for (unsigned int j = 0; j < host_index[this_deme].size(); ++j) {
                int this_host = host_index[this_deme][j];
                int this_age = host_pop[this_host].get_age(t);
                if ((this_age >= params->sweep_age_min[i]) && (this_age <= params->sweep_age_max[i])) {
                  sweep_output_today += host_pop[this_host].get_output(params->sweep_measure[i],
                                                                       params->sweep_state[i],
                                                                       params->sweep_diagnostic[i],
                                                                       params->sweep_age_min[i],
                                                                       params->sweep_age_max[i],
                                                                       t);
                }
              }
            }
            
            // divide through by denominator and change units
            if (params->sweep_measure[i] == Measure_prevalence) {
              if (this_deme == -1) {
                sweep_output_today *= 100.0 / double(sum(H));
              } else {
                sweep_output_today *= 100.0 / double(H[this_deme]);
              }
            } else if (params->sweep_measure[i] == Measure_incidence) {
              if (this_deme == -1) {
                sweep_output_today *= 365.0 / double(sum(H));
              } else {
                sweep_output_today *= 365.0 / double(H[this_deme]);
              }
            }
            
            // store values
            sweep_output[i] = sweep_output_today;
            
          }
        }
        
      }
    }  // end sweep outputs
    
    // survey outputs
    if (params->any_survey_outputs) {
      
      // all required surveys for this day
      while (params->surveys_expanded_sampling_time[surveys_expanded_index] == (t + 1)) {
        
        // get reporting time
        int reporting_time = params->surveys_expanded_reporting_time[surveys_expanded_index];
        
        // get corresponding row in survey dataframe
        int survey_index = params->surveys_expanded_study_ID[surveys_expanded_index] - 1;
        
        // get survey parameters
        int surveys_t_start = params->surveys_t_start[survey_index];
        int surveys_t_end = params->surveys_t_end[survey_index];
        int surveys_deme = params->surveys_deme[survey_index];
        Measure surveys_measure = params->surveys_measure[survey_index];
        Sampling surveys_sampling = params->surveys_sampling[survey_index];
        Diagnostic surveys_diagnostic = params->surveys_diagnostic[survey_index];
        int surveys_age_min = params->surveys_age_min[survey_index];
        int surveys_age_max = params->surveys_age_max[survey_index];
        double surveys_sample_size = params->surveys_sample_size[survey_index];
        bool surveys_sample_proportion = (surveys_sample_size <= 1.0);
        int surveys_n_days_remaining = surveys_t_end - (t + 1) + 1;
        
        // if first day of survey then finalise sample sizes
        if (surveys_t_start == (t + 1)) {
          if (surveys_measure == Measure_prevalence) {
            if (surveys_sample_proportion) {
              if (surveys_deme == -1) {
                surveys_sample_size_remaining[survey_index] = round(surveys_sample_size * sum(H));
              } else {
                surveys_sample_size_remaining[survey_index] = round(surveys_sample_size * H[surveys_deme]);
              }
            } else {
              surveys_sample_size_remaining[survey_index] = (int)surveys_sample_size;
            }
          } else if (surveys_measure == Measure_incidence) {
            if (!surveys_sample_proportion) {
              surveys_sample_size_remaining[survey_index] = (int)surveys_sample_size;
            }
          }
        }
        
        // if prevalence study then spread sampling out evenly over days. If
        // incidence study then sampling continues on any given day until we
        // reach the total
        int max_samples_today = surveys_sample_size_remaining[survey_index];
        if (surveys_measure == Measure_prevalence) {
          max_samples_today = round(surveys_sample_size_remaining[survey_index] / double(surveys_n_days_remaining));
        }
        
        // hosts will be sampled by looping through a list of target host
        // indices and querying whether they meet study criteria. Get this list
        // differently depending on whether sampling one deme or all demes
        vector<int> query_indices;
        if (surveys_deme == -1) {
          query_indices = seq_int(0, sum(H) - 1);
        } else {
          query_indices = host_index[surveys_deme];
        }
        reshuffle(query_indices);
        
        // loop through query list until sampled enough hosts
        int n_sampled = 0;
        for (int i = 0; i < query_indices.size(); ++i) {
          
          // break if no remaining samples required
          if ((surveys_measure == Measure_prevalence) || 
              ((surveys_measure == Measure_incidence) && !surveys_sample_proportion)) {
            if (n_sampled == max_samples_today) {
              break;
            }
          }
          
          // skip if host already sampled in this study
          int this_host_ID = host_pop[query_indices[i]].host_ID;
          if (study_sampled_IDs[survey_index].count(this_host_ID)) {
            continue;
          }
          
          // establish whether host meets study criteria and simultaneously get
          // malaria positive status by various diagnostics
          bool true_positive = false;
          bool microscopy_positive = false;
          bool PCR_positive = false;
          bool include_host = host_pop[query_indices[i]].get_indiv_output(surveys_measure,
                                                                          surveys_sampling,
                                                                          surveys_diagnostic,
                                                                          surveys_age_min,
                                                                          surveys_age_max,
                                                                          t,
                                                                          true_positive,
                                                                          microscopy_positive,
                                                                          PCR_positive);
          
          // if sampling a proportion of incident cases then draw whether to
          // keep or discard this sample
          if ((surveys_measure == Measure_incidence) && surveys_sample_proportion) {
            if (!rbernoulli1(surveys_sample_size)) {
              continue;
            }
          }
          
          // if host to be included
          if (include_host) {
            Host this_host = host_pop[query_indices[i]];
            
            // save individual-level output to list
            Rcpp::List this_sample;
            this_sample["study_ID"] = survey_index + 1;
            this_sample["sampling_time"] = t + 1;
            this_sample["reporting_time"] = reporting_time;
            this_sample["sample_ID"] = this_host.host_ID;
            this_sample["home_deme"] = this_host.home_deme + 1;
            this_sample["sampled_deme"] = this_host.deme + 1;
            this_sample["true_positive"] = true_positive;
            this_sample["PCR_positive"] = PCR_positive;
            this_sample["microscopy_positive"] = microscopy_positive;
            this_sample["age"] = this_host.get_age(t);
            surveys_indlevel_output.push_back(this_sample);
            
            // store infection IDs in separate object
            vector<int> infection_IDs;
            for (int j = 0; j < this_host.inoc_ID_vec.size(); ++j) {
              if ((this_host.inoc_state_asexual[j] == Acute_asexual) || (this_host.inoc_state_asexual[j] == Chronic_asexual)) {
                infection_IDs.push_back(this_host.inoc_ID_vec[j]);
              }
            }
            surveys_indlevel_output_infection_IDs.push_back(infection_IDs);
            
            // add to list of sampled hosts
            study_sampled_IDs[survey_index].insert(this_host_ID);
            
            // update number of samples taken and remaining
            n_sampled++;
            surveys_sample_size_remaining[survey_index]--;
          }
          
        }  // end loop through query list
        
        // increment row of expanded surveys dataframe
        surveys_expanded_index++;
      }
    }  // end survey outputs
    
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
