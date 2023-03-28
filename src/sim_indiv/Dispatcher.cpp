
#include "Dispatcher.h"

using namespace std;

// unique IDs for hosts, mosquitoes and inoculations. Global variables should
// normally be avoided, but in this case these IDs should only ever increment
// and nothing depends on their value so code is arguable clearer with these
// defined globally rather than passing around
int next_host_ID = 1;
int next_mosq_ID = 1;
int next_inoc_ID = 1;

//------------------------------------------------
// initialise
void Dispatcher::init(Parameters &params_) {
  
  // store pointer to parameters
  params = &params_;
  
  // open filestream to write transmission record to file
  open_trans_record();
  
  // initialise population of human hosts over all demes and seed initial
  // infections
  host_pop.init(*params, transmission_record);
  host_pop.seed_infections();
  
  // initialise population of mosquitoes per deme
  mosq_pop = vector<Mosquito_pop>(params->n_demes);
  for (int k = 0; k < params->n_demes; ++k) {
    mosq_pop[k].init(*params, params->M[k]);
  }
  
  // initialise sampling objects and objects for storing results
  if (params->any_daily_outputs) {
    daily.init(*params, params->daily_df);
    daily_output = vector<vector<double>>(params->max_time, vector<double>(daily.n_rows));
  }
  if (params->any_survey_outputs) {
    survey.init(*params, params->surveys_df);
  }
  
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
    cpp11::stop("unable to create transmission record at specified location. Check that the path exists and that you have write access");
  }
  
  // write header line of transmission record
  transmission_record << "time,event,human_ID,mosquito_ID,child_inoc_ID,parent_inoc_ID,deme\n";
  
}

//------------------------------------------------
// run main simulation
void Dispatcher::run_simulation(cpp11::list &args_progress) {
  
  // start message
  if (!params->silent) {
    print("Running simulation");
  }
  
  // copy update_progress R function
  cpp11::function update_progress = cpp11::package("SIMPLEGEN")["update_progress"];
  
  //-----------------------------
  // main simulation loop
  
  // loop through daily time steps
  for (int t = 0; t < params->max_time; ++t) {
    
    //-----------------------------
    // 1) UPDATE MOSQUITOES
    // - shift values in ringbuffer
    // - loop through mosquito population and check for death or emergence. Change Sv, Ev, Iv as needed.
    
    for (int k = 0; k < params->n_demes; ++k) {
      mosq_pop[k].apply_ringbuffer_death();
      mosq_pop[k].update_pop(t);
    }
    
    //-----------------------------
    // 2) UPDATE HOSTS
    // - apply migration
    // - loop through and update life events
    
    host_pop.update_hosts(t);
    
    //-----------------------------
    // 3) INFECT MOSQUITOES
    // - loop through hosts, get infectivity, infect a binomial number of susceptible mosquitoes
    // - draw whether dies in latent phase. If yes, increase ringbuffer. If no, activate in population and schedule transition from latent and death.
    
    for (int k = 0; k < params->n_demes; ++k) {
      mosq_pop[k].draw_new_infections(host_pop, t, k);
    }
    
    //-----------------------------
    // 4) INFECT HOSTS
    // - loop through mosquito pop, binomial draw how many hosts each mosquito bites
    // - draw host at random and implement infectious bite
    
    for (int k = 0; k < params->n_demes; ++k) {
      host_pop.draw_new_infections(mosq_pop[k], t, k);
    }
    
    //-----------------------------
    // 5) STORE OUTPUT
    // - get host states
    
    // store daily output
    if (params->any_daily_outputs) {
      daily.calculate(daily_output[t], host_pop, mosq_pop, t);
    }
    
    // store survey output
    if (params->any_survey_outputs) {
      survey.take_survey(host_pop);
    }
    
    
    // update progress bar
    //update_progress(args_progress, "pb_sim", t, params->max_time, true);
    
  } // end t loop
  
  // store output
  cpp11::writable::doubles_matrix<> daily_output_final = copy_mat(daily_output);
  
  cpp11::writable::list survey_inoc_ID;
  for (size_t i = 0; i < survey.inoc_ID.size(); ++i) {
    survey_inoc_ID.push_back({cpp11::literals::operator""_nm("ID", 0) = survey.inoc_ID[i]});
  }
  
  tmp_output.push_back({cpp11::literals::operator""_nm("daily", 0) = daily_output_final});
  tmp_output.push_back({cpp11::literals::operator""_nm("inoc_ID", 0) = survey_inoc_ID});
  tmp_output.push_back({cpp11::literals::operator""_nm("host_ID", 0) = survey.host_ID});
}
