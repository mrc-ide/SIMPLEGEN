
#include "Host_pop.h"
#include "Mosquito_pop.h"

using namespace std;

//------------------------------------------------
// initialise
void Host_pop::init(Parameters &params_, ofstream &transmission_record_) {
  
  // store pointers
  params = &params_;
  transmission_record = &transmission_record_;
  
  // get total number of hosts in each deme
  H = params->H_init;
  
  // initialise single vector of human hosts over all demes. This is preferable
  // to using separate vectors of hosts for each deme, as this would mean moving
  // hosts around due to migration. With a single population we can simply
  // change the "deme" attribute of a host to represent migration (and change
  // the host_index, see objects below)
  host_vec = vector<Host>(sum(H));
  
  // for each deme, store the integer index of all hosts in that deme. The
  // actual host object can be found by refering to this index within host_pop
  host_index = vector<vector<int>>(params->n_demes);
  int deme_start = 0;
  for (int k = 0; k < params->n_demes; ++k) {
    host_index[k] = seq_int(deme_start, deme_start + H[k] - 1);
    deme_start += H[k];
  }
  
  // TODO - more efficient to store a large vector and manage memory as in mosq_pop?
  
  // initialise all hosts
  int i2 = 0;
  for (int k = 0; k < params->n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      host_vec[i2].init(*params, *transmission_record, k);
      i2++;
    }
  }
  
}

//------------------------------------------------
// seed infections
void Host_pop::seed_infections() {
  
  // sample hosts without replacement within each deme
  for (int k = 0; k < params->n_demes; ++k) {
    vector<int> samp_replace = sample4(params->seed_infections[k], 0, H[k] - 1);
    for (int i = 0; i < params->seed_infections[k]; ++i) {
      int this_host = host_index[k][samp_replace[i]];
      host_vec[this_host].infect(0, -1, -1);
    }
  }
  
}

//------------------------------------------------
// update hosts over all demes given life events
void Host_pop::update_hosts(int t) {
  
  // check for death in all hosts. This is done first, before updating hosts,
  // because death can result in a host moving deme back to their home deme
  // meaning some hosts could miss their update if done altogether
  for (int k = 0; k < params->n_demes; ++k) {
    for (int i = 0; i < H[k]; ++i) {
      int this_host = host_index[k][i];
      if (host_vec[this_host].time_death == t) {
        
        // move back to home deme
        if (host_vec[this_host].deme != host_vec[this_host].home_deme) {
          quick_erase(host_index[k], i);
          int home_deme = host_vec[this_host].home_deme;
          host_index[home_deme].push_back(this_host);
          H[k]--;
          H[home_deme]++;
        }
        
        // apply death event
        host_vec[this_host].death(t);
      }
    }
  }
  
  // now update all hosts with other events
  for (size_t i = 0; i < host_vec.size(); ++i) {
    host_vec[i].update(t);
  }
}

//------------------------------------------------
// infect hosts from mosquitoes
void Host_pop::draw_new_infections(Mosquito_pop &mosq_pop, int t, int k) {
  
  // calculate probability of a single host being bitten by a single infectious
  // mosquito
  double prob_infectious_bite = params->a / double(H[k]);
  
  // loop through all mosquitoes int this deme. Focus on infectious mosquitoes
  for (int i = 0; i < mosq_pop.mosq_vec_n; ++i) {
    if (mosq_pop.mosq_vec[i].infectious_on) {
      
      // draw binomial number of hosts bitten by this mosquito
      int n_bites = rbinom1(H[k], prob_infectious_bite);
      
      // for each bite, sample hosts at random
      for (int j = 0; j < n_bites; ++j) {
        int samp1 = sample2(0, H[k] - 1);
        int this_host = host_index[k][samp1];
        
        // draw whether infection successful
        double prob_infection = host_vec[this_host].get_prob_infection(t);
        if (rbernoulli1(prob_infection)){
          
          // if successful then flush mosquito buffer to transmission record
          mosq_pop.mosq_vec[i].write_buffer(*transmission_record);
          
          // apply infection
          host_vec[this_host].infect(t, mosq_pop.mosq_vec[i].mosq_ID, mosq_pop.mosq_vec[i].inoc_ID);
        }
      }
    }
  }
  
  
}
