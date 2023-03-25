
#include "Mosquito.h"

using namespace std;

// link to globals
extern int next_mosq_ID;
extern int next_inoc_ID;

//------------------------------------------------
// default constructor
Mosquito::Mosquito() {
  
  // set dummy values
  mosq_ID = -1;
  inoc_ID = -1;
  infectious_on = false;
  time_infectious = -1;
  time_death = -1;
  source_time = -1;
  source_deme = -1;
  source_host_ID = -1;
}

//------------------------------------------------
// initialise using values from human host. Overwrites existing values if
// present
void Mosquito::init_from_host(Parameters &params, int t, int t_death, Host &host) {
  
  // set IDs and increment
  mosq_ID = next_mosq_ID++;
  inoc_ID = next_inoc_ID++;
  
  // time of death and of transition from latent state
  infectious_on = false;
  time_infectious = t + params.v + 1;
  time_death = t_death;
  
  // record properties of host for writing to transmission record later
  source_time = t;
  source_deme = host.deme;
  source_host_ID = host.host_ID;
  
  // copy over all inoc IDs in sexual stage
  source_inoc_ID.clear();
  for (size_t i = 0; i < host.inoc.size(); ++i) {
    if ((host.inoc[i].state_sexual == Acute_sexual) || (host.inoc[i].state_sexual == Chronic_sexual)) {
      source_inoc_ID.push_back(host.inoc[i].ID);
    }
  }
}

//------------------------------------------------
// initialise by copying over values from an existing mosquito
void Mosquito::init_from_mosq(Mosquito &mosq)  {
  
  // copy over all values
  mosq_ID = mosq.mosq_ID;
  inoc_ID = mosq.inoc_ID;
  infectious_on = mosq.infectious_on;
  time_infectious = mosq.time_infectious;
  time_death = mosq.time_death;
  source_time = mosq.source_time;
  source_deme = mosq.source_deme;
  source_host_ID = mosq.source_host_ID;
  source_inoc_ID = mosq.source_inoc_ID;
}

//------------------------------------------------
// transition from latent to infetious state
void Mosquito::Ev_to_Iv() {
  infectious_on = true;
}

//------------------------------------------------
// flush buffer to transmission record
void Mosquito::write_buffer(std::ofstream &transmission_record) {
  
  // do nothing if buffer already flushed
  if (source_inoc_ID.size() == 0) {
    return;
  }
  
  // write to transmission record
  transmission_record << source_time << ",M," << source_host_ID << "," << mosq_ID << "," << inoc_ID << "," << source_inoc_ID[0];
  for (int i = 1; i < source_inoc_ID.size(); ++i) {
    transmission_record << ";" << source_inoc_ID[i];
  }
  transmission_record << "," << source_deme + 1 << "\n";
  
  // clear source_inoc_ID to record that buffer has been flushed
  source_inoc_ID.clear();
}

//------------------------------------------------
// print summary
void Mosquito::print_summary() {
  
  print("---------------------------");
  
  print("mosq_ID:", mosq_ID);
  print("inoc_ID:", inoc_ID);
  print("infectious_on:", infectious_on);
  print("time_infectious:", time_infectious);
  print("time_death:", time_death);
  print("source_time:", source_time);
  print("source_deme:", source_deme);
  print("source_host_ID:", source_host_ID);
  print("source_inoc_ID:");
  print_vector(source_inoc_ID);
  
  print("");
}
