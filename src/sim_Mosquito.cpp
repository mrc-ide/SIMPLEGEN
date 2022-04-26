
#include "sim_Mosquito.h"

using namespace std;

//------------------------------------------------
// infection
void sim_Mosquito::infection(int t, int &next_infection_ID, sim_Host &host) {
  
  // set infection ID and increment
  infection_ID = next_infection_ID++;
  
  // record properties of host for writing to transmission record later
  source_time = t;
  source_deme = host.deme;
  source_host_ID = host.host_ID;
  
  if (source_infection_ID_vec.size() != 0) {
    Rcpp::stop("error in sim_Mosquito::infection(), source_infection_ID_vec already populated");
  }
  
  for (int i = 0; i < host.infection_ID_vec.size(); ++i) {
    if (host.infection_active[i]) {
      source_infection_ID_vec.push_back(host.infection_ID_vec[i]);
    }
  }
  
}

//------------------------------------------------
// write buffered info to transmission record and clear buffer
void sim_Mosquito::write_buffer(ofstream &transmission_record) {
  
  if (source_infection_ID_vec.size() != 0) {
    
    transmission_record << source_time << ",2," << source_host_ID << "," << mosquito_ID << "," << infection_ID << "," << source_infection_ID_vec[0];
    for (int i = 1; i < source_infection_ID_vec.size(); ++i) {
      transmission_record << ";" << source_infection_ID_vec[i];
    }
    transmission_record << "," << source_deme + 1 << "\n";
    source_infection_ID_vec.clear();
    
  }
}

//------------------------------------------------
// set mosquito ID and increment
void sim_Mosquito::set_mosquito_ID(int &mosquito_ID) {
  this->mosquito_ID = mosquito_ID++;
}

//------------------------------------------------
// print status
void sim_Mosquito::print_status() {
  print("mosquito_ID:", mosquito_ID);
  print("infection_ID:", infection_ID);
  print("source_time:", source_time);
  print("source_deme:", source_deme);
  print("source_host_ID:", source_host_ID);
  print("source_infection_ID_vec:");
  print_vector(source_infection_ID_vec);
}
