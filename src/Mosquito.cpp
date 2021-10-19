
#include "Mosquito.h"

using namespace std;

//------------------------------------------------
// infection
void Mosquito::infection(int t, int &next_infection_ID, Host &host) {
  
  // set infection ID and increment
  infection_ID = next_infection_ID++;
  
  // record properties of host for writing to transmission record later
  source_time = t;
  source_deme = host.deme;
  source_host_ID = host.host_ID;
  
  if (source_infection_ID_vec.size() != 0) {
    Rcpp::stop("error in Mosquito::infection(), source_infection_ID_vec already populated");
  }
  
  for (int i = 0; i < host.inoc_ID_vec.size(); ++i) {
    if (host.inoc_active[i]) {
      source_infection_ID_vec.push_back(host.inoc_ID_vec[i]);
    }
  }
  
}

//------------------------------------------------
// write buffered info to transmission record and clear buffer
void Mosquito::write_buffer(ofstream &transmission_record) {
  
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
void Mosquito::set_mosquito_ID(int &mosquito_ID) {
  this->mosquito_ID = mosquito_ID++;
}
