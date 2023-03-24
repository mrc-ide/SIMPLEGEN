
#include "sim_indiv.hpp"

using namespace std;

//------------------------------------------------
cpp11::list sim_indiv(cpp11::list args, cpp11::list args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define parameters object and load values from arguments
  Parameters params;
  params.load_params(args);
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.init(params);
  dispatcher.run_simulation(args_progress);
  
  // end timer
  chrono_timer(t1);
  
  // return
  return dispatcher.tmp_output;
  
  //cpp11::writable::list tmp_output;
  //tmp_output.push_back({cpp11::literals::operator""_nm("foo", 0) = -9});
  //return tmp_output;
}
