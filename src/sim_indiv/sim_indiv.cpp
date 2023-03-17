
#include "sim_indiv.hpp"

using namespace std;

//------------------------------------------------
void sim_indiv(cpp11::list args, cpp11::list args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define parameters object and load values from arguments
  Parameters params;
  params.load_params(args);
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.init(params);
  dispatcher.run_simulation(args_progress);
  /*
  // copy update_progress function from package
  cpp11::function update_progress = cpp11::package("SIMPLEGEN")["update_progress"];
  
  
  double z = 100.0;
  for (int i = 0; i < max_time; ++i) {
    update_progress(args_progress, "pb_sim", i, max_time, true);
    for (int j = 0; j < 1e8; ++j) {
      z = pow(z, 0.9);
    }
  }
  */
  // end timer
  chrono_timer(t1);
  
}
