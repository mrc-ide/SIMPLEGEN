
#include "Dispatcher.h"

using namespace std;

//------------------------------------------------
// constructor
Dispatcher::Dispatcher() {
  
  // objects for sampling from probability distributions
  //sampler_age_stable = Sampler(param_ptr->age_stable, 1000);
  //sampler_age_death = Sampler(param_ptr->age_death, 1000);
  //sampler_duration_infection = Sampler(param_ptr->duration_infection, 1000);
  
}

//------------------------------------------------
// run main simulation
void Dispatcher::run_simulation() {
  
  
  print_vector(prob_acute);
  foobar();
  
}

