
#include "Dispatcher.h"

using namespace std;

//------------------------------------------------
// constructor
Dispatcher::Dispatcher() {
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(age_stable, 1000);
  sampler_age_death = Sampler(age_death, 1000);
  //sampler_duration_infection = Sampler(duration_infection, 1000);
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  /*
  schedule_death = vector<set<int>>(max_time+1);
  schedule_Eh_to_Ih = vector<vector<pair<int, int>>>(max_time+1);
  schedule_Ih_to_Sh = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_recovery = vector<vector<pair<int, int>>>(max_time+1);
  */
}

//------------------------------------------------
// run main simulation
void Dispatcher::run_simulation() {
  
  
  print_vector(prob_acute);
  foobar();
  
}

