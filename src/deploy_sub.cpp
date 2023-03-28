
#include <cpp11.hpp>

// This file contains a series of functions that simply redirect to functions
// within sub-directories. This is to help keep code tidy within the src folder

// include various files within sub-directories. Note that we are including the
// body files here and not the header files
#include "sim_indiv/sim_indiv.cpp"
#include "sim_indiv/Parameters.cpp"
#include "sim_indiv/Dispatcher.cpp"
#include "sim_indiv/Host.cpp"
#include "sim_indiv/Host_inoc.cpp"
#include "sim_indiv/Mosquito.cpp"
#include "sim_indiv/Mosquito_pop.cpp"
#include "sim_indiv/Host_pop.cpp"
#include "sim_indiv/Sweep.cpp"
#include "sim_indiv/Survey.cpp"
#include "prune_record/prune_main.cpp"

using namespace std;

//------------------------------------------------
[[cpp11::register]]
cpp11::list sim_indiv_deploy(cpp11::list args, cpp11::list args_progress) {
  return sim_indiv(args, args_progress);
}

//------------------------------------------------
[[cpp11::register]]
cpp11::list prune_transmission_record_deploy(cpp11::list args) {
  return prune_transmission_record(args);
}
