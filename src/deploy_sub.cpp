#include <cpp11.hpp>

// This file contains a series of functions that simply redirect to functions
// within sub-directories. This is to help keep code tidy within the src folder

// include various files within sub-directories. Note that we are including the
// body files here and not the header files
#include "indiv_sim/sim_indiv.cpp"

using namespace std;

//------------------------------------------------
[[cpp11::register]]
void sim_indiv_deploy(cpp11::list args_progress) {
  sim_indiv(args_progress);
}
