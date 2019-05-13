
#include "main.h"
#include "misc_v5.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// main function (when not run using Rcpp)
#ifndef RCPP_ACTIVE
int main(int argc, const char * argv[]) {
  
  // run simulation
  indiv_sim_cpp();
  
}
#endif

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // do something
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("foo") = -9);
}
#else
int indiv_sim_cpp() {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // do something
  
  // end timer
  chrono_timer(t1);
  
  return 0;
}
#endif
