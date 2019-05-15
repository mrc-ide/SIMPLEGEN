
#pragma once

#include "misc_v6.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <list>

//------------------------------------------------
// class defining all model parameters
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scalar parameters
  double a, p, mu, prob_AC;
  int u, v, g, max_innoculations;
  
  // distributions
  std::vector<double> prob_acute;
  std::vector<double> prob_infection;
  std::vector<double> infectivity_acute;
  std::vector<double> infectivity_chronic;
  std::vector<std::vector<double>> duration_acute;
  std::vector<std::vector<double>> duration_chronic;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  
  // methods
  void load_values(double a, double p, double mu, double prob_AC,
                   int u, int v, int g, int max_innoculations,
                   std::vector<double> prob_acute,
                   std::vector<double> prob_infection,
                   std::vector<double> infectivity_acute,
                   std::vector<double> infectivity_chronic,
                   std::vector<std::vector<double>> duration_acute,
                   std::vector<std::vector<double>> duration_chronic);
  void summary();
  
};

