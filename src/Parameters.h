
#pragma once

#include "misc_v6.h"

#include <list>

//------------------------------------------------
// class defining all model parameters
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scalar epi parameters
  static double a, p, mu, prob_AC;
  static int u, v, g, max_innoculations;
  
  // epi distributions
  static std::vector<double> prob_acute;
  static std::vector<double> prob_infection;
  static std::vector<double> infectivity_acute;
  static std::vector<double> infectivity_chronic;
  static std::vector<std::vector<double>> duration_acute;
  static std::vector<std::vector<double>> duration_chronic;
  
  // deme parameters
  static std::vector<int> H;
  static std::vector<int> seed_infections;
  static std::vector<int> M;
  
  // demog parameters
  static std::vector<double> life_table;
  static std::vector<double> age_death;
  static std::vector<double> age_stable;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  
  // methods
  void load_epi_params(double a, double p, double mu, double prob_AC,
                       int u, int v, int g, int max_innoculations,
                       std::vector<double> prob_acute,
                       std::vector<double> prob_infection,
                       std::vector<double> infectivity_acute,
                       std::vector<double> infectivity_chronic,
                       std::vector<std::vector<double>> duration_acute,
                       std::vector<std::vector<double>> duration_chronic);
  
  void load_deme_params(std::vector<int> H,
                        std::vector<int> seed_infections,
                        std::vector<int> M);
  
  void load_demog_params(std::vector<double> life_table,
                         std::vector<double> age_death,
                         std::vector<double> age_stable);
  
  void summary();
  
};

