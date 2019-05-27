
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
  static int n_prob_infection;
  static std::vector<std::vector<double>> infectivity_acute;
  static std::vector<std::vector<double>> infectivity_chronic;
  static double max_infectivity;
  static std::vector<std::vector<double>> duration_acute;
  static std::vector<std::vector<double>> duration_chronic;
  
  // deme parameters
  static std::vector<int> H_init;
  static std::vector<int> seed_infections;
  static std::vector<int> M;
  static int n_demes;
  
  // demog parameters
  static std::vector<double> life_table;
  static std::vector<double> age_death;
  static std::vector<double> age_stable;
  
  // run parameters
  static int max_time;
  static bool output_daily_counts, output_age_distributions, output_infection_history, silent;
  static std::vector<int> output_age_times;
  
  // misc parameters
  static double prob_mosq_death;  // daily probability of mosquito death
  //double prob_survive_extrinsic;  // probability of mosquito surviving the extrinsic incubation period
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  
  // methods
  void load_epi_params(double a, double p, double mu, double prob_AC,
                       int u, int v, int g, int max_innoculations,
                       std::vector<double> prob_acute,
                       std::vector<double> prob_infection,
                       std::vector<std::vector<double>> infectivity_acute,
                       std::vector<std::vector<double>> infectivity_chronic,
                       std::vector<std::vector<double>> duration_acute,
                       std::vector<std::vector<double>> duration_chronic);
  
  void load_deme_params(std::vector<int> H_init,
                        std::vector<int> seed_infections,
                        std::vector<int> M);
  
  void load_demog_params(std::vector<double> life_table,
                         std::vector<double> age_death,
                         std::vector<double> age_stable);
  
  void load_run_params(int max_time,
                       bool output_daily_counts,
                       bool output_age_distributions,
                       bool output_infection_history,
                       bool silent,
                       std::vector<int> output_age_times);
  
  void summary();
  
};

