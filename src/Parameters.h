
#pragma once

#include "misc_v8.h"

#include <list>

//------------------------------------------------
// class defining all model parameters
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scalar epi parameters
  static double a, p, mu, treatment_seeking_alpha, treatment_seeking_beta;
  static int u, v, g, max_inoculations;
  
  // epi distributions
  static std::vector<double> prob_infection;
  static int n_prob_infection;
  static double max_prob_infection;
  static std::vector<double> prob_acute;
  static int n_prob_acute;
  static std::vector<double> prob_AC;
  static int n_prob_AC;
  static std::vector<std::vector<double>> duration_acute;
  static int n_duration_acute;
  static std::vector<std::vector<double>> duration_chronic;
  static int n_duration_chronic;
  static std::vector<std::vector<double>> detectability_microscopy_acute;
  static int n_detectability_microscopy_acute;
  static std::vector<std::vector<double>> detectability_microscopy_chronic;
  static int n_detectability_microscopy_chronic;
  static std::vector<std::vector<double>> detectability_PCR_acute;
  static int n_detectability_PCR_acute;
  static std::vector<std::vector<double>> detectability_PCR_chronic;
  static int n_detectability_PCR_chronic;
  static std::vector<std::vector<double>> time_treatment_acute;
  static int n_time_treatment_acute;
  static std::vector<std::vector<double>> time_treatment_chronic;
  static int n_time_treatment_chronic;
  static std::vector<double> duration_prophylactic;
  static int n_duration_prophylactic;
  static std::vector<std::vector<double>> infectivity_acute;
  static int n_infectivity_acute;
  static std::vector<std::vector<double>> infectivity_chronic;
  static int n_infectivity_chronic;
  static double max_infectivity;
  
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
  static bool save_transmission_record, output_daily_counts, output_age_distributions, silent;
  static std::string transmission_record_location;
  static std::vector<int> output_age_times;
  
  // misc parameters
  static double prob_mosq_death;  // daily probability of mosquito death
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  
  // methods
  void load_epi_params(double a, double p, double mu,
                       int u, int v, int g,
                       std::vector<double> prob_infection,
                       std::vector<double> prob_acute,
                       std::vector<double> prob_AC,
                       std::vector<std::vector<double>> duration_acute,
                       std::vector<std::vector<double>> duration_chronic,
                       std::vector<std::vector<double>> detectability_microscopy_acute,
                       std::vector<std::vector<double>> detectability_microscopy_chronic,
                       std::vector<std::vector<double>> detectability_PCR_acute,
                       std::vector<std::vector<double>> detectability_PCR_chronic,
                       std::vector<std::vector<double>> time_treatment_acute,
                       std::vector<std::vector<double>> time_treatment_chronic,
                       double treatment_seeking_alpha,
                       double treatment_seeking_beta,
                       std::vector<double> duration_prophylactic,
                       std::vector<std::vector<double>> infectivity_acute,
                       std::vector<std::vector<double>> infectivity_chronic,
                       int max_inoculations);
  
  void load_deme_params(std::vector<int> H_init,
                        std::vector<int> seed_infections,
                        std::vector<int> M);
  
  void load_demog_params(std::vector<double> life_table,
                         std::vector<double> age_death,
                         std::vector<double> age_stable);
  
  void load_run_params(int max_time,
                       bool save_transmission_record,
                       std::string transmission_record_location,
                       bool output_daily_counts,
                       bool output_age_distributions,
                       std::vector<int> output_age_times,
                       bool silent);
  
  void summary();
  
};

