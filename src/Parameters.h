
#pragma once

#include "misc_v8.h"

#include <list>

//------------------------------------------------
// enumerate sampling strategy methods
enum Case_detection {active, passive};
enum Diagnosis {microscopy, PCR};

//------------------------------------------------
// class defining all model parameters
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scalar epi parameters
  static double a, p, mu;
  static int u, v, g, max_inoculations;
  
  // state transition probabilities and durations
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
  
  // detectability
  static std::vector<std::vector<double>> detectability_microscopy_acute;
  static int n_detectability_microscopy_acute;
  static std::vector<std::vector<double>> detectability_microscopy_chronic;
  static int n_detectability_microscopy_chronic;
  static std::vector<std::vector<double>> detectability_PCR_acute;
  static int n_detectability_PCR_acute;
  static std::vector<std::vector<double>> detectability_PCR_chronic;
  static int n_detectability_PCR_chronic;
  
  // treatment
  static double treatment_seeking_mean, treatment_seeking_sd;
  static std::vector<std::vector<double>> time_treatment_acute;
  static int n_time_treatment_acute;
  static std::vector<std::vector<double>> time_treatment_chronic;
  static int n_time_treatment_chronic;
  static std::vector<double> duration_prophylactic;
  static int n_duration_prophylactic;
  
  // infectivity
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
  
  // sampling strategy parameters
  static bool obtain_samples;
  static std::vector<int> ss_time;
  static std::vector<int> ss_deme;
  static std::vector<Case_detection> ss_case_detection;
  static std::vector<Diagnosis> ss_diagnosis;
  static std::vector<int> ss_n;
  
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
  void load_epi_params(Rcpp::List args);
  void load_deme_params(Rcpp::List args);
  void load_demog_params(Rcpp::List args);
  void load_sampling_params(Rcpp::List args);
  void load_run_params(Rcpp::List args);
  
  void summary();
  
};

