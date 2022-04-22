
#pragma once

#include "misc_v14.h"

#include <list>

//------------------------------------------------
// enumerate sampling strategy methods
enum Model_state {Model_S, Model_E, Model_A, Model_C, Model_P, Model_H, Model_Sv, Model_Ev, Model_Iv, Model_M};
enum Measure {Measure_count, Measure_prevalence, Measure_incidence, Measure_EIR};
enum Sampling {Sampling_none, Sampling_ACD, Sampling_PCD};
enum Diagnostic {Diagnostic_true, Diagnostic_microscopy, Diagnostic_PCR};
enum Case_detection {Case_active, Case_passive};

//------------------------------------------------
// class defining all model parameters
class sim_Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scalar epi parameters
  double a, p, mu;
  int u, v, g, max_infections;
  
  // state transition probabilities and durations
  std::vector<double> prob_infection;
  std::vector<double> prob_acute;
  std::vector<double> prob_AC;
  std::vector<std::vector<double>> duration_acute;
  std::vector<std::vector<double>> duration_chronic;
  int n_prob_infection;
  int n_prob_acute;
  int n_prob_AC;
  int n_duration_acute;
  int n_duration_chronic;
  double max_prob_infection;
  
  // detectability
  std::vector<std::vector<double>> detectability_microscopy_acute;
  std::vector<std::vector<double>> detectability_microscopy_chronic;
  std::vector<std::vector<double>> detectability_PCR_acute;
  std::vector<std::vector<double>> detectability_PCR_chronic;
  int n_detectability_microscopy_acute;
  int n_detectability_microscopy_chronic;
  int n_detectability_PCR_acute;
  int n_detectability_PCR_chronic;
  
  // treatment
  double treatment_seeking_mean;
  double treatment_seeking_sd;
  std::vector<std::vector<double>> time_treatment_acute;
  std::vector<std::vector<double>> time_treatment_chronic;
  std::vector<std::vector<double>> duration_prophylactic;
  int n_time_treatment_acute;
  int n_time_treatment_chronic;
  int n_duration_prophylactic;
  
  // onward infectivity
  std::vector<std::vector<double>> infectivity_acute;
  std::vector<std::vector<double>> infectivity_chronic;
  int n_infectivity_acute;
  int n_infectivity_chronic;
  double max_infectivity;
  
  // misc epi parameters
  double prob_mosq_death;  // daily probability of mosquito death
  
  // deme parameters
  std::vector<int> H_init;
  std::vector<int> seed_infections;
  std::vector<int> M;
  int n_demes;
  
  // migration parameters
  std::vector<std::vector<double>> mig_mat;
  
  // demog parameters
  std::vector<double> life_table;
  int n_life_table;
  std::vector<double> age_death;
  std::vector<double> age_stable;
  
  // sampling parameters
  bool any_daily_outputs;
  std::map<std::pair<int, int>, std::vector<int>> daily_map;
  std::vector<bool> daily_flag_deme;
  std::vector<int> daily_deme;
  std::vector<Measure> daily_measure;
  std::vector<Model_state> daily_state;
  std::vector<Diagnostic> daily_diagnostic;
  std::vector<int> daily_age_min;
  std::vector<int> daily_age_max;
  int n_daily_outputs;
  
  bool any_sweep_outputs;
  std::vector<int> sweep_time;
  std::vector<int> sweep_time_ordered;
  std::vector<int> sweep_deme;
  std::vector<Measure> sweep_measure;
  std::vector<Model_state> sweep_state;
  std::vector<Diagnostic> sweep_diagnostic;
  std::vector<int> sweep_age_min;
  std::vector<int> sweep_age_max;
  int n_sweep_outputs;
  
  bool any_survey_outputs;
  std::vector<int> surveys_t_start;
  std::vector<int> surveys_t_end;
  std::vector<int> surveys_deme;
  std::vector<Measure> surveys_measure;
  std::vector<Sampling> surveys_sampling;
  std::vector<Diagnostic> surveys_diagnostic;
  std::vector<int> surveys_age_min;
  std::vector<int> surveys_age_max;
  std::vector<double> surveys_sample_size;
  std::vector<int> surveys_n_days;
  int n_survey_outputs;
  
  std::vector<int> surveys_expanded_study_ID;
  std::vector<int> surveys_expanded_sampling_time;
  std::vector<int> surveys_expanded_reporting_time;
  
  // run parameters
  int max_time;
  bool save_transmission_record;
  std::string transmission_record_location;
  bool pb_markdown;
  bool silent;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  sim_Parameters() {};
  
  // methods
  void load_params(Rcpp::List args);
  void load_model_params(Rcpp::List args);
  void load_deme_params(Rcpp::List args);
  void load_migration_params(Rcpp::List args);
  void load_demog_params(Rcpp::List args);
  void load_sampling_params(Rcpp::List args);
  void load_sampling_params_daily(Rcpp::List args);
  void load_sampling_params_sweep(Rcpp::List args);
  void load_sampling_params_survey(Rcpp::List args);
  void load_run_params(Rcpp::List args);
  
  void summary();
  void print_daily_maps();
  
};

