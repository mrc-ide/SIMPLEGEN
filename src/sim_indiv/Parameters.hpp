
#pragma once

#include <cpp11.hpp>
//#include <list>
#include "../misc_v17.hpp"
#include "../Sampler_v6.hpp"
#include "../utils.hpp"


//------------------------------------------------
// enumerate sampling strategy methods
//enum Model_state {Model_S, Model_E, Model_A, Model_C, Model_P, Model_H, Model_Sv, Model_Ev, Model_Iv, Model_M};
//enum Measure {Measure_count, Measure_prevalence, Measure_incidence, Measure_EIR};
//enum Sampling {Sampling_none, Sampling_ACD, Sampling_PCD};
//enum Diagnostic {Diagnostic_true, Diagnostic_microscopy, Diagnostic_PCR};
//enum Case_detection {Case_active, Case_passive};

//------------------------------------------------
// class defining all transmission model parameters. Also contains all Sampler
// objects, and provides methods for taking draws from these distributions
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // scalar epi parameters
  double a, p;
  int u, v, g, max_inoculations;
  
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
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  std::vector<Sampler> sampler_duration_acute;
  std::vector<Sampler> sampler_duration_chronic;
  std::vector<Sampler> sampler_time_treatment_acute;
  std::vector<Sampler> sampler_time_treatment_chronic;
  std::vector<Sampler> sampler_duration_prophylactic;
  
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
  int max_age;
  
  
  // sampling parameters
  bool any_daily_outputs;
  bool any_sweep_outputs;
  bool any_survey_outputs;
  cpp11::data_frame daily_df;
  cpp11::data_frame sweeps_df;
  cpp11::data_frame surveys_df;
  
  // run parameters
  int max_time;
  bool save_transmission_record;
  std::string transmission_record_location;
  bool pb_markdown;
  bool silent;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  
  // methods
  void load_params(cpp11::list args);
  void load_model_params(cpp11::list args);
  void load_deme_params(cpp11::list args);
  void load_migration_params(cpp11::list args);
  void load_demog_params(cpp11::list args);
  void load_sampling_params(cpp11::list args);
  void load_run_params(cpp11::list args);
  
  void define_samplers();
  int draw_age_stable();
  int draw_age_death();
  int draw_duration_acute(int n);
  int draw_duration_chronic(int n);
  int draw_time_treatment_acute(int n);
  int draw_time_treatment_chronic(int n);
  int draw_duration_prophylactic(int n);
  double get_prob_infection(int n);
  double get_prob_acute(int n);
  double get_prob_AC(int n);
  
  double get_detectability_microscopy_acute(int t, int t0, int n);
  double get_detectability_microscopy_chronic(int t, int t0, int n);
  double get_detectability_PCR_acute(int t, int t0, int n);
  double get_detectability_PCR_chronic(int t, int t0, int n);
  double get_infectivity_acute(int t, int t0, int n);
  double get_infectivity_chronic(int t, int t0, int n);
  
  void summary();
};

