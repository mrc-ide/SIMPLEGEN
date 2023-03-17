
#include "Parameters.hpp"

using namespace std;


//------------------------------------------------
// load all parameter values
void Parameters::load_params(cpp11::list args) {
  
  // load all parameter values
  load_model_params(args);
  load_deme_params(args);
  load_migration_params(args);
  load_demog_params(args);
  //load_sampling_params_daily(args);
  //load_sampling_params_sweep(args);
  //load_sampling_params_survey(args);
  load_run_params(args);
  
}

//------------------------------------------------
// load epi model parameter values
void Parameters::load_model_params(cpp11::list args) {
  
  // scalar epi parameters
  a = cpp_to_double(args["a"]);
  p = cpp_to_double(args["p"]);
  mu = cpp_to_double(args["mu"]);
  u = cpp_to_int(args["u"]);
  v = cpp_to_int(args["v"]);
  g = cpp_to_int(args["g"]);
  max_inoculations = cpp_to_int(args["max_inoculations"]);
  
  // state transition probabilities and durations
  prob_infection = cpp_to_vector_double(args["prob_infection"]);
  prob_acute = cpp_to_vector_double(args["prob_acute"]);
  prob_AC = cpp_to_vector_double(args["prob_AC"]);
  duration_acute = cpp_to_matrix_double(args["duration_acute"]);
  duration_chronic = cpp_to_matrix_double(args["duration_chronic"]);
  n_prob_infection = int(prob_infection.size());
  n_prob_acute = int(prob_acute.size());
  n_prob_AC = int(prob_AC.size());
  n_duration_acute = int(duration_acute.size());
  n_duration_chronic = int(duration_chronic.size());
  
  // detectability
  detectability_microscopy_acute = cpp_to_matrix_double(args["detectability_microscopy_acute"]);
  detectability_microscopy_chronic = cpp_to_matrix_double(args["detectability_microscopy_chronic"]);
  detectability_PCR_acute = cpp_to_matrix_double(args["detectability_PCR_acute"]);
  detectability_PCR_chronic = cpp_to_matrix_double(args["detectability_PCR_chronic"]);
  n_detectability_microscopy_acute = int(detectability_microscopy_acute.size());
  n_detectability_microscopy_chronic = int(detectability_microscopy_chronic.size());
  n_detectability_PCR_acute = int(detectability_PCR_acute.size());
  n_detectability_PCR_chronic = int(detectability_PCR_chronic.size());
  
  // treatment
  treatment_seeking_mean = cpp_to_double(args["treatment_seeking_mean"]);
  treatment_seeking_sd = cpp_to_double(args["treatment_seeking_sd"]);
  time_treatment_acute = cpp_to_matrix_double(args["time_treatment_acute"]);
  time_treatment_chronic = cpp_to_matrix_double(args["time_treatment_chronic"]);
  duration_prophylactic = cpp_to_matrix_double(args["duration_prophylactic"]);
  n_time_treatment_acute = int(time_treatment_acute.size());
  n_time_treatment_chronic = int(time_treatment_chronic.size());
  n_duration_prophylactic = int(duration_prophylactic.size());
  
  // infectivity
  infectivity_acute = cpp_to_matrix_double(args["infectivity_acute"]);
  infectivity_chronic = cpp_to_matrix_double(args["infectivity_chronic"]);
  n_infectivity_acute = int(infectivity_acute.size());
  n_infectivity_chronic = int(infectivity_chronic.size());
  
  // get max prob_infection
  max_prob_infection = max_vec(prob_infection);
  
  // get max infectivity over all distributions, both acute and chronic
  max_infectivity = std::max(max_mat(infectivity_acute),
                             max_mat(infectivity_chronic));
}

//------------------------------------------------
// load deme parameter values
void Parameters::load_deme_params(cpp11::list args) {
  
  // distributions
  H_init = cpp_to_vector_int(args["H"]);
  seed_infections = cpp_to_vector_int(args["seed_infections"]);
  M = cpp_to_vector_int(args["M"]);
  n_demes = int(H_init.size());
  
}

//------------------------------------------------
// load migration parameter values
void Parameters::load_migration_params(cpp11::list args) {
  
  // migration matrix
  mig_mat = cpp_to_matrix_double(args["mig_mat"]);
  
}

//------------------------------------------------
// load demography parameter values
void Parameters::load_demog_params(cpp11::list args) {
  
  // distributions
  life_table = cpp_to_vector_double(args["life_table"]);
  n_life_table = int(life_table.size());
  age_death = cpp_to_vector_double(args["age_death"]);
  age_stable = cpp_to_vector_double(args["age_stable"]);
  
}
/*
//------------------------------------------------
// load daily sampling strategy parameter values
void sim_Parameters::load_sampling_params_daily(Rcpp::List args) {
  
  // return if no daily outputs
  any_daily_outputs = args["any_daily_outputs"];
  if (!any_daily_outputs) {
    n_daily_outputs = 0;
    return;
  }
  
  // extract dataframe
  Rcpp::List daily_df = args["daily"];
  
  // extract columns into vectors
  daily_deme = cpp_to_vector_int(daily_df["deme"]);
  n_daily_outputs = daily_deme.size();
  
  vector<int> daily_measure_int = cpp_to_vector_int(daily_df["measure"]);
  daily_measure = vector<Measure>(daily_measure_int.size());
  for (size_t i = 0; i < daily_measure_int.size(); ++i) {
    daily_measure[i] = static_cast<Measure>(daily_measure_int[i]);
  }
  
  vector<int> daily_state_int = cpp_to_vector_int(daily_df["state"]);
  daily_state = vector<Model_state>(daily_state_int.size());
  for (size_t i = 0; i < daily_state_int.size(); ++i) {
    daily_state[i] = static_cast<Model_state>(daily_state_int[i]);
  }
  
  vector<int> daily_diagnostic_int = cpp_to_vector_int(daily_df["diagnostic"]);
  daily_diagnostic = vector<Diagnostic>(daily_diagnostic_int.size());
  for (size_t i = 0; i < daily_diagnostic_int.size(); ++i) {
    daily_diagnostic[i] = static_cast<Diagnostic>(daily_diagnostic_int[i]);
  }
  
  daily_age_min = cpp_to_vector_int(daily_df["age_min"]);
  daily_age_max = cpp_to_vector_int(daily_df["age_max"]);
  
  // Create a map to assist in working out if a host is required to produce
  // daily output. The key to the map is a deme-&-age combination (1-year age
  // groups). The value associated with this key is a list of all output indices
  // (i.e. rows of the daily dataframe) that apply to this group.
  //
  // For example, if output row 2 needs to be the prevalence in 0-5 year olds
  // in deme 1, then the map at key {1,0} and {1,1} etc. up to {1,5} will all
  // contain the value 2 (the output index), along with any other output indices
  // that apply to this group.
  //
  // A simpler vector object also exists just for checking if a deme is required
  // in any outputs. This avoids looping through hosts in demes that are not
  // required in any output.
  //
  // The function print_daily_maps() can be used to look at these objects.
  
  daily_flag_deme = vector<bool>(n_demes, false);
  for (int i = 0; i < n_daily_outputs; i++) {
    if (daily_deme[i] == -1) {
      for (int k = 0; k < n_demes; ++k) {
        daily_flag_deme[k] = true;
        for (int j = daily_age_min[i]; j < (daily_age_max[i] + 1); ++j) {
          daily_map[make_pair(k, j)].push_back(i);
        }
      }
    } else {
      int this_deme = daily_deme[i];
      daily_flag_deme[this_deme] = true;
      for (int j = daily_age_min[i]; j < (daily_age_max[i] + 1); ++j) {
        daily_map[make_pair(this_deme, j)].push_back(i);
      }
    }
  }
  
}

//------------------------------------------------
// print maps that assist in working out which hosts are required to produce
// daily output
void sim_Parameters::print_daily_maps() {
  
  // individual-level map
  print("DAILY MAP:");
  for (auto it = daily_map.cbegin(); it != daily_map.cend(); ++it) {
    pair<int, int> map_key = it->first;
    Rcpp::Rcout << "deme" << map_key.first << " age" << map_key.second << ": ";
    for (size_t i = 0; i < it->second.size(); ++i) {
      Rcpp::Rcout << it->second[i] << " ";
    }
    Rcpp::Rcout << "\n";
  }
  
  // deme-level mep
  print("DAILY DEME MAP:");
  print_vector(daily_flag_deme);
  
}

//------------------------------------------------
// load population sweep sampling strategy parameter values
void sim_Parameters::load_sampling_params_sweep(Rcpp::List args) {
  
  // return if no sweep outputs
  any_sweep_outputs = args["any_sweep_outputs"];
  if (!any_sweep_outputs) {
    n_sweep_outputs = 0;
    return;
  }
  
  // extract dataframe
  Rcpp::List sweep_df = args["sweeps"];
  
  // extract columns into vectors
  sweep_time = cpp_to_vector_int(sweep_df["time"]);
  n_sweep_outputs = sweep_time.size();
  sweep_time_ordered = cpp_to_vector_int(args["sweep_time_ordered"]);
  sweep_deme = cpp_to_vector_int(sweep_df["deme"]);
  
  vector<int> sweep_measure_int = cpp_to_vector_int(sweep_df["measure"]);
  sweep_measure = vector<Measure>(sweep_measure_int.size());
  for (size_t i = 0; i < sweep_measure_int.size(); ++i) {
    sweep_measure[i] = static_cast<Measure>(sweep_measure_int[i]);
  }
  
  vector<int> sweep_state_int = cpp_to_vector_int(sweep_df["state"]);
  sweep_state = vector<Model_state>(sweep_state_int.size());
  for (size_t i = 0; i < sweep_state_int.size(); ++i) {
    sweep_state[i] = static_cast<Model_state>(sweep_state_int[i]);
  }
  
  vector<int> sweep_diagnostic_int = cpp_to_vector_int(sweep_df["diagnostic"]);
  sweep_diagnostic = vector<Diagnostic>(sweep_diagnostic_int.size());
  for (size_t i = 0; i < sweep_diagnostic_int.size(); ++i) {
    sweep_diagnostic[i] = static_cast<Diagnostic>(sweep_diagnostic_int[i]);
  }
  
  sweep_age_min = cpp_to_vector_int(sweep_df["age_min"]);
  sweep_age_max = cpp_to_vector_int(sweep_df["age_max"]);
  
}

//------------------------------------------------
// load survey sampling strategy parameter values
void sim_Parameters::load_sampling_params_survey(Rcpp::List args) {
  
  // return if no survey outputs
  any_survey_outputs = args["any_survey_outputs"];
  if (!any_survey_outputs) {
    n_survey_outputs = 0;
    return;
  }
  
  // extract dataframes
  Rcpp::List surveys_df = args["surveys"];
  Rcpp::List surveys_expanded_df = args["surveys_expanded"];
  
  // extract surveys info
  surveys_t_start = cpp_to_vector_int(surveys_df["t_start"]);
  surveys_t_end = cpp_to_vector_int(surveys_df["t_end"]);
  surveys_deme = cpp_to_vector_int(surveys_df["deme"]);
  n_survey_outputs = surveys_deme.size();
  
  vector<int> surveys_measure_int = cpp_to_vector_int(surveys_df["measure"]);
  surveys_measure = vector<Measure>(surveys_measure_int.size());
  for (size_t i = 0; i < n_survey_outputs; ++i) {
    surveys_measure[i] = static_cast<Measure>(surveys_measure_int[i]);
  }
  
  vector<int> surveys_sampling_int = cpp_to_vector_int(surveys_df["sampling"]);
  surveys_sampling = vector<Sampling>(surveys_sampling_int.size());
  for (size_t i = 0; i < n_survey_outputs; ++i) {
    surveys_sampling[i] = static_cast<Sampling>(surveys_sampling_int[i]);
  }
  
  vector<int> surveys_diagnostic_int = cpp_to_vector_int(surveys_df["diagnostic"]);
  surveys_diagnostic = vector<Diagnostic>(surveys_diagnostic_int.size());
  for (size_t i = 0; i < n_survey_outputs; ++i) {
    surveys_diagnostic[i] = static_cast<Diagnostic>(surveys_diagnostic_int[i]);
  }
  
  surveys_age_min = cpp_to_vector_int(surveys_df["age_min"]);
  surveys_age_max = cpp_to_vector_int(surveys_df["age_max"]);
  surveys_sample_size = cpp_to_vector_double(surveys_df["sample_size"]);
  surveys_n_days = cpp_to_vector_int(surveys_df["n_days"]);
  
  // extract expanded surveys info
  surveys_expanded_study_ID = cpp_to_vector_int(surveys_expanded_df["study_ID"]);
  surveys_expanded_sampling_time = cpp_to_vector_int(surveys_expanded_df["sampling_time"]);
  surveys_expanded_reporting_time = cpp_to_vector_int(surveys_expanded_df["reporting_time"]);
  
}*/

//------------------------------------------------
// load run parameter values
void Parameters::load_run_params(cpp11::list args) {
  
  // load values
  max_time = cpp_to_int(args["max_time"]);
  save_transmission_record = cpp_to_bool(args["save_transmission_record"]);
  transmission_record_location = cpp_to_string(args["transmission_record_location"]);
  silent = cpp_to_bool(args["silent"]);
  pb_markdown = cpp_to_bool(args["pb_markdown"]);
  
}

//------------------------------------------------
// create all sampler objects
void Parameters::define_samplers() {
  
  // define number of random values to generate in one go (see Sampler class
  // header for details)
  int sampler_draws = 1000;
  
  // create samplers
  sampler_age_stable = Sampler(age_stable, sampler_draws);
  sampler_age_death = Sampler(age_death, sampler_draws);
  sampler_duration_acute = make_sampler_vec(duration_acute, sampler_draws);
  sampler_duration_chronic = make_sampler_vec(duration_chronic, sampler_draws);
  sampler_time_treatment_acute = make_sampler_vec(time_treatment_acute, sampler_draws);
  sampler_time_treatment_chronic = make_sampler_vec(time_treatment_chronic, sampler_draws);
  sampler_duration_prophylactic = make_sampler_vec(duration_prophylactic, sampler_draws);
  
}

//------------------------------------------------
// print summary of parameters
void Parameters::summary() {
  
  // print epi scalars
  print("a:", a);
  print("p:", p);
  print("mu:", mu);
  print("u:", u);
  print("v:", v);
  print("g:", g);
  print("treatment_seeking_mean:", treatment_seeking_mean);
  print("treatment_seeking_sd:", treatment_seeking_sd);
  print("max_inoculations:", max_inoculations);
  
  // print epi distributions
  print("prob_infection:");
  print_vector(prob_infection);
  print("prob_acute:");
  print_vector(prob_acute);
  print("prob_AC:");
  print_vector(prob_AC);
  print("duration_acute:");
  print_matrix(duration_acute);
  print("duration_chronic:");
  print_matrix(duration_chronic);
  print("detectability_microscopy_acute:");
  print_matrix(detectability_microscopy_acute);
  print("detectability_microscopy_chronic:");
  print_matrix(detectability_microscopy_chronic);
  print("detectability_PCR_acute:");
  print_matrix(detectability_PCR_acute);
  print("detectability_PCR_chronic:");
  print_matrix(detectability_PCR_chronic);
  print("time_treatment_acute:");
  print_matrix(time_treatment_acute);
  print("time_treatment_chronic:");
  print_matrix(time_treatment_chronic);
  print("duration_prophylactic:");
  print_matrix(duration_prophylactic);
  print("infectivity_acute:");
  print_matrix(infectivity_acute);
  print("infectivity_chronic:");
  print_matrix(infectivity_chronic);
  
  // print deme parameters
  print("H_init:");
  print_vector(H_init);
  print("seed_infections:");
  print_vector(seed_infections);
  print("M:");
  print_vector(M);
  
  // print demog parameters
  print("life_table:");
  print_vector(life_table);
  print("age_death:");
  print_vector(age_death);
  print("age_stable:");
  print_vector(age_stable);
  
  // print run scalars
  print("max_time:", max_time);
  print("save_transmission_record:", save_transmission_record);
  print("transmission_record_location:", transmission_record_location);
  print("pb_markdown:", pb_markdown);
  print("silent:", silent);
  
}
