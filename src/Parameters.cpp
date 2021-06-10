
#include "Parameters.h"

using namespace std;

//------------------------------------------------
// load all parameter values
void Parameters::load_params(Rcpp::List args) {
  
  // load all parameter values
  load_model_params(args);
  load_deme_params(args);
  load_migration_params(args);
  load_demog_params(args);
  load_sampling_params(args);
  load_run_params(args);
  
}

//------------------------------------------------
// load epi model parameter values
void Parameters::load_model_params(Rcpp::List args) {
  
  // scalar epi parameters
  a = rcpp_to_double(args["a"]);
  p = rcpp_to_double(args["p"]);
  mu = rcpp_to_double(args["mu"]);
  u = rcpp_to_int(args["u"]);
  v = rcpp_to_int(args["v"]);
  g = rcpp_to_int(args["g"]);
  max_inoculations = rcpp_to_int(args["max_inoculations"]);
  
  // state transition probabilities and durations
  prob_infection = rcpp_to_vector_double(args["prob_infection"]);
  n_prob_infection = int(prob_infection.size());
  prob_acute = rcpp_to_vector_double(args["prob_acute"]);
  n_prob_acute = int(prob_acute.size());
  prob_AC = rcpp_to_vector_double(args["prob_AC"]);
  n_prob_AC = int(prob_AC.size());
  duration_acute = rcpp_to_matrix_double(args["duration_acute"]);
  n_duration_acute = int(duration_acute.size());
  duration_chronic = rcpp_to_matrix_double(args["duration_chronic"]);
  n_duration_chronic = int(duration_chronic.size());
  
  // detectability
  detectability_microscopy_acute = rcpp_to_matrix_double(args["detectability_microscopy_acute"]);
  n_detectability_microscopy_acute = int(detectability_microscopy_acute.size());
  detectability_microscopy_chronic = rcpp_to_matrix_double(args["detectability_microscopy_chronic"]);
  n_detectability_microscopy_chronic = int(detectability_microscopy_chronic.size());
  detectability_PCR_acute = rcpp_to_matrix_double(args["detectability_PCR_acute"]);
  n_detectability_PCR_acute = int(detectability_PCR_acute.size());
  detectability_PCR_chronic = rcpp_to_matrix_double(args["detectability_PCR_chronic"]);
  n_detectability_PCR_chronic = int(detectability_PCR_chronic.size());
  
  // treatment
  treatment_seeking_mean = rcpp_to_double(args["treatment_seeking_mean"]);
  treatment_seeking_sd = rcpp_to_double(args["treatment_seeking_sd"]);
  time_treatment_acute = rcpp_to_matrix_double(args["time_treatment_acute"]);
  n_time_treatment_acute = int(time_treatment_acute.size());
  time_treatment_chronic = rcpp_to_matrix_double(args["time_treatment_chronic"]);
  n_time_treatment_chronic = int(time_treatment_chronic.size());
  duration_prophylactic = rcpp_to_matrix_double(args["duration_prophylactic"]);
  n_duration_prophylactic = int(duration_prophylactic.size());
  
  // infectivity
  infectivity_acute = rcpp_to_matrix_double(args["infectivity_acute"]);
  n_infectivity_acute = int(infectivity_acute.size());
  infectivity_chronic = rcpp_to_matrix_double(args["infectivity_chronic"]);
  n_infectivity_chronic = int(infectivity_chronic.size());
  
  // get max prob_infection
  max_prob_infection = max(prob_infection);
  
  // get max infectivity over all distributions
  max_infectivity = 0.0;
  for (unsigned int i = 0; i < infectivity_acute.size(); ++i) {
    max_infectivity = (max(infectivity_acute[i]) > max_infectivity) ? max(infectivity_acute[i]) : max_infectivity;
  }
  for (unsigned int i = 0; i < infectivity_chronic.size(); ++i) {
    max_infectivity = (max(infectivity_chronic[i]) > max_infectivity) ? max(infectivity_chronic[i]) : max_infectivity;
  }
  
  // misc parameters
  prob_mosq_death = 1 - exp(-mu);  // daily probability of mosquito death
  
}

//------------------------------------------------
// load deme parameter values
void Parameters::load_deme_params(Rcpp::List args) {
  
  // distributions
  H_init = rcpp_to_vector_int(args["H"]);
  seed_infections = rcpp_to_vector_int(args["seed_infections"]);
  M = rcpp_to_vector_int(args["M"]);
  n_demes = int(H_init.size());
  
}

//------------------------------------------------
// load migration parameter values
void Parameters::load_migration_params(Rcpp::List args) {
  
  // migration matrix
  mig_mat = rcpp_to_matrix_double(args["mig_mat"]);
  
}

//------------------------------------------------
// load demography parameter values
void Parameters::load_demog_params(Rcpp::List args) {
  
  // distributions
  life_table = rcpp_to_vector_double(args["life_table"]);
  n_life_table = int(life_table.size());
  age_death = rcpp_to_vector_double(args["age_death"]);
  age_stable = rcpp_to_vector_double(args["age_stable"]);
  
}

//------------------------------------------------
// load sampling strategy parameter values
void Parameters::load_sampling_params(Rcpp::List args) {
  load_sampling_params_daily(args);
  load_sampling_params_sweep(args);
}

//------------------------------------------------
// load daily sampling strategy parameter values
void Parameters::load_sampling_params_daily(Rcpp::List args) {
  
  // return if no daily outputs
  any_daily_outputs = args["any_daily_outputs"];
  if (!any_daily_outputs) {
    n_daily_outputs = 0;
    return;
  }
  
  // extract dataframe
  Rcpp::List daily_df = args["daily"];
  
  // extract columns into vectors
  daily_deme = Rcpp::as<vector<int>>(daily_df["deme"]);
  n_daily_outputs = daily_deme.size();
  
  vector<int> daily_measure_int = Rcpp::as<vector<int>>(daily_df["measure"]);
  daily_measure = std::vector<Measure>(daily_measure_int.size());
  for (unsigned int i = 0; i < daily_measure_int.size(); ++i) {
    daily_measure[i] = static_cast<Measure>(daily_measure_int[i]);
  }
  
  vector<int> daily_state_int = Rcpp::as<vector<int>>(daily_df["state"]);
  daily_state = std::vector<Model_state>(daily_state_int.size());
  for (unsigned int i = 0; i < daily_state_int.size(); ++i) {
    daily_state[i] = static_cast<Model_state>(daily_state_int[i]);
  }
  
  vector<int> daily_diagnostic_int = Rcpp::as<vector<int>>(daily_df["diagnostic"]);
  daily_diagnostic = std::vector<Diagnostic>(daily_diagnostic_int.size());
  for (unsigned int i = 0; i < daily_diagnostic_int.size(); ++i) {
    daily_diagnostic[i] = static_cast<Diagnostic>(daily_diagnostic_int[i]);
  }
  
  daily_age_min = Rcpp::as<vector<int>>(daily_df["age_min"]);
  daily_age_max = Rcpp::as<vector<int>>(daily_df["age_max"]);
  daily_inoculations = Rcpp::as<vector<int>>(daily_df["inoculations"]);
  
  // create a map to assist in working out if a host is required to produce
  // daily output. The key to the map is a deme-&-age combination (1-year age
  // groups). The value associated with this key is a list of all output indices
  // that apply to this group.
  //
  // For example, if output field 2 needs to be the prevalence in 0-5 year olds
  // in deme 1, then the map at key {1,0} and {1,1} etc. up to {1,5} will all
  // contain the value 2 (the output index), along with any other output indices
  // that apply to this group.
  //
  // a simpler vector object also exists just for checking if a deme is required
  // in any outputs. This avoids looping through hosts in demes that are not
  // required in any output.
  
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
  
  /*
  // uncomment this block of code to print out the full daily map
  print("DAILY MAP:");
  for (auto it = daily_map.cbegin(); it != daily_map.cend(); ++it) {
    pair<int, int> map_key = it->first;
    cout << "deme" << map_key.first << " age" << map_key.second << ": ";
    for (int i = 0; i < it->second.size(); ++i) {
      cout << it->second[i] << " ";
    }
    cout << "\n";
  }
  print("DAILY DEME MAP:");
  print_vector(daily_flag_deme);
  Rcpp::stop("foobar");
  */
}

//------------------------------------------------
// load population sweep sampling strategy parameter values
void Parameters::load_sampling_params_sweep(Rcpp::List args) {
  
  // return if no sweep outputs
  any_sweep_outputs = args["any_sweep_outputs"];
  if (!any_sweep_outputs) {
    n_sweep_outputs = 0;
    return;
  }
  
  // extract dataframe
  Rcpp::List sweep_df = args["sweeps"];
  
  // extract columns into vectors
  sweep_time = Rcpp::as<vector<int>>(sweep_df["time"]);
  n_sweep_outputs = sweep_time.size();
  
  sweep_time_ordered = Rcpp::as<vector<int>>(args["sweep_time_ordered"]);;
  
  sweep_deme = Rcpp::as<vector<int>>(sweep_df["deme"]);
  
  vector<int> sweep_measure_int = Rcpp::as<vector<int>>(sweep_df["measure"]);
  sweep_measure = std::vector<Measure>(sweep_measure_int.size());
  for (unsigned int i = 0; i < sweep_measure_int.size(); ++i) {
    sweep_measure[i] = static_cast<Measure>(sweep_measure_int[i]);
  }
  
  vector<int> sweep_state_int = Rcpp::as<vector<int>>(sweep_df["state"]);
  sweep_state = std::vector<Model_state>(sweep_state_int.size());
  for (unsigned int i = 0; i < sweep_state_int.size(); ++i) {
    sweep_state[i] = static_cast<Model_state>(sweep_state_int[i]);
  }
  
  vector<int> sweep_diagnostic_int = Rcpp::as<vector<int>>(sweep_df["diagnostic"]);
  sweep_diagnostic = std::vector<Diagnostic>(sweep_diagnostic_int.size());
  for (unsigned int i = 0; i < sweep_diagnostic_int.size(); ++i) {
    sweep_diagnostic[i] = static_cast<Diagnostic>(sweep_diagnostic_int[i]);
  }
  
  sweep_age_min = Rcpp::as<vector<int>>(sweep_df["age_min"]);
  sweep_age_max = Rcpp::as<vector<int>>(sweep_df["age_max"]);
  sweep_inoculations = Rcpp::as<vector<int>>(sweep_df["inoculations"]);
  
}

//------------------------------------------------
// load run parameter values
void Parameters::load_run_params(Rcpp::List args) {
  
  // load values
  max_time = rcpp_to_int(args["max_time"]);
  output_format = rcpp_to_int(args["output_format"]);
  save_transmission_record = rcpp_to_bool(args["save_transmission_record"]);
  transmission_record_location = rcpp_to_string(args["transmission_record_location"]);
  silent = rcpp_to_bool(args["silent"]);
  pb_markdown = rcpp_to_bool(args["pb_markdown"]);
  
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

