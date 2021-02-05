
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
  duration_prophylactic = rcpp_to_vector_double(args["duration_prophylactic"]);
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
  /*
  // whether samples are being obtained
  obtain_samples = rcpp_to_bool(args["obtain_samples"]);
  
  // if so, load other parameters
  if (obtain_samples) {
    ss_time = rcpp_to_vector_int(args["ss_time"]);
    ss_deme = rcpp_to_vector_int(args["ss_deme"]);
    int n_row = int(ss_time.size());
    
    vector<string> ss_case_detection_string = rcpp_to_vector_string(args["ss_case_detection"]);
    ss_case_detection = vector<Case_detection>(n_row);
    for (unsigned int i = 0; i < ss_case_detection_string.size(); ++i) {
      if (ss_case_detection_string[i] == "active") {
        ss_case_detection[i] = active;
      } else if (ss_case_detection_string[i] == "passive") {
        ss_case_detection[i] = passive;
      } else {
        Rcpp::stop("error in Parameters::load_sampling_params(): case detection method not recognised");
      }
    }
    
    vector<string> ss_diagnosis_string = rcpp_to_vector_string(args["ss_diagnosis"]);
    ss_diagnosis = vector<Diagnosis>(n_row);
    for (unsigned int i = 0; i < ss_diagnosis_string.size(); ++i) {
      if (ss_diagnosis_string[i] == "microscopy") {
        ss_diagnosis[i] = microscopy;
      } else if (ss_diagnosis_string[i] == "PCR") {
        ss_diagnosis[i] = PCR;
      } else {
        Rcpp::stop("error in Parameters::load_sampling_params(): diagnosis method not recognised");
      }
    }
    
    ss_n = rcpp_to_vector_int(args["ss_n"]);
    
  }
  */
}

//------------------------------------------------
// load run parameter values
void Parameters::load_run_params(Rcpp::List args) {
  
  // load values
  max_time = rcpp_to_int(args["max_time"]);
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
  print_vector(duration_prophylactic);
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

