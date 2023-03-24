
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
  load_sampling_params(args);
  load_run_params(args);
  define_samplers();
  
}

//------------------------------------------------
// load epi model parameter values
void Parameters::load_model_params(cpp11::list args) {
  
  // scalar epi parameters
  a = cpp_to_double(args["a"]);
  p = cpp_to_double(args["p"]);
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
  max_age = int(age_stable.size()) - 1;
}

//------------------------------------------------
// load sampling data.frames
void Parameters::load_sampling_params(cpp11::list args) {
  
  // daily outputs
  any_daily_outputs = cpp_to_bool(args["any_daily_outputs"]);
  if (any_daily_outputs) {
    daily_df = args["daily"];
  }
  
  // sweep outputs
  any_sweep_outputs = cpp_to_bool(args["any_sweep_outputs"]);
  if (any_sweep_outputs) {
    sweeps_df = args["sweeps"];
  }
  
  // survey outputs
  any_survey_outputs = cpp_to_bool(args["any_survey_outputs"]);
  if (any_survey_outputs) {
    surveys_df = args["surveys"];
  }
}

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
// draw from stable demography distribution
int Parameters::draw_age_stable() {
  return sampler_age_stable.draw();
}

//------------------------------------------------
// draw from time to death distribution
int Parameters::draw_age_death() {
  return sampler_age_death.draw();
}

//------------------------------------------------
// draw duration of acute bloodstage infection
int Parameters::draw_duration_acute(int n) {
  return sampler_duration_acute[min(n, n_duration_acute - 1)].draw() + 1;
}

//------------------------------------------------
// draw duration of chronic bloodstage infection
int Parameters::draw_duration_chronic(int n) {
  return sampler_duration_chronic[min(n, n_duration_chronic - 1)].draw() + 1;
}

//------------------------------------------------
// draw time to treatment of acute infection
int Parameters::draw_time_treatment_acute(int n) {
  return sampler_time_treatment_acute[min(n, n_time_treatment_acute - 1)].draw() + 1;
}

//------------------------------------------------
// draw time to treatment of chronic infection
int Parameters::draw_time_treatment_chronic(int n) {
  return sampler_time_treatment_chronic[min(n, n_time_treatment_chronic - 1)].draw() + 1;
}

//------------------------------------------------
// draw duration of prophylactic period
int Parameters::draw_duration_prophylactic(int n) {
  return sampler_duration_prophylactic[min(n, n_duration_prophylactic - 1)].draw() + 1;
}

//------------------------------------------------
// return probability of becoming infected upon bite
double Parameters::get_prob_infection(int n) {
  return prob_infection[min(n, n_prob_infection - 1)];
}

//------------------------------------------------
// return probability of transitioning from liverstage to acute bloodstage infection
double Parameters::get_prob_acute(int n) {
  return prob_acute[min(n, n_prob_acute - 1)];
}

//------------------------------------------------
// return probability of transitioning from acute to chronic bloodstage infection
double Parameters::get_prob_AC(int n) {
  return prob_AC[min(n, n_prob_AC - 1)];
}

//------------------------------------------------
// return probability of detection by microscopy in acute phase at time t, given
// that entered phase at time t0
double Parameters::get_detectability_microscopy_acute(int t, int t0, int n) {
  int n_final = min(n, n_detectability_microscopy_acute - 1);
  int t_final = min(t - t0, int(detectability_microscopy_acute[n_final].size()) - 1);
  return detectability_microscopy_acute[n_final][t_final];
}

//------------------------------------------------
// return probability of detection by microscopy in chronic phase at time t,
// given that entered phase at time t0
double Parameters::get_detectability_microscopy_chronic(int t, int t0, int n) {
  int n_final = min(n, n_detectability_microscopy_chronic - 1);
  int t_final = min(t - t0, int(detectability_microscopy_chronic[n_final].size()) - 1);
  return detectability_microscopy_chronic[n_final][t_final];
}

//------------------------------------------------
// return probability of detection by PCR in acute phase at time t, given that
// entered phase at time t0
double Parameters::get_detectability_PCR_acute(int t, int t0, int n) {
  int n_final = min(n, n_detectability_PCR_acute - 1);
  int t_final = min(t - t0, int(detectability_PCR_acute[n_final].size()) - 1);
  return detectability_PCR_acute[n_final][t_final];
}

//------------------------------------------------
// return probability of detection by PCR in chronic phase at time t, given that
// entered phase at time t0
double Parameters::get_detectability_PCR_chronic(int t, int t0, int n) {
  int n_final = min(n, n_detectability_PCR_chronic - 1);
  int t_final = min(t - t0, int(detectability_PCR_chronic[n_final].size()) - 1);
  return detectability_PCR_chronic[n_final][t_final];
}

//------------------------------------------------
// return probability of onward infection in acute sexual phase at time t, given
// that entered phase at time t0
double Parameters::get_infectivity_acute(int t, int t0, int n) {
  int n_final = min(n, n_infectivity_acute - 1);
  int t_final = min(t - t0, int(infectivity_acute[n_final].size()) - 1);
  return infectivity_acute[n_final][t_final];
}

//------------------------------------------------
// return probability of onward infection in chronic sexual phase at time t,
// given that entered phase at time t0
double Parameters::get_infectivity_chronic(int t, int t0, int n) {
  int n_final = min(n, n_infectivity_chronic - 1);
  int t_final = min(t - t0, int(infectivity_chronic[n_final].size()) - 1);
  return infectivity_chronic[n_final][t_final];
}

//------------------------------------------------
// print summary of parameters
void Parameters::summary() {
  
  // print epi scalars
  print("a:", a);
  print("p:", p);
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
