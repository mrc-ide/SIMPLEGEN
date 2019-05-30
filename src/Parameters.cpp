
#include "Parameters.h"

using namespace std;

//------------------------------------------------
// declare static member variables

// scalar epi parameters
double Parameters::a;
double Parameters::p;
double Parameters::mu;
int Parameters::u;
int Parameters::v;
int Parameters::g;
int Parameters::max_inoculations;

// epi distributions
vector<double> Parameters::prob_infection;
int Parameters::n_prob_infection;
double Parameters::max_prob_infection;
vector<double> Parameters::prob_acute;
int Parameters::n_prob_acute;
vector<double> Parameters::prob_AC;
int Parameters::n_prob_AC;
vector<vector<double>> Parameters::duration_acute;
int Parameters::n_duration_acute;
vector<vector<double>> Parameters::duration_chronic;
int Parameters::n_duration_chronic;
vector<vector<double>> Parameters::time_treatment_acute;
int Parameters::n_time_treatment_acute;
vector<vector<double>> Parameters::time_treatment_chronic;
int Parameters::n_time_treatment_chronic;
vector<vector<double>> Parameters::infectivity_acute;
int Parameters::n_infectivity_acute;
vector<vector<double>> Parameters::infectivity_chronic;
int Parameters::n_infectivity_chronic;
double Parameters::max_infectivity;

// deme parameters
vector<int> Parameters::H_init;
vector<int> Parameters::seed_infections;
vector<int> Parameters::M;
int Parameters::n_demes;

// demog parameters
vector<double> Parameters::life_table;
vector<double> Parameters::age_death;
vector<double> Parameters::age_stable;

// run parameters
int Parameters::max_time;
bool Parameters::output_daily_counts;
bool Parameters::output_age_distributions;
bool Parameters::output_infection_history;
bool Parameters::silent;
vector<int> Parameters::output_age_times;

// misc parameters
double Parameters::prob_mosq_death;

//------------------------------------------------
// load epi parameter values
void Parameters::load_epi_params(double a, double p, double mu,
                                 int u, int v, int g,
                                 vector<double> prob_infection,
                                 vector<double> prob_acute,
                                 vector<double> prob_AC,
                                 vector<vector<double>> duration_acute,
                                 vector<vector<double>> duration_chronic,
                                 vector<vector<double>> time_treatment_acute,
                                 vector<vector<double>> time_treatment_chronic,
                                 vector<vector<double>> infectivity_acute,
                                 vector<vector<double>> infectivity_chronic,
                                 int max_inoculations) {
  
  // define scalars
  this->a = a;
  this->p = p;
  this->mu = mu;
  this->u = u;
  this->v = v;
  this->g = g;
  this->max_inoculations = max_inoculations;
  
  // distributions
  this->prob_infection = prob_infection;
  n_prob_infection = int(prob_infection.size());
  this->prob_acute = prob_acute;
  n_prob_acute = int(prob_acute.size());
  this->prob_AC = prob_AC;
  n_prob_AC = int(prob_AC.size());
  this->duration_acute = duration_acute;
  n_duration_acute = int(duration_acute.size());
  this->duration_chronic = duration_chronic;
  n_duration_chronic = int(duration_chronic.size());
  this->time_treatment_acute = time_treatment_acute;
  n_time_treatment_acute = int(time_treatment_acute.size());
  this->time_treatment_chronic = time_treatment_chronic;
  n_time_treatment_chronic = int(time_treatment_chronic.size());
  this->infectivity_acute = infectivity_acute;
  n_infectivity_acute = int(infectivity_acute.size());
  this->infectivity_chronic = infectivity_chronic;
  n_infectivity_chronic = int(infectivity_chronic.size());
  
  // get max prob_infection
  max_prob_infection = max(prob_infection);
  
  // get max infectivity over all distributions
  max_infectivity = 0.0;
  for (int i = 0; i < int(infectivity_acute.size()); ++i) {
    max_infectivity = (max(infectivity_acute[i]) > max_infectivity) ? max(infectivity_acute[i]) : max_infectivity;
  }
  for (int i = 0; i < int(infectivity_chronic.size()); ++i) {
    max_infectivity = (max(infectivity_chronic[i]) > max_infectivity) ? max(infectivity_chronic[i]) : max_infectivity;
  }
  
  // misc parameters
  prob_mosq_death = 1 - exp(-mu);  // daily probability of mosquito death
  
}

//------------------------------------------------
// load deme parameter values
void Parameters::load_deme_params(vector<int> H_init,
                                  vector<int> seed_infections,
                                  vector<int> M) {
  
  // distributions
  this->H_init = H_init;
  this->seed_infections = seed_infections;
  this->M = M;
  n_demes = int(H_init.size());
}

//------------------------------------------------
// load demography parameter values
void Parameters::load_demog_params(vector<double> life_table,
                                   vector<double> age_death,
                                   vector<double> age_stable) {
  
  // distributions
  this->life_table = life_table;
  this->age_death = age_death;
  this->age_stable = age_stable;
}

//------------------------------------------------
// load run parameter values
void Parameters::load_run_params(int max_time,
                                 bool output_daily_counts,
                                 bool output_age_distributions,
                                 bool output_infection_history,
                                 bool silent,
                                 std::vector<int> output_age_times) {
  
  // load values
  this->max_time = max_time;
  this->output_daily_counts = output_daily_counts;
  this->output_age_distributions = output_age_distributions;
  this->output_infection_history = output_infection_history;
  this->silent = silent;
  this->output_age_times = output_age_times;
  
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
  print("output_daily_counts:", output_daily_counts);
  print("output_age_distributions:", output_age_distributions);
  print("output_infection_history:", output_infection_history);
  print("silent:", silent);
  
  // print run vectors
  print("output_age_times:");
  print_vector(output_age_times);
}

