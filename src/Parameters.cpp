
#include "Parameters.h"

using namespace std;

//------------------------------------------------
// declare static member variables

// scalar epi parameters
double Parameters::a;
double Parameters::p;
double Parameters::mu;
double Parameters::prob_AC;
int Parameters::u;
int Parameters::v;
int Parameters::g;
int Parameters::max_innoculations;

// epi distributions
vector<double> Parameters::prob_acute;
vector<double> Parameters::prob_infection;
vector<double> Parameters::infectivity_acute;
vector<double> Parameters::infectivity_chronic;
vector<vector<double>> Parameters::duration_acute;
vector<vector<double>> Parameters::duration_chronic;

// deme parameters
vector<int> Parameters::H;
vector<int> Parameters::seed_infections;
vector<int> Parameters::M;

// demog parameters
vector<double> Parameters::age_death;
vector<double> Parameters::age_stable;

// run parameters
int Parameters::max_time;
bool Parameters::output_daily_counts;
bool Parameters::output_age_distributions;
bool Parameters::output_infection_history;
bool Parameters::silent;
vector<int> Parameters::output_age_times;

//------------------------------------------------
// load epi parameter values
void Parameters::load_epi_params(double a, double p, double mu, double prob_AC,
                                 int u, int v, int g, int max_innoculations,
                                 vector<double> prob_acute,
                                 vector<double> prob_infection,
                                 vector<double> infectivity_acute,
                                 vector<double> infectivity_chronic,
                                 vector<vector<double>> duration_acute,
                                 vector<vector<double>> duration_chronic) {
  
  // define scalars
  this->a = a;
  this->p = p;
  this->mu = mu;
  this->prob_AC = prob_AC;
  this->u = u;
  this->v = v;
  this->g = g;
  this->max_innoculations = max_innoculations;
  
  // distributions
  this->prob_acute = prob_acute;
  this->prob_infection = prob_infection;
  this->infectivity_acute = infectivity_acute;
  this->infectivity_chronic = infectivity_chronic;
  this->duration_acute = duration_acute;
  this->duration_chronic = duration_chronic;
}

//------------------------------------------------
// load deme parameter values
void Parameters::load_deme_params(vector<int> H,
                                  vector<int> seed_infections,
                                  vector<int> M) {
  
  // distributions
  this->H = H;
  this->seed_infections = seed_infections;
  this->M = M;
}

//------------------------------------------------
// load demography parameter values
void Parameters::load_demog_params(vector<double> age_death,
                                   vector<double> age_stable) {
  
  // distributions
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
  print("prob_AC:", prob_AC);
  print("u:", u);
  print("v:", v);
  print("g:", g);
  print("max_innoculations:", max_innoculations);
  
  // print epi distributions
  print("prob_acute:");
  print_vector(prob_acute);
  print("prob_infection:");
  print_vector(prob_infection);
  print("infectivity_acute:");
  print_vector(infectivity_acute);
  print("infectivity_chronic:");
  print_vector(infectivity_chronic);
  print("duration_acute:");
  print_matrix(duration_acute);
  print("duration_chronic:");
  print_matrix(duration_chronic);
  
  // print deme parameters
  print("H:");
  print_vector(H);
  print("seed_infections:");
  print_vector(seed_infections);
  print("M:");
  print_vector(M);
  
  // print demog parameters
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

