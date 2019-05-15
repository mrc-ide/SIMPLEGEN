
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
vector<double> Parameters::life_table;
vector<double> Parameters::age_death;
vector<double> Parameters::age_stable;

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
void Parameters::load_demog_params(vector<double> life_table,
                       vector<double> age_death,
                       vector<double> age_stable) {
  
  // distributions
  this->life_table = life_table;
  this->age_death = age_death;
  this->age_stable = age_stable;
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
  print("life_table:");
  print_vector(life_table);
  print("age_death:");
  print_vector(age_death);
  print("age_stable:");
  print_vector(age_stable);
}

