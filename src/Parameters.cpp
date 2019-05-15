
#include "Parameters.h"

using namespace std;

//------------------------------------------------
// load parameter values
void Parameters::load_values(double a, double p, double mu, double prob_AC,
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
// print summary of parameters
void Parameters::summary() {
  
  // print scalars
  print("a:", a);
  print("p:", p);
  print("mu:", mu);
  print("prob_AC:", prob_AC);
  print("u:", u);
  print("v:", v);
  print("g:", g);
  print("max_innoculations:", max_innoculations);
  
  // print distributions
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
  
}

