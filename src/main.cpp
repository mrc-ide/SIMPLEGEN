
#include "main.h"
#include "misc_v6.h"
#include "Parameters.h"
#include "Dispatcher.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // split args into sublists
  Rcpp::List args_epi_parameters = args["epi_parameters"];
  
  // define parameters object and load values
  Parameters params;
  params.load_values(rcpp_to_double(args_epi_parameters["a"]),
                     rcpp_to_double(args_epi_parameters["p"]),
                     rcpp_to_double(args_epi_parameters["mu"]),
                     rcpp_to_double(args_epi_parameters["prob_AC"]),
                     rcpp_to_int(args_epi_parameters["u"]),
                     rcpp_to_int(args_epi_parameters["v"]),
                     rcpp_to_int(args_epi_parameters["g"]),
                     rcpp_to_int(args_epi_parameters["max_innoculations"]),
                     rcpp_to_vector_double(args_epi_parameters["prob_acute"]),
                     rcpp_to_vector_double(args_epi_parameters["prob_infection"]),
                     rcpp_to_vector_double(args_epi_parameters["infectivity_acute"]),
                     rcpp_to_vector_double(args_epi_parameters["infectivity_chronic"]),
                     rcpp_to_matrix_double(args_epi_parameters["duration_acute"]),
                     rcpp_to_matrix_double(args_epi_parameters["duration_chronic"]));
  //params.summary();
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.run_simulation();
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("foo") = -9);
}
#else
int main(int argc, const char * argv[]) {  // main function when not using Rcpp
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define paths to arg files
  string file_path = "/Users/rverity/Dropbox/Bob/Work/My Programs/Simulation/SIMPLEGEN/R_ignore/SIMPLEGEN_Xcode/args/";
  string scalar_path = file_path + "scalar.txt";
  string prob_acute_path = file_path + "prob_acute.txt";
  string prob_infection_path = file_path + "prob_infection.txt";
  string infectivity_acute_path = file_path + "infectivity_acute.txt";
  string infectivity_chronic_path = file_path + "infectivity_chronic.txt";
  string duration_acute_path = file_path + "duration_acute.txt";
  string duration_chronic_path = file_path + "duration_chronic.txt";
  
  // read in scalar arguments from file
  vector<double> scalar = file_to_vector_double(scalar_path);
  
  // define parameters object and load values
  Parameters params;
  params.load_values(scalar[0], scalar[1], scalar[2], scalar[3],
                     scalar[4], scalar[5], scalar[6], scalar[7],
                     file_to_vector_double(prob_acute_path),
                     file_to_vector_double(prob_infection_path),
                     file_to_vector_double(infectivity_acute_path),
                     file_to_vector_double(infectivity_chronic_path),
                     file_to_matrix_double(duration_acute_path),
                     file_to_matrix_double(duration_chronic_path));
  //params.summary();
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.run_simulation();
  
  // end timer
  chrono_timer(t1);
  
  return 0;
}
#endif
