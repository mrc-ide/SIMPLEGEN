
#include "main.h"
#include "Parameters.h"
#include "Dispatcher.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
Rcpp::List indiv_sim_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define parameters object and load values
  Parameters params;
  params.load_epi_params(rcpp_to_double(args["a"]),
                         rcpp_to_double(args["p"]),
                         rcpp_to_double(args["mu"]),
                         rcpp_to_int(args["u"]),
                         rcpp_to_int(args["v"]),
                         rcpp_to_int(args["g"]),
                         rcpp_to_vector_double(args["prob_infection"]),
                         rcpp_to_vector_double(args["prob_acute"]),
                         rcpp_to_vector_double(args["prob_AC"]),
                         rcpp_to_matrix_double(args["duration_acute"]),
                         rcpp_to_matrix_double(args["duration_chronic"]),
                         rcpp_to_matrix_double(args["detectability_microscopy_acute"]),
                         rcpp_to_matrix_double(args["detectability_microscopy_chronic"]),
                         rcpp_to_matrix_double(args["detectability_PCR_acute"]),
                         rcpp_to_matrix_double(args["detectability_PCR_chronic"]),
                         rcpp_to_matrix_double(args["time_treatment_acute"]),
                         rcpp_to_matrix_double(args["time_treatment_chronic"]),
                         rcpp_to_double(args["treatment_seeking_alpha"]),
                         rcpp_to_double(args["treatment_seeking_beta"]),
                         rcpp_to_vector_double(args["duration_prophylactic"]),
                         rcpp_to_matrix_double(args["infectivity_acute"]),
                         rcpp_to_matrix_double(args["infectivity_chronic"]),
                         rcpp_to_int(args["max_inoculations"]));
  
  params.load_deme_params(rcpp_to_vector_int(args["H"]),
                          rcpp_to_vector_int(args["seed_infections"]),
                          rcpp_to_vector_int(args["M"]));
  
  params.load_demog_params(rcpp_to_vector_double(args["life_table"]),
                           rcpp_to_vector_double(args["age_death"]),
                           rcpp_to_vector_double(args["age_stable"]));
  
  params.load_run_params(rcpp_to_int(args["max_time"]),
                         rcpp_to_bool(args["save_transmission_record"]),
                         rcpp_to_string(args["transmission_record_location"]),
                         rcpp_to_bool(args["output_daily_counts"]),
                         rcpp_to_bool(args["output_age_distributions"]),
                         rcpp_to_vector_int(args["output_age_times"]),
                         rcpp_to_bool(args["silent"]));
  
  //params.summary();
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.run_simulation();
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("daily_values") = dispatcher.daily_values);
}
#else
int main(int argc, const char * argv[]) {  // main function when not using Rcpp
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define paths to arg files
  string file_path = "/Users/rverity/Dropbox/Bob/Work/My Programs/Simulation/SIMPLEGEN/R_ignore/SIMPLEGEN_Xcode/args/";
  
  string epi_scalar_path = file_path + "epi_scalar.txt";
  string prob_acute_path = file_path + "prob_acute.txt";
  string prob_infection_path = file_path + "prob_infection.txt";
  string infectivity_acute_path = file_path + "infectivity_acute.txt";
  string infectivity_chronic_path = file_path + "infectivity_chronic.txt";
  string duration_acute_path = file_path + "duration_acute.txt";
  string duration_chronic_path = file_path + "duration_chronic.txt";
  
  string H_path = file_path + "H.txt";
  string seed_infections_path = file_path + "seed_infections.txt";
  string M_path = file_path + "M.txt";
  
  string life_table_path = file_path + "life_table.txt";
  string age_death_path = file_path + "age_death.txt";
  string age_stable_path = file_path + "age_stable.txt";
  
  string run_scalar_path = file_path + "run_scalar.txt";
  string output_age_times_path = file_path + "output_age_times.txt";
  
  // read in scalar arguments from file
  vector<double> epi_scalar = file_to_vector_double(epi_scalar_path);
  vector<double> run_scalar = file_to_vector_double(run_scalar_path);
  int max_time = int(run_scalar[0]);
  bool output_daily_counts = int(run_scalar[1]);
  bool output_age_distributions = int(run_scalar[2]);
  bool output_infection_history = int(run_scalar[3]);
  bool silent = int(run_scalar[4]);
  
  // define parameters object and load values
  Parameters params;
  params.load_epi_params(epi_scalar[0], epi_scalar[1], epi_scalar[2], epi_scalar[3],
                         epi_scalar[4], epi_scalar[5], epi_scalar[6], epi_scalar[7],
                         file_to_vector_double(prob_acute_path),
                         file_to_vector_double(prob_infection_path),
                         file_to_matrix_double(infectivity_acute_path),
                         file_to_matrix_double(infectivity_chronic_path),
                         file_to_matrix_double(duration_acute_path),
                         file_to_matrix_double(duration_chronic_path));
  
  params.load_deme_params(file_to_vector_int(H_path),
                          file_to_vector_int(seed_infections_path),
                          file_to_vector_int(M_path));
  
  params.load_demog_params(file_to_vector_double(life_table_path),
                           file_to_vector_double(age_death_path),
                           file_to_vector_double(age_stable_path));
  
  params.load_run_params(max_time, output_daily_counts, output_age_distributions,
                         output_infection_history, silent,
                         file_to_vector_int(output_age_times_path));
  
  
  //params.summary();
  
  // create dispatcher object and run simulations
  Dispatcher dispatcher;
  dispatcher.run_simulation();
  
  // end timer
  chrono_timer(t1);
  
  return 0;
}
#endif
