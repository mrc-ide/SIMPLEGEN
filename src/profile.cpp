
//------------------------------------------------
// version of main function used in profiling
#ifdef RCPP_ACTIVE
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
