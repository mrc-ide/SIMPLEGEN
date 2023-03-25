
#pragma once

#include "Parameters.h"
#include "Host.h"
#include "Mosquito.h"
#include "../misc_v17.h"
#include "Host_pop.h"
#include "Mosquito_pop.h"

//------------------------------------------------
// enumerate sampling elements
enum Measure {Measure_count, Measure_prevalence, Measure_incidence_active, Measure_incidence_passive, Measure_EIR};
enum Sample_state {Samp_S, Samp_E, Samp_A, Samp_C, Samp_P, Samp_H, Samp_Sv, Samp_Ev, Samp_Iv, Samp_M};
enum Diagnostic {Diagnostic_true, Diagnostic_microscopy, Diagnostic_PCR};

//------------------------------------------------
// class for sampling from human and mosquito populations in either a "daily" or
// a "sweep" format
class Sweep {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters* params;
  
  // extracted columns
  int n_rows;
  std::vector<Measure> measure;
  std::vector<Sample_state> state;
  std::vector<Diagnostic> diagnostic;
  std::vector<int> deme;
  std::vector<int> age_min;
  std::vector<int> age_max;
  std::vector<int> time;
  
  // key specifying which rows of sweep_df apply to each [deme][1-year age
  // group] combination
  std::vector<std::vector<std::vector<int>>> samp_key;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Sweep() {};
  
  // main methods
  void init(Parameters &params_, cpp11::data_frame sweep_df);
  void calculate(std::vector<double> &ret, Host_pop &host_pop, std::vector<Mosquito_pop> &mosq_pop, int t);
  void print_df();
  void print_samp_key();
};
