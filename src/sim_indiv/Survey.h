
#pragma once

#include "Parameters.h"
#include "../misc_v17.h"
#include "Host_pop.h"
//#include "Mosquito_pop.h"

//------------------------------------------------
// enumerate sampling strategy methods
//enum Model_state {Model_S, Model_E, Model_A, Model_C, Model_P, Model_H, Model_Sv, Model_Ev, Model_Iv, Model_M};
//enum Measure {Measure_count, Measure_prevalence, Measure_incidence, Measure_EIR};
//enum Sampling {Sampling_none, Sampling_ACD, Sampling_PCD};
//enum Diagnostic {Diagnostic_true, Diagnostic_microscopy, Diagnostic_PCR};
//enum Case_detection {Case_active, Case_passive};

//------------------------------------------------
// class for sampling from human and mosquito populations in survey format
class Survey {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointer to parameters
  Parameters* params;
  
  int time_start;
  int sample_size;
  int deme;
  
  std::vector<int> host_ID;
  std::vector<std::vector<int>> inoc_ID;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Survey() {};
  
  // main methods
  void init(Parameters &params_, cpp11::data_frame &survey_df);
  void take_survey(Host_pop &host_pop);
};
