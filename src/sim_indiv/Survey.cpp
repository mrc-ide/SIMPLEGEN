
#include "Survey.h"

using namespace std;

//------------------------------------------------
// initialise
void Survey::init(Parameters &params_, cpp11::data_frame &survey_df) {
  
  // store pointer to parameters
  params = &params_;
  
  vector<int> time_start_vec = cpp_to_vector_int(survey_df["time_start"]);
  time_start = time_start_vec[0];
  vector<int> sample_size_vec = cpp_to_vector_int(survey_df["sample_size"]);
  sample_size = sample_size_vec[0];
  vector<int> deme_vec = cpp_to_vector_int(survey_df["deme"]);
  deme = deme_vec[0];
  
}

//------------------------------------------------
// take sample from population
void Survey::take_survey(Host_pop &host_pop) {
  
  // check sample size not too large
  if ((deme == -1) && (sample_size > host_pop.host_vec.size())) {
    cpp11::stop("sample size exceeds population size");
  } else if (sample_size > host_pop.H[deme]) {
    cpp11::stop("sample size exceeds population size");
  }
  
  // sample at random from population
  vector<int> samp_index;
  if (deme == -1) {
    samp_index = sample4(sample_size, 0, host_pop.host_vec.size() - 1);
  } else {
    samp_index = sample4(sample_size, 0, host_pop.H[deme] - 1);
    for (int i = 0; i < sample_size; ++i) {
      samp_index[i] = host_pop.host_index[deme][samp_index[i]];
    }
  }
  
  // get IDs
  host_ID = vector<int>(sample_size);
  inoc_ID = vector<vector<int>>(sample_size);
  for (int i = 0; i < sample_size; ++i) {
    host_ID[i] = host_pop.host_vec[samp_index[i]].host_ID;
    for (int j = 0; j < params->max_inoculations; ++j) {
      if (host_pop.host_vec[samp_index[i]].inoc[j].active) {
        inoc_ID[i].push_back(host_pop.host_vec[samp_index[i]].inoc[j].ID);
      }
    }
  }
  
}
