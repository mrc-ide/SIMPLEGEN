
#include "Sweep.h"

using namespace std;

//------------------------------------------------
// initialise
void Sweep::init(Parameters &params_, cpp11::data_frame &sweep_df) {
  
  // store pointer to parameters
  params = &params_;
  
  // extract raw values from sweep_df data.frame
  vector<int> measure_raw = cpp_to_vector_int(sweep_df["measure"]);
  n_rows = int(measure_raw.size());
  measure = vector<Measure>(n_rows);
  vector<int> state_raw = cpp_to_vector_int(sweep_df["state"]);
  state = vector<Sample_state>(n_rows);
  vector<int> diagnostic_raw = cpp_to_vector_int(sweep_df["diagnostic"]);
  diagnostic = vector<Diagnostic>(n_rows);
  for (int i = 0; i < n_rows; ++i) {
    measure[i] = static_cast<Measure>(measure_raw[i]);
    state[i] = static_cast<Sample_state>(state_raw[i]);
    diagnostic[i] = static_cast<Diagnostic>(diagnostic_raw[i]);
  }
  deme = cpp_to_vector_int(sweep_df["deme"]);
  age_min = cpp_to_vector_int(sweep_df["age_min"]);
  age_max = cpp_to_vector_int(sweep_df["age_max"]);
  time = cpp_to_vector_int(sweep_df["time"]);
  
  // create key specifying which rows of sweep_df apply to each [deme][1-year
  // age group] combination
  samp_key = vector<vector<vector<int>>>(params->n_demes, vector<vector<int>>(params->max_age + 1));
  for (int k = 0; k < params->n_demes; ++k) {
    for (int i = 0; i < (params->max_age + 1); ++i) {
      for (int j = 0; j < n_rows; ++j) {
        if (((deme[j] == k) || (deme[j] == -1)) && (age_min[j] <= i) && (age_max[j] >= i)) {
          samp_key[k][i].push_back(j);
        }
      }
    }
  }
}

//------------------------------------------------
// calculate values on host and vector populations
void Sweep::calculate(vector<double> &ret, Host_pop &host_pop, vector<Mosquito_pop> &mosq_pop, int t) {
  
  // for prevalence calculations we will need a denominator to go with each
  // value in ret. This will be divided through at the end
  vector<int> denom(ret.size());
  
  // loop through all individuals and work out which rows in the original
  // sweep_df to this individual using the key
  for (size_t i = 0; i < host_pop.host_vec.size(); ++i) {
    int this_deme = host_pop.host_vec[i].deme;
    int this_age = host_pop.host_vec[i].get_age_years(t);
    for (size_t j = 0; j < samp_key[this_deme][this_age].size(); ++j) {
      int this_row = samp_key[this_deme][this_age][j];
      
      // counts and prevalence calculated the same way initially
      if ((measure[this_row] == Measure_count) || (measure[this_row] == Measure_prevalence)) {
        if (state[this_row] == Samp_S) {
          denom[this_row]++;
          if (host_pop.host_vec[i].get_host_state() == Host_Sh) {
            ret[this_row] += 1.0;
          }
        } else if (state[this_row] == Samp_E) {
          denom[this_row]++;
          if (host_pop.host_vec[i].get_host_state() == Host_Eh) {
            ret[this_row] += 1.0;
          }
        } else if (state[this_row] == Samp_P) {
          denom[this_row]++;
          if (host_pop.host_vec[i].get_host_state() == Host_Ph) {
            ret[this_row] += 1.0;
          }
        } else if (state[this_row] == Samp_H) {
          denom[this_row]++;
          ret[this_row] += 1.0;
        } else if (state[this_row] == Samp_A) {
          denom[this_row]++;
          if (host_pop.host_vec[i].get_host_state() == Host_Ah) {
            if (diagnostic[this_row] == Diagnostic_true) {
              ret[this_row] += 1.0;
            } else if (diagnostic[this_row] == Diagnostic_microscopy) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_microscopy(t);
            } else if (diagnostic[this_row] == Diagnostic_PCR) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_PCR(t);
            }
          }
        } else if (state[this_row] == Samp_C) {
          denom[this_row]++;
          if (host_pop.host_vec[i].get_host_state() == Host_Ch) {
            if (diagnostic[this_row] == Diagnostic_true) {
              ret[this_row] += 1.0;
            } else if (diagnostic[this_row] == Diagnostic_microscopy) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_microscopy(t);
            } else if (diagnostic[this_row] == Diagnostic_PCR) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_PCR(t);
            }
          }
        }
      }
      
      // incidence by passive case detection at time of treatment
      if (measure[this_row] == Measure_incidence_passive) {
        if ((host_pop.host_vec[i].get_time_treatment() == t) && (host_pop.host_vec[i].get_time_host_state_change() == t) && (host_pop.host_vec[i].get_host_state() == Host_Ah) && ((host_pop.host_vec[i].get_host_state_previous() == Host_Eh) || (host_pop.host_vec[i].get_host_state_previous() == Host_Ch))) {
          if (diagnostic[this_row] == Diagnostic_true) {
            ret[this_row] += 1.0;
          } else if (diagnostic[this_row] == Diagnostic_microscopy) {
            ret[this_row] += host_pop.host_vec[i].get_detectability_microscopy(t);
          } else if (diagnostic[this_row] == Diagnostic_PCR) {
            ret[this_row] += host_pop.host_vec[i].get_detectability_PCR(t);
          }
        }
      }
      
      // incidence by active case detection (100% of cases detected)
      if (measure[this_row] == Measure_incidence_active) {
        if (state[this_row] == Samp_A) {
          if ((host_pop.host_vec[i].get_time_host_state_change() == t) && (host_pop.host_vec[i].get_host_state() == Host_Ah) && ((host_pop.host_vec[i].get_host_state_previous() == Host_Eh) || (host_pop.host_vec[i].get_host_state_previous() == Host_Ch))) {
            if (diagnostic[this_row] == Diagnostic_true) {
              ret[this_row] += 1.0;
            } else if (diagnostic[this_row] == Diagnostic_microscopy) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_microscopy(t);
            } else if (diagnostic[this_row] == Diagnostic_PCR) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_PCR(t);
            }
          }
        } else if (state[this_row] == Samp_C) {
          if ((host_pop.host_vec[i].get_time_host_state_change() == t) && (host_pop.host_vec[i].get_host_state() == Host_Ch) && (host_pop.host_vec[i].get_host_state_previous() == Host_Eh)) {
            if (diagnostic[this_row] == Diagnostic_true) {
              ret[this_row] += 1.0;
            } else if (diagnostic[this_row] == Diagnostic_microscopy) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_microscopy(t);
            } else if (diagnostic[this_row] == Diagnostic_PCR) {
              ret[this_row] += host_pop.host_vec[i].get_detectability_PCR(t);
            }
          }
        }
      }
    }
  }
  
  // divide through by denominator for prevalence
  for (int i = 0; i < n_rows; ++i) {
    if (measure[i] == Measure_prevalence) {
      ret[i] /= double(denom[i]);
    }
  }
  
  // calculate EIR
  for (int i = 0; i < n_rows; ++i) {
    if (measure[i] == Measure_EIR) {
      
      // weighted average over all demes, or calculate in single deme
      if (deme[i] == -1) {
        int Iv_total = 0;
        int H_total = 0;
        for (int k = 0; k < params->n_demes; ++k) {
          Iv_total += mosq_pop[k].Iv;
          H_total += host_pop.H[k];
        }
        ret[i] = params->a * double(Iv_total) / double(H_total);
      } else {
        ret[i] = params->a * double(mosq_pop[deme[i]].Iv) / double(host_pop.H[deme[i]]);
      }
    }
  }
  
}

//------------------------------------------------
// print raw input data.frame
void Sweep::print_df() {
  print("measure", "state", "diagnostic", "deme", "age_min", "age_max", "time");
  for (int i = 0; i < n_rows; ++i) {
    print(measure[i], state[i], diagnostic[i], deme[i], age_min[i], age_max[i], time[i]);
  }
}

//------------------------------------------------
// print age key
void Sweep::print_samp_key() {
  print_array(samp_key);
}
