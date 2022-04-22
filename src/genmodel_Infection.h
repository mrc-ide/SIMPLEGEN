
#pragma once

#include "Sampler_v5.h"
#include "genmodel_Mosquito.h"
#include <vector>

//------------------------------------------------
// class defining infection within human host for genetic model
class genmodel_Infection {
  
public:
  
  // PUBLIC OBJECTS
  int n_haplos;
  std::vector<int> haplo_IDs;
  std::vector<double> haplo_densities;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  genmodel_Infection() {};
  
  // main events
  void denovo(int this_row, int &next_haplo_ID, double alpha,
              std::vector<int> &output_record_row,
              std::vector<int> &output_child_haplo_ID,
              std::vector<std::vector<int>> &output_parent_haplo_ID);
  int draw_haplo();
  void sample_mosquito(int this_row, Sampler &sampler_hepatocyte, int &next_haplo_ID, double alpha,
                       genmodel_Mosquito &m, std::vector<int> &parent_infection_IDs,
                       std::vector<int> &output_record_row,
                       std::vector<int> &output_child_haplo_ID,
                       std::vector<std::vector<int>> &output_parent_haplo_ID);
  void print_status();
  
};
