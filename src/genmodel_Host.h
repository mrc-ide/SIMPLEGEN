
#pragma once

#include "Sampler_v5.h"
#include "genmodel_Infection.h"
#include <vector>
#include <map>

// forward-declare mosquito class
class genmodel_Mosquito;

//------------------------------------------------
// class defining human host for genetic model
class genmodel_Host {
  
public:
  
  // PUBLIC OBJECTS
  int n_infections;
  std::map<int, genmodel_Infection> infections_map;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  genmodel_Host();
  
  // main events
  void denovo_infection(int this_row, int infection_ID, int &next_haplo_ID, double alpha,
                        std::vector<int> &output_record_row,
                        std::vector<int> &output_child_haplo_ID,
                        std::vector<std::vector<int>> &output_parent_haplo_ID);
  void infection(int this_row, int infection_ID, Sampler &sampler_hepatocyte, int &next_haplo_ID, double alpha,
                 genmodel_Mosquito &m, std::vector<int> &parent_infection_IDs,
                 std::vector<int> &output_row,
                 std::vector<int> &output_child_haplo_ID,
                 std::vector<std::vector<int>> &output_parent_haplo_ID);
  int draw_haplo(const std::vector<int> &infection_IDs, const std::vector<double> &infection_densities);
  void print_status();
  
};
