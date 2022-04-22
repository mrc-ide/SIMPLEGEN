
#pragma once

#include "Sampler_v5.h"
#include <vector>

// forward-declare host class
class genmodel_Host;

//------------------------------------------------
// class defining mosquito for genetic model
class genmodel_Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  int infection_row;
  int infection_ID;
  int n_oocysts;
  std::vector<std::pair<int, int>> oocysts_parents;
  std::vector<int> oocyst_counts;
  std::vector<std::vector<int>> oocysts_products;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  genmodel_Mosquito() {};
  
  // main events
  void infection(int infection_row, int infection_ID, Sampler &sampler_oocyst, genmodel_Host &host,
                 std::vector<int> &parent_infection_IDs, std::vector<double> &parent_infection_densities);
  void draw_haplo(int this_row, int n_hepatocytes, int &next_haplo_ID,
                  std::vector<int> &haplo_IDs, std::vector<double> &haplo_densities, double alpha,
                  std::vector<int> &output_record_row, std::vector<int> &output_child_haplo_ID,
                  std::vector<std::vector<int>> &output_parent_haplo_ID);
  void print_status();
  
};
