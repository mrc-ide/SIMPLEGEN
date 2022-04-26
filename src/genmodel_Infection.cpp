
#include "genmodel_Infection.h"
#include "probability_v17.h"
#include "misc_v14.h"

using namespace std;


//------------------------------------------------
// create haplotypes de novo
void genmodel_Infection::denovo(int this_row, int &next_haplo_ID, double alpha,
                                vector<int> &output_record_row,
                                vector<int> &output_child_haplo_ID,
                                vector<vector<int>> &output_parent_haplo_ID) {
  
  // single haplotype with corresponding density
  n_haplos = 1;
  haplo_IDs = vector<int>(1, next_haplo_ID++);
  haplo_densities = vector<double>(1, rgamma1(alpha, 1.0));
  
  // store output
  output_record_row.push_back(this_row);
  output_child_haplo_ID.push_back(haplo_IDs[0]);
  output_parent_haplo_ID.push_back({-1});
}

//------------------------------------------------
// draw haplotype ID at random proportional to haplotype densities
int genmodel_Infection::draw_haplo() {
  int i = sample1(haplo_densities, sum(haplo_densities));
  return haplo_IDs[i];
}

//------------------------------------------------
// sample haplotypes from mosquito
void genmodel_Infection::sample_mosquito(int this_row, Sampler &sampler_hepatocyte, int &next_haplo_ID, double alpha,
                                         genmodel_Mosquito &m, vector<int> &parent_infection_IDs,
                                         vector<int> &output_record_row,
                                         vector<int> &output_child_haplo_ID,
                                         vector<vector<int>> &output_parent_haplo_ID) {
  
  // check that parent infection ID is singular (no mosquito superinfection),
  // and present in mosquito
  if (parent_infection_IDs.size() > 1) {
    Rcpp::stop("human infections (event 1) must have a single parent_infection_ID because mosquito superinfection is not allowed under this model");
  }
  if (m.infection_ID != parent_infection_IDs[0]) {
    Rcpp::stop("requesting infection ID that is not present in mosquito");
  }
  
  // draw number of hepatocytes
  int n_hepatocytes = sampler_hepatocyte.draw() + 1;
  
  // sample haplos from mosquito
  m.draw_haplo(this_row, n_hepatocytes, next_haplo_ID, haplo_IDs, haplo_densities, alpha,
               output_record_row, output_child_haplo_ID, output_parent_haplo_ID);
  
  // store number of haplos
  n_haplos = haplo_IDs.size();
}
//------------------------------------------------
// print class
void genmodel_Infection::print_status() {
  print("n_haplos:", n_haplos);
  print("haplo_IDs:");
  print_vector(haplo_IDs);
  print("haplo_densities:");
  print_vector(haplo_densities);
}
