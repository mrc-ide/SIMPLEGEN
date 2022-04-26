
#include "genmodel_Host.h"
#include "genmodel_Mosquito.h"
#include "probability_v17.h"
#include "misc_v14.h"

using namespace std;

//------------------------------------------------
// constructor
genmodel_Host::genmodel_Host() {
  n_infections = 0;
}

//------------------------------------------------
// create haplotypes de novo
void genmodel_Host::denovo_infection(int this_row, int infection_ID, int &next_haplo_ID, double alpha, vector<int> &output_record_row,
                                     vector<int> &output_child_haplo_ID, vector<vector<int>> &output_parent_haplo_ID) {
  
  // single infection with haplotypes created de novo
  n_infections = 1;
  infections_map[infection_ID] = genmodel_Infection();
  infections_map[infection_ID].denovo(this_row, next_haplo_ID, alpha, output_record_row,
                                      output_child_haplo_ID, output_parent_haplo_ID);
}

//------------------------------------------------
// draw haplotypes from infective mosquito
void genmodel_Host::infection(int this_row, int infection_ID, Sampler &sampler_hepatocyte, int &next_haplo_ID, double alpha,
                              genmodel_Mosquito &m, vector<int> &parent_infection_IDs,
                              vector<int> &output_row,
                              vector<int> &output_child_haplo_ID,
                              vector<vector<int>> &output_parent_haplo_ID) {
  
  // new infection from mosquito
  n_infections++;
  infections_map[infection_ID] = genmodel_Infection();
  infections_map[infection_ID].sample_mosquito(this_row, sampler_hepatocyte, next_haplo_ID, alpha, m, parent_infection_IDs,
                                               output_row, output_child_haplo_ID, output_parent_haplo_ID);
  
}

//------------------------------------------------
// draw a single haplotype ID
int genmodel_Host::draw_haplo(const vector<int> &infection_IDs, const vector<double> &infection_densities) {
  
  // check that all infection IDs present in map
  for (int i = 0; i < infection_IDs.size(); ++i) {
    if (infections_map.count(infection_IDs[i]) == 0) {
      Rcpp::stop("requested infection IDs do not match those in host object");
    }
  }
  
  // draw an infection proportionally to infection densities
  int i = sample1(infection_densities, sum(infection_densities));
  int this_infection_ID = infection_IDs[i];
  
  // draw a haplotype from this infection and return
  return infections_map[this_infection_ID].draw_haplo();
}

//------------------------------------------------
// print class
void genmodel_Host::print_status() {
  print("\nHost status:");
  print("n_infections:", n_infections);
  for (auto &x : infections_map) {
    print("infection_ID:", x.first);
    x.second.print_status();
  }
}
