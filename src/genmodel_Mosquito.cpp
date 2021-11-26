
#include "genmodel_Mosquito.h"
#include "genmodel_Host.h"
#include "misc_v11.h"
#include "probability_v14.h"

using namespace std;

//------------------------------------------------
// pass haplos from host to mosquito
void genmodel_Mosquito::infection(int infection_row, int infection_ID,
                                  Sampler &sampler_oocyst,
                                  genmodel_Host &host,
                                  vector<int> &parent_infection_IDs,
                                  vector<double> &parent_infection_densities) {
  
  // store infection time and ID
  this->infection_row = infection_row;
  this->infection_ID = infection_ID;
  
  // draw number of oocysts
  n_oocysts = sampler_oocyst.draw() + 1;
  
  // make a map in which the two parental haplo IDs are the key, and the number
  // of times this combination is observed is the value
  map<pair<int, int>, int> oocysts_parents_map;
  
  // sample parental gametocytes from host in proportion to infection densities
  // followed by haplotype densities
  for (int i = 0; i < n_oocysts; ++i) {
    int parent1 = host.draw_haplo(parent_infection_IDs, parent_infection_densities);
    int parent2 = host.draw_haplo(parent_infection_IDs, parent_infection_densities);
    
    // avoid redundant ordering of parents
    if (parent1 > parent2) {
      int tmp = parent1;
      parent1 = parent2;
      parent2 = tmp;
    }
    
    // add to map or increment count if already present
    if (oocysts_parents_map.count(make_pair(parent1, parent2)) == 0) {
      oocysts_parents_map[make_pair(parent1, parent2)] = 1;
    } else {
      oocysts_parents_map[make_pair(parent1, parent2)]++;
    }
  }
  
  // change format from map to vectors to allow efficient sampling
  int n_map = oocysts_parents_map.size();
  oocysts_parents = vector<pair<int, int>>(n_map);
  oocyst_counts = vector<int>(n_map);
  int i = 0;
  for (auto const& x : oocysts_parents_map) {
    oocysts_parents[i].first = x.first.first;
    oocysts_parents[i].second = x.first.second;
    oocyst_counts[i] = x.second;
    i++;
  }
  
  // sanity check
  if (sum(oocyst_counts) != n_oocysts) {
    Rcpp::stop("sum of oocyst_counts does not equal n_oocysts");
  }
  
  // initialise object for storing recombinant products. These will be given
  // unique values the first time they are accessed. -1 indicates elements that
  // have not yet been accessed.
  oocysts_products = vector<vector<int>>(n_map, vector<int>(4, -1));
  
}

//------------------------------------------------
// draw haplo ID
void genmodel_Mosquito::draw_haplo(int this_row, int n_hepatocytes, int &next_haplo_ID,
                                   vector<int> &haplo_IDs, vector<double> &haplo_densities, double alpha,
                                   vector<int> &output_record_row, vector<int> &output_child_haplo_ID,
                                   vector<vector<int>> &output_parent_haplo_ID) {
  
  // determine number of times each recombinant product is drawn
  map<pair<int, int>, int> draws_map;
  for (int k = 0; k < n_hepatocytes; ++k) {
    
    // draw oocyst
    int i = sample1(oocyst_counts, n_oocysts);
    
    // if parents are clonal then return first recombinant product, otherwise
    // return one of the four recombinant products with equal probability
    int j = 0;
    if (oocysts_parents[i].first != oocysts_parents[i].second) {
      j = sample2(0, 3);
    }
    
    // add to map or increment count if already present
    if (draws_map.count(make_pair(i, j)) == 0) {
      draws_map[make_pair(i, j)] = 1;
    } else {
      draws_map[make_pair(i, j)]++;
    }
    
  }
  
  // get final haplos
  int n_haplos = draws_map.size();
  haplo_IDs = vector<int>(n_haplos);
  haplo_densities = vector<double>(n_haplos);
  int k = 0;
  for (auto const& x : draws_map) {
    
    // extract values from map
    int i = x.first.first;
    int j = x.first.second;
    int n = x.second;
    
    // if the chosen element has not been accessed then fill in with new
    // haplotype ID and record event
    if (oocysts_products[i][j] == -1) {
      oocysts_products[i][j] = next_haplo_ID++;
      
      // store event
      output_record_row.push_back(infection_row);
      output_child_haplo_ID.push_back(oocysts_products[i][j]);
      output_parent_haplo_ID.push_back({oocysts_parents[i].first, oocysts_parents[i].second});
    }
    int this_product = oocysts_products[i][j];
    
    // draw corresponding density, taking into account number of times this
    // haplo is sampled. Makes use of the fact that the sum of gamma random
    // variables with same rate is also gamma
    double this_density = rgamma1(n*alpha, 1.0);
    
    // save haplo ID and density
    haplo_IDs[k] = next_haplo_ID++;
    haplo_densities[k] = this_density;
    
    // store event
    output_record_row.push_back(this_row);
    output_child_haplo_ID.push_back(haplo_IDs[k]);
    output_parent_haplo_ID.push_back({this_product});
    
    k++;
  }
  
}


//------------------------------------------------
// print class
void genmodel_Mosquito::print_status() {
  print("\nMosquito status:");
  print("infection_row:", infection_row);
  print("infection_ID:", infection_ID);
  print("n_oocysts:", n_oocysts);
  print("oocyst parents:");
  for (int i = 0; i < oocysts_parents.size(); ++i) {
    print(oocysts_parents[i].first, oocysts_parents[i].second);
  }
  print("oocyst products:");
  print_matrix(oocysts_products);
}
