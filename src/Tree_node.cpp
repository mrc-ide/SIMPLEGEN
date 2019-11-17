
#include "Tree_node.h"
#include "misc_v9.h"
#include "probability_v9.h"

using namespace std;


//------------------------------------------------
// constructor calling init()
Tree_node::Tree_node(int t, const std::vector<int> &contig_lengths, Sampler &sampler_oocyst,
                     Sampler &sampler_hepatocyte, map<int, Tree_node> &tree) {
  init(t, contig_lengths, sampler_oocyst, sampler_hepatocyte, tree);
}

//------------------------------------------------
// initialise with same haplo_ID covering every chromosome
void Tree_node::init(int t, const vector<int> &contig_lengths, Sampler &sampler_oocyst,
                     Sampler &sampler_hepatocyte, map<int, Tree_node> &tree) {
  
  // store basic properties
  this->t = t;
  this->contig_lengths = contig_lengths;
  n_contigs = int(contig_lengths.size());
  
  // store pointers
  sampler_oocyst_ptr = &sampler_oocyst;
  sampler_hepatocyte_ptr = &sampler_hepatocyte;
  tree_ptr = &tree;
  
}

//------------------------------------------------
// initialise with same haplo_ID covering every chromosome
void Tree_node::draw_haplotypes_denovo(int &haplo_ID) {
  
  // initialise with a single haplotype
  n_haplotypes = 1;
  haplo_ID_vec = vector<int>(1, haplo_ID++);
  haplo_density = vector<double>(1);
  haplo_density[0] = rgamma1(1.0, 1.0);
  
  // initialise intervals with the same single haplo_ID over every chromosome
  intervals = vector<vector<vector<vector<int>>>>(1, vector<vector<vector<int>>>(n_contigs));
  for (int i = 0; i < n_contigs; ++i) {
    intervals[0][i].push_back({0, contig_lengths[i], -1});
  }
  
}

//------------------------------------------------
// initialise from ancestral inoculation IDs
//void Tree_node::draw_haplotypes_recombine(int &haplo_ID, const std::vector<int> &inoc_IDs,
//                                          double r, double alpha) {
void Tree_node::draw_haplotypes_recombine(int &haplo_ID, const std::vector<int> &inoc_IDs,
                                          double r, double alpha, map<int, Tree_node> &tree) {
  
  // get the complete vector of haplotype IDs and densities by concatenating
  // over the ancestral inoculations
  vector<int> parental_haplo_IDs;
  vector<double> parental_haplo_densities;
  for (int i = 2; i < int(inoc_IDs.size()); ++i) {  // start at i=2 because first two values represent the key ID and the time
    //push_back_multiple(parental_haplo_IDs, (*tree_ptr)[inoc_IDs[i]].haplo_ID_vec);
    //push_back_multiple(parental_haplo_densities, (*tree_ptr)[inoc_IDs[i]].haplo_density);
    
    push_back_multiple(parental_haplo_IDs, tree[inoc_IDs[i]].haplo_ID_vec);
    push_back_multiple(parental_haplo_densities, tree[inoc_IDs[i]].haplo_density);
  }
  double sum_parental_haplo_densities = sum(parental_haplo_densities);
  
  // draw number of oocysts
  int n_oocysts = sampler_oocyst_ptr->draw() + 1;
  
  // draw the parental pairs that make up each oocyst. Each parent is a haplo_ID
  // taken from one of the ancestral inoculations. haplo_IDs are drawn at random
  // with probaiblities given by the relative densities of each haplotype.
  vector<pair<int, int>> oocyst_parents(n_oocysts);
  for (int i = 0; i < n_oocysts; ++i) {
    int samp1 = sample1(parental_haplo_densities, sum_parental_haplo_densities) - 1;
    int samp2 = sample1(parental_haplo_densities, sum_parental_haplo_densities) - 1;
    oocyst_parents[i] = make_pair(parental_haplo_IDs[samp1], parental_haplo_IDs[samp2]);
  }
  
  // draw number of infected hepatocytes
  int n_hepatocytes = sampler_hepatocyte_ptr->draw() + 1;
  
  // draw the number of infected hepatocytes that originate from each oocysts
  // product
  vector<int> n_products(4*n_oocysts);
  int N = n_hepatocytes;
  int n_allocated = 0;
  for (int i = 0; i < 4*n_oocysts; ++i) {
    n_products[i] = rbinom1(N - n_allocated, 1/double(4*n_oocysts - i));
    n_allocated += n_products[i];
  }
  
  // group clonal products together
  for (int i = 0; i < n_oocysts; ++i) {
    if (oocyst_parents[i].first == oocyst_parents[i].second) {
      for (int j = 1; j < 4; ++j) {
        n_products[4*i] += n_products[4*i+j];
        n_products[4*i+j] = 0;
      }
    }
  }
  
  // for each oocysts product that leads to infected hepatocytes, create a new
  // interval map by recombination and draw the density of this haplotype. Then
  // push to invervals object
  n_haplotypes = 0;
  for (int i = 0; i < 4*n_oocysts; ++i) {
    if (n_products[i] != 0) {
      
      // update ID and density vector
      n_haplotypes++;
      haplo_ID_vec.push_back(haplo_ID++);
      haplo_density.push_back(rgamma1(alpha*n_products[i], 1.0));
      
      // draw recombinant products and push to intervals object
      draw_intervals(oocyst_parents[i % n_oocysts].first, oocyst_parents[i % n_oocysts].second, r);
    }
  }
  
}

//------------------------------------------------
// create intervals by recombination of two parental IDs, and push to intervals
// array
void Tree_node::draw_intervals(int parent0, int parent1, double r) {
  return;
  // loop through chromosomes
  vector<vector<vector<int>>> x(n_contigs);
  for (int i = 0; i < n_contigs; ++i) {
    
    // deal with clonal case
    if (parent0 == parent1) {
      x[i].push_back({0, contig_lengths[i], parent0});
      continue;
    }
    
    // draw starting parent with equal probability
    int current_parent = rbernoulli1(0.5) ? parent1 : parent0;
    
    // draw number of breakpoints. Note, the lines that follow work even if
    // n_breaks is zero
    int n_breaks = rpois1(r*contig_lengths[i]);
    
    // draw location of breakpoints ensuring no duplicates
    vector<int> break_pos(n_breaks);
    for (int j = 0; j < n_breaks; ++j) {
      break_pos[j] = sample2(1, contig_lengths[i]);
    }
    sort(break_pos.begin(), break_pos.end());
    remove_duplicates(break_pos);
    
    // create intervals from breakpoints
    int interval_start = 0;
    for (int j = 0; j < n_breaks; ++j) {
      x[i].push_back({interval_start, break_pos[j]-1, current_parent});
      interval_start = break_pos[j];
      current_parent = (current_parent == parent0) ? parent1 : parent0;
    }
    x[i].push_back({interval_start, contig_lengths[i], current_parent});
    
  }  // end i loop
  
  // push new haplotype
  intervals.push_back(x);
  
}

//------------------------------------------------
// print coords. Setting chrom = 0 prints all chromosomes
void Tree_node::print_intervals(int haplo, int chrom) {
  
  if (haplo > n_haplotypes) {
    Rcpp::stop("error in Tree_node::print_intervals, haplo greater than n_haplotypes");
  }
  if (chrom == 0) {
    for (int i = 0; i < int(intervals[0].size()); ++i) {
      print_intervals(haplo, i+1);
    }
  } else {
    Rcpp::Rcout << "Chrom" << chrom << ":\n";
    for (const auto & x : intervals[haplo][chrom-1]) {
      Rcpp::Rcout << "[" << x[0] << ", " << x[1] << "] : " << x[2] << "\n";
    }
  }
  
}
