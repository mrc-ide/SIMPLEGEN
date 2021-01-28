
#include "Tree_node.h"
#include "misc_v9.h"
#include "probability_v10.h"

using namespace std;


//------------------------------------------------
// constructor calling init()
Tree_node::Tree_node(int t, const std::vector<int> &contig_lengths, Sampler &sampler_oocyst,
                     Sampler &sampler_hepatocyte, map<int, Tree_node> &tree) {
  init(t, contig_lengths, sampler_oocyst, sampler_hepatocyte, tree);
}

//------------------------------------------------
// initialise node by copying over basic properties
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
// draw lineages de novo (without ancestors)
void Tree_node::draw_lineages_denovo(int &lineage_ID, double alpha) {
  
  // initialise with a single lineage
  n_lineages = 1;
  lineage_ID_vec = vector<int>(1, lineage_ID++);
  lineage_density = vector<double>(1);
  lineage_density[0] = rgamma1(alpha, 1.0);
  
  // initialise intervals with no ancestry (-1) over every contig
  intervals = vector<vector<vector<vector<int>>>>(1, vector<vector<vector<int>>>(n_contigs));
  for (int i = 0; i < n_contigs; ++i) {
    intervals[0][i].push_back({0, contig_lengths[i] - 1, -1});
  }
  
}

//------------------------------------------------
// draw lineages via recombination using ancestral inoculation IDs
void Tree_node::draw_lineages_recombine(int &lineage_ID, const std::vector<int> &inoc_IDs,
                                        double r, double alpha) {
  
  // get the complete vector of lineage IDs and densities by concatenating
  // over the ancestral inoculations
  vector<int> parental_lineage_IDs;
  vector<double> parental_lineage_densities;
  for (int i = 2; i < int(inoc_IDs.size()); ++i) {  // start at i=2 because first two values are the inoculation ID of this node and the timing
    push_back_multiple(parental_lineage_IDs, (*tree_ptr)[inoc_IDs[i]].lineage_ID_vec);
    push_back_multiple(parental_lineage_densities, (*tree_ptr)[inoc_IDs[i]].lineage_density);
  }
  double sum_parental_lineage_densities = sum(parental_lineage_densities);
  
  // draw number of oocysts
  int n_oocysts = sampler_oocyst_ptr->draw() + 1;
  
  // draw the parental pairs that make up each oocyst. Each parent is a lineage_ID
  // taken from one of the ancestral inoculations. lineage_IDs are drawn at random
  // with probaiblities given by their relative densities.
  vector<pair<int, int>> oocyst_parents(n_oocysts);
  for (int i = 0; i < n_oocysts; ++i) {
    int samp1 = sample1(parental_lineage_densities, sum_parental_lineage_densities);
    int samp2 = sample1(parental_lineage_densities, sum_parental_lineage_densities);
    oocyst_parents[i] = make_pair(parental_lineage_IDs[samp1], parental_lineage_IDs[samp2]);
  }
  
  // draw number of infected hepatocytes
  int n_hepatocytes = sampler_hepatocyte_ptr->draw() + 1;
  
  // draw the number of infected hepatocytes that originate from each oocysts
  // product. This is a simple multinomial draw, implemented as a series of
  // binomials.
  vector<vector<int>> n_products(n_oocysts, vector<int>(4));
  int N = n_hepatocytes;
  int n_allocated = 0;
  int n_trials_remaining = 4*n_oocysts;
  for (int i = 0; i < n_oocysts; ++i) {
    for (int j = 0; j < 4; ++j) {
      n_products[i][j] = rbinom1(N - n_allocated, 1/double(n_trials_remaining));
      n_allocated += n_products[i][j];
      n_trials_remaining--;
    }
  }
  
  // if an oocyst is clonal then group all products together
  for (int i = 0; i < n_oocysts; ++i) {
    if (oocyst_parents[i].first == oocyst_parents[i].second) {
      n_products[i][0] = sum(n_products[i]);
      for (int j = 1; j < 4; ++j) {
        n_products[i][j] = 0;
      }
    }
  }
  
  // if oocysts have identical clonal parents then group oocysts together
  for (int i = 1; i < n_oocysts; ++i) {
    if (oocyst_parents[i].first == oocyst_parents[i].second) {
      for (int j = 0; j < i; ++j) {
        if (oocyst_parents[i].first == oocyst_parents[j].first && oocyst_parents[i].second == oocyst_parents[j].second) {
          n_products[j][0] += n_products[i][0];
          n_products[i][0] = 0;
        }
      }
    }
  }
  
  // for each oocysts product that leads to infected hepatocytes, create a new
  // interval map by recombination and draw the density of this lineage. Then
  // push to invervals object
  n_lineages = 0;
  for (int i = 0; i < n_oocysts; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (n_products[i][j] == 0) {
        continue;
      }
      
      // update ID and density vector
      n_lineages++;
      lineage_ID_vec.push_back(lineage_ID++);
      lineage_density.push_back(rgamma1(alpha*n_products[i][j], 1.0));
      
      // draw recombinant products and push to intervals object
      draw_intervals(oocyst_parents[i].first, oocyst_parents[i].second, r);
    }
  }
  
}

//------------------------------------------------
// create intervals by recombination of two parental IDs, and push to intervals
// array
void Tree_node::draw_intervals(int parent0, int parent1, double r) {
  
  // loop through contigs
  vector<vector<vector<int>>> x(n_contigs);
  for (int i = 0; i < n_contigs; ++i) {
    
    // deal with clonal case
    if (parent0 == parent1) {
      x[i].push_back({0, contig_lengths[i] - 1, parent0});
      continue;
    }
    
    // draw starting parent with equal probability
    int current_parent = rbernoulli1(0.5) ? parent1 : parent0;
    
    // draw number of breakpoints. Note, the following method works even if
    // n_breaks is zero
    int n_breaks = rpois1(r*contig_lengths[i]);
    
    // draw location of breakpoints ensuring no duplicates
    vector<int> break_pos(n_breaks);
    for (int j = 0; j < n_breaks; ++j) {
      break_pos[j] = sample2(1, contig_lengths[i] - 1);
    }
    sort(break_pos.begin(), break_pos.end());
    remove_duplicates(break_pos);
    
    // create intervals from breakpoints
    int interval_start = 0;
    for (int j = 0; j < int(break_pos.size()); ++j) {
      x[i].push_back({interval_start, break_pos[j]-1, current_parent});
      interval_start = break_pos[j];
      current_parent = (current_parent == parent0) ? parent1 : parent0;
    }
    x[i].push_back({interval_start, contig_lengths[i] - 1, current_parent});
    
  }  // end i loop
  
  // push new lineage
  intervals.push_back(x);
  
}

//------------------------------------------------
// print lineage intervals. Setting contig = 0 prints all contigs
void Tree_node::print_intervals(int lineage, int contig) {
  
  if (lineage > n_lineages) {
    Rcpp::stop("error in Tree_node::print_intervals, lineage greater than n_lineages");
  }
  if (contig == 0) {
    for (int i = 0; i < int(intervals[0].size()); ++i) {
      print_intervals(lineage, i + 1);
    }
  } else {
    Rcpp::Rcout << "Contig" << contig << ":\n";
    for (const auto & x : intervals[lineage][contig - 1]) {
      Rcpp::Rcout << "[" << x[0] << ", " << x[1] << "] : " << x[2] << "\n";
    }
  }
  
}
