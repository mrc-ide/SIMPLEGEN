
#pragma once

#include "Sampler_v4.h"

#include <map>
#include <vector>

//------------------------------------------------
// class defining a single node in the inoculation tree
class Tree_node {
  
public:
  
  // PUBLIC OBJECTS
  
  // basic properties
  int t;
  std::vector<int> contig_lengths;
  int n_contigs;
  int n_lineages;
  std::vector<int> lineage_ID_vec;
  std::vector<double> lineage_density;
  
  // pointers to external objects
  Sampler *sampler_oocyst_ptr;
  Sampler *sampler_hepatocyte_ptr;
  std::map<int, Tree_node> *tree_ptr;
  
  // main object, specifying ancestry of each lineage broken down into
  // intervals. Nested vectors are over:
  //     1) number of lineages in this inoculation
  //     2) number of contigs
  //     3) number of disjoint intervals within a contig
  //     4) three values giving interval start, interval end, and ancestral
  //        lineage_ID that applies throughout this interval
  std::vector<std::vector<std::vector<std::vector<int>>>> intervals;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Tree_node() {};
  Tree_node(int t, const std::vector<int> &contig_lengths, Sampler &sampler_oocyst,
            Sampler &sampler_hepatocyte, std::map<int, Tree_node> &tree);
  
  // other methods
  void init(int t, const std::vector<int> &contig_lengths, Sampler &sampler_oocyst,
            Sampler &sampler_hepatocyte, std::map<int, Tree_node> &tree);
  void draw_lineages_denovo(int &lineage_ID, double alpha);
  void draw_lineages_recombine(int &lineage_ID, const std::vector<int> &inoc_IDs, double r, double alpha);
  void draw_intervals(int parent0, int parent1, double r);
  void print_intervals(int lineage = 0, int contig = 0);
  
};
