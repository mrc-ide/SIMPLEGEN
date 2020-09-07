
#pragma once

#include "Sampler_v3.h"

#include <map>
#include <vector>

//------------------------------------------------
// class defining a single node in the pruned transmission record
class Haplo_node {
  
public:
  
  // PUBLIC OBJECTS
  
  // basic properties
  int t;
  std::vector<int> contig_lengths;
  int n_contigs;
  int n_haplotypes;
  std::vector<int> haplo_ID_vec;
  std::vector<double> haplo_density;
  
  // pointers to external objects
  Sampler *sampler_oocyst_ptr;
  Sampler *sampler_hepatocyte_ptr;
  std::map<int, Haplo_node> *haplo_tree_ptr;
  
  // main object, specifying ancestry of each haplotype broken down into
  // intervals. Nested vectors are over: 1) number of haplotypes, 2) number of
  // chromosomes, 3) number of intervals within a chromosome, 4) three values
  // giving interval start, interval end, and ancestry index (haplo_ID) that
  // applies throughout this interval
  std::vector<std::vector<std::vector<std::vector<int>>>> intervals;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Haplo_node() {};
  Haplo_node(int t, const std::vector<int> &contig_lengths, Sampler &sampler_oocyst,
             Sampler &sampler_hepatocyte, std::map<int, Haplo_node> &haplo_tree);
  
  // other methods
  void init(int t, const std::vector<int> &contig_lengths, Sampler &sampler_oocyst,
            Sampler &sampler_hepatocyte, std::map<int, Haplo_node> &haplo_tree);
  void draw_haplotypes_denovo(int &haplo_ID, double alpha);
  void draw_haplotypes_recombine(int &haplo_ID, const std::vector<int> &inoc_IDs, double r, double alpha);
  void draw_intervals(int parent0, int parent1, double r);
  void print_intervals(int haplo = 0, int chrom = 0);
  
};
