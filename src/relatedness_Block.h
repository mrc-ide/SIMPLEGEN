
#pragma once

#include "relatedness_Coaltracker.h"
#include <vector>
#include <set>
#include <map>

//------------------------------------------------
// class defining a set of recombinant blocks for calculating coalescence times
class relatedness_Block {
  
public:
  
  // PUBLIC OBJECTS
  int haplo_ID;
  int contig;
  int generation;
  
  std::vector<int> left;
  std::vector<int> right;
  std::vector<int> child_haplo_ID;
  
  // pointer to map of relatedness blocks
  std::map<std::pair<int, int>, relatedness_Block>* block_map_ptr;
  
  // pointer to object for tracking coalescent events
  std::vector<relatedness_Coaltracker>* coaltracker_ptr;
  
  // pointer to output record
  Rcpp::List* store_output_ptr;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  relatedness_Block() {};
  relatedness_Block(int haplo_ID, int contig, int generation, int left, int right, int child_haplo_ID,
                    std::map<std::pair<int, int>, relatedness_Block>* block_map_ptr,
                    std::vector<relatedness_Coaltracker>* coaltracker_ptr,
                    Rcpp::List* store_output_ptr);
  
  // main events
  void add_block(int left, int right, int child_haplo_ID);
  void get_overlap(int block_left, int block_right, int block_parent_haplo_ID, int block_coal_time);
  void map_parent(int block_left, int block_right, int block_parent_haplo_ID, int block_descendant_haplo_ID, int block_coal_time);
  void check_coalescence(int block_left, int block_right, int block_descendant_haplo_ID, int block_coal_time);
  void coalescence(int parent_index, int block_left, int block_right, int block_descendant_haplo_ID, int block_coal_time);
  void get_blocks_in_order();
  void print_status();
  
};
