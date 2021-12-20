
#pragma once

#include "misc_v12.h"

#include <vector>
#include <set>
#include <map>

//------------------------------------------------
// class that keeps track of number of remaining uncoalesced lineages at each
// point along a contig
class relatedness_Coaltracker {
  
public:
  
  // PUBLIC OBJECTS
  int contig;
  std::vector<int> left;
  std::vector<int> right;
  std::vector<int> lineages_remaining;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  relatedness_Coaltracker() {};
  
  // main events
  void init(int contig, int left, int right, int lineages_remaining);
  void add_block(int left, int right, int lineages_remaining);
  void get_overlap(int block_left, int block_right, std::vector<int> &res_left, std::vector<int> &res_right);
  void coalescence(int block_left, int block_right);
  void tidy_blocks();
  void drop_block(int i);
  void print_status();
  
};
