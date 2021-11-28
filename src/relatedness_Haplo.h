
#pragma once

#include <vector>
#include <map>

//------------------------------------------------
// class defining a single haplotype for relatedness calculation
class relatedness_Haplo {
  
public:
  
  // PUBLIC OBJECTS
  std::vector<int> descendant_IDs;
  std::vector<int> generations;
  std::vector<int> recombinations;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  relatedness_Haplo() {};
  relatedness_Haplo(int haplo_ID);
  
  // main events
  void initialise(int haplo_ID);
  void merge(relatedness_Haplo m);
  void increment_generations();
  void increment_recombinations();
  void print_status();
  
};
