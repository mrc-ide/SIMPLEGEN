
#pragma once

#include <vector>

//------------------------------------------------
// Sampler class allowing fast random draws from any discrete distribution
class Sampler {
  
public:
  
  // PUBLIC OBJECTS
  
  int n;
  int index;
  std::vector<double> p;
  double sum_p;
  std::vector<int> x;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Sampler() {};
  Sampler(std::vector<double> &p, int n);
  
  // methods
  void reset();
  int draw();
};
