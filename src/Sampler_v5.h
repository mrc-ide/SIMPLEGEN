
#pragma once

#include <vector>

//------------------------------------------------
// Sampler class allowing fast random draws from any discrete distribution
// 
// We can sample from any discrete distribution by drawing a uniform value
// between 0 and 1, and then looping through the cumulative distribution until
// we hit this value. This provides a straightfoward way to sample from any
// arbitrary distribution, but it can be quite inefficient - particularly if the
// number of possible values is large. This sampler class provides a more
// efficient method by drawing a large number of samples in one go. For an
// arbitrarily large number of samples we only need to loop through the
// cumulative distribution once, counting the binomial number of observations of
// each value. These values are stored within a vector in the sampler class, and
// can be returned one at a time as needed. Once we reach the end of the vector
// we can re-draw another batch of values.
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
  void print_summary();
};
