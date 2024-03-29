
#include "Sampler_v6.h"

using namespace std;

//------------------------------------------------
// constructor
Sampler::Sampler(vector<double> &p, int n) {
  
  // initialise objects
  this->n = n;
  index = 0;
  this->p = p;
  sum_p = sum(p);
  x = vector<int>(n);
  
  // checks on inputs
  for (size_t i = 0; i < p.size(); ++i) {
    if (p[i] < 0.0) {
      cpp11::stop("error when initialising sampler: p values less than 0");
    }
  }
  if (sum_p == 0) {
    cpp11::stop("error when initialising sampler: sum over p is zero");
  }
  
  // populate x
  reset();
}

//------------------------------------------------
// reset x with fresh set of draws
void Sampler::reset() {
  
  // loop through p distribution, draw how many of each value are observed. This
  // results in a series of draws from p, sorted in increasing order
  int n_cum = 0;
  int n_remaining = n;
  double p_remaining = sum_p;
  for (size_t i = 0; i < p.size(); ++i) {
    double q = p[i] / p_remaining;
    if (q > 1.0) {  // deal with rounding issue causing probabilities > 1
      q = 1.0;
    }
    int n_i = rbinom1(n_remaining, q);
    fill(x.begin() + n_cum, x.begin() + n_cum + n_i, i);
    n_cum += n_i;
    n_remaining -= n_i;
    p_remaining -= p[i];
    if (n_remaining == 0) {
      break;
    }
  }
  
  // shuffle x to break order
  reshuffle(x);
}

//------------------------------------------------
// draw value
int Sampler::draw() {
  
  // reset x if needed
  if (index == (n - 1)) {
    reset();
    index = 0;
  }
  
  // get next value
  return x[index++];
}

//------------------------------------------------
// print summary of values
void Sampler::print_summary() {
  
  print("n =", n);
  print("index =", index);
  print("sum_p =", sum_p);
  print("p:");
  print_vector(p);
}
