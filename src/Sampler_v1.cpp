
#include "Sampler_v1.h"
#include "misc_v6.h"
#include "probability_v2.h"

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
  for (int i=0; i<int(p.size()); ++i) {
    int n_i = rbinom1(n_remaining, p[i]/p_remaining);
    fill(x.begin()+n_cum, x.begin()+n_cum+n_i, i);
    n_cum += n_i;
    n_remaining -= n_i;
    p_remaining -= p[i];
    if (n_remaining == 0) {
      break;
    }
  }
  
  // re-shuffle x to break order
  reshuffle(x);
  
}

//------------------------------------------------
// draw value
int Sampler::draw() {
  
  // reset x if needed
  if (index == (n-1)) {
    reset();
    index = 0;
  }
  
  // get next value
  return x[index++];
  
}
