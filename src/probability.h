
#pragma once

//#include <vector>
//#include <random>
//#include <math.h>

#include <dust/random/random.hpp>
#include "misc_v17.h"


// #####################################
// #        DISCRETE UNIVARIATE        #
// #####################################

//------------------------------------------------
// sample single value from given probability vector (that sums to p_sum).
// Starting in probability_v10 the first value returned from this vector is 0
// rather than 1 (i.e. moving to C++-style zero-based indexing)
int sample1(const std::vector<double> &p, double p_sum = 1.0);
int sample1(const std::vector<int> &p, int p_sum);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal
// probability
int sample2(int a, int b);

//------------------------------------------------
// draw from Bernoulli(p) distribution. NB, found to be considerably faster than
// doing a single binomial draw
bool rbernoulli1(double p);

//------------------------------------------------
// draw from binomial(N, p) distribution
int rbinom1(int N, double p);

//------------------------------------------------
// draw from Poisson distribution with rate lambda
int rpois1(double lambda);

//------------------------------------------------
// draw from Geometric(p) distribution, with mean (1-p)/p
int rgeom1(double p);


// #######################################
// #        DISCRETE MULTIVARIATE        #
// #######################################

//------------------------------------------------
// equivalent to sample2, but draws n values without replacement. Values are
// returned in order of increasing value
std::vector<int> sample4(int n, int a, int b);



// #######################################
// #        CONTINUOUS UNIVARIATE        #
// #######################################

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1]
// (NB, need to work out whether 0 and 1 exactly are possible)
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b]
// (NB, need to work out whether limits a and b exactly are possible)
double runif1(double a, double b);

//------------------------------------------------
// draw from Beta(alpha, beta) distribution
double rbeta1(double alpha, double beta);

//------------------------------------------------
// draw from Gamma distribution with shape alpha and rate beta
double rgamma1(double alpha, double beta);


// ######################
// #        MISC        #
// ######################

//------------------------------------------------
template<class TYPE>
void reshuffle(std::vector<TYPE> &x) {
  int rnd1;
  TYPE tmp1;
  int n = int(x.size());
  for (int i = 0; i < n; ++i) {
    
    // draw random index from i to end of vector
    rnd1 = floor(runif1(i, n));
    
    // swap for value at position i
    tmp1 = x[rnd1];
    x[rnd1] = x[i];
    x[i] = tmp1;
  }
}


