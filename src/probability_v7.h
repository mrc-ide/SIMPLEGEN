
#pragma once

#include <vector>
#include <random>
#include <math.h>

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(double a, double b);

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(double p);

//------------------------------------------------
// draw from binomial(N,p) distribution
int rbinom1(int N, double p);

//------------------------------------------------
// draw from multinomial(N,p) distribution, where p sums to p_sum
std::vector<int> rmultinom1(int N, const std::vector<double> &p, double p_sum = 1.0);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean = 0.0, double sd = 1.0);

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b);

//------------------------------------------------
// draw from multivariate normal distribution with mean mu and
// variance/covariance matrix sigma*scale^2. The inputs consist of mu,
// sigma_chol, and scale, where sigma_chol is the Cholesky decomposition of
// sigma. Output values are stored in x.
void rmnorm1(std::vector<double> &x, std::vector<double> &mu, std::vector<std::vector<double>> &sigma_chol, double scale = 1.0);

//------------------------------------------------
// resample a vector without replacement
template<class TYPE>
void reshuffle(std::vector<TYPE> &x) {
  int rnd1;
  TYPE tmp1;
  int n = int(x.size());
  for (int i=0; i<n; i++) {
    
    // draw random index from i to end of vector. Note that although runif
    // returns a double, by forcing to int we essentially round this value down
    // to nearest int.
    rnd1 = runif1(i, n);
    
    // swap for value at position i
    tmp1 = x[rnd1];
    x[rnd1] = x[i];
    x[i] = tmp1;
  }
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal
// probability
int sample2(int a, int b);

//------------------------------------------------
// sample single value from given probability vector (that sums to p_sum)
int sample1(std::vector<double> &p, double p_sum = 1.0);

//------------------------------------------------
// sample a series of values from given probability vector (that sums to p_sum).
// Results are stored in ret (passed in by reference), and the number of draws
// is dictated by the length of this vector. Option to return vector in shuffled
// order
void sample3(std::vector<int> &ret, std::vector<double> &p, double p_sum = 1.0, bool return_shuffled = true);

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(shape1,shape2) distribution
double rbeta1(double shape1, double shape2);

//------------------------------------------------
// draw from zero-truncated Poisson distribution with rate lambda
// mean = lambda/(1 -exp(-lambda))
int rztpois1(double lambda);

//------------------------------------------------
// probability mass of Poisson distribution
double dpois1(int n, double lambda, bool return_log);

//------------------------------------------------
// draw from symmetric dichlet(alpha) distribution of length n
std::vector<double> rdirichlet1(double alpha, int n);

//------------------------------------------------
// draw from Geometric(p) distribution, with mean (1-p)/p
int rgeom1(const double p);

//------------------------------------------------
// draw from exponential(r) distribution
double rexp1(const double r);

//------------------------------------------------
// binomial coeffiient n choose k
int choose(int n, int k);

//------------------------------------------------
// binomial coeffiient n choose k, returned in log space
double lchoose(int n, int k);
