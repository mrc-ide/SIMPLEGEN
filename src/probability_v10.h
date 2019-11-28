
#pragma once

#include <vector>
#include <random>
#include <math.h>

//------------------------------------------------
double runif_0_1();

//------------------------------------------------
double runif1(double a, double b);

//------------------------------------------------
bool rbernoulli1(double p);

//------------------------------------------------
int rbinom1(int N, double p);

//------------------------------------------------
std::vector<int> rmultinom1(int N, const std::vector<double> &p, double p_sum = 1.0);

//------------------------------------------------
double rnorm1(double mean = 0.0, double sd = 1.0);

//------------------------------------------------
double rnorm1_interval(double mean, double sd, double a, double b);

//------------------------------------------------
void rmnorm1(std::vector<double> &x, const std::vector<double> &mu,
             const std::vector<std::vector<double>> &sigma_chol, double scale = 1.0);

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

//------------------------------------------------
int sample1(const std::vector<double> &p, double p_sum = 1.0);

//------------------------------------------------
int sample2(int a, int b);

//------------------------------------------------
void sample3(std::vector<int> &ret, const std::vector<double> &p,
             double p_sum = 1.0, bool return_shuffled = true);

//------------------------------------------------
std::vector<int> sample4(int n, int a, int b);

//------------------------------------------------
double rgamma1(double shape, double rate);

//------------------------------------------------
double rbeta1(double shape1, double shape2);

//------------------------------------------------
int rpois1(double lambda);

//------------------------------------------------
int rztpois1(double lambda);

//------------------------------------------------
double dpois1(int n, double lambda, bool return_log);

//------------------------------------------------
std::vector<double> rdirichlet1(double alpha, int n);

//------------------------------------------------
int rgeom1(const double p);

//------------------------------------------------
double rexp1(const double r);

//------------------------------------------------
int choose(int n, int k);

//------------------------------------------------
double lchoose(int n, int k);
