
#pragma once

//#include <vector>
//#include <random>
//#include <math.h>

#include <dust/r/random.hpp>

#include "misc_v17.hpp"


// #####################################
// #        DISCRETE UNIVARIATE        #
// #####################################
/*
//------------------------------------------------
int sample1(const std::vector<double> &p, double p_sum = 1.0);
int sample1(const std::vector<int> &p, int p_sum);

//------------------------------------------------
int sample2(int a, int b);

//------------------------------------------------
bool rbernoulli1(double p);

//------------------------------------------------
int rbinom1(int N, double p);

//------------------------------------------------
int rhyper1(int m, int n, int k);

//------------------------------------------------
double dhyper1(double x, int m, int n, int k, bool return_log = true);

//------------------------------------------------
int rpois1(double lambda);

//------------------------------------------------
int rztpois1(double lambda);

//------------------------------------------------
double dpois1(int n, double lambda, bool return_log = true);

//------------------------------------------------
int rgeom1(const double p);


// #######################################
// #        DISCRETE MULTIVARIATE        #
// #######################################

//------------------------------------------------
void sample3(std::vector<int> &ret, const std::vector<double> &p,
             double p_sum = 1.0, bool return_shuffled = true);

//------------------------------------------------
std::vector<int> sample4(int n, int a, int b);

//------------------------------------------------
std::vector<int> rmultinom1(int N, const std::vector<double> &p, double p_sum = 1.0);

//------------------------------------------------
double dmultinom1(const std::vector<int> &x, int x_sum, const std::vector<double> &p, double p_sum);


// #######################################
// #        CONTINUOUS UNIVARIATE        #
// #######################################

//------------------------------------------------
double runif_0_1();

//------------------------------------------------
double runif1(double a, double b);

//------------------------------------------------
double rnorm1(double mean = 0.0, double sd = 1.0);

//------------------------------------------------
double dnorm1(double x, double mean = 0.0, double sd = 1.0, bool return_log = true);

//------------------------------------------------
double rnorm1_interval(double mean, double sd, double a, double b);

//------------------------------------------------
double rlnorm1(double meanlog, double sdlog);

//------------------------------------------------
double dlnorm1(double x, double meanlog, double sdlog, bool return_log = true);

//------------------------------------------------
double rgamma1(double shape, double rate);

//------------------------------------------------
double dgamma1(double x, double shape, double rate, bool return_log = true);

//------------------------------------------------
double rinvgamma1(double shape, double scale);

//------------------------------------------------
double dinvgamma1(double x, double shape, double scale, bool return_log = true);

//------------------------------------------------
double rbeta1(double shape1, double shape2);

//------------------------------------------------
double dbeta1(double x, double shape1, double shape2, bool return_log = true);

//------------------------------------------------
double rexp1(const double r);


// #########################################
// #        CONTINUOUS MULTIVARIATE        #
// #########################################

//------------------------------------------------
void rmnorm1(std::vector<double> &x, const std::vector<double> &mu,
             const std::vector<std::vector<double>> &sigma_chol, double scale = 1.0);

//------------------------------------------------
double dmnorm1(const std::vector<double> &x,
               const std::vector<double> &mu,
               double logdet,
               const std::vector<std::vector<double>> &chol_inverse);

//------------------------------------------------
double dmnorm2(const std::vector<double> &x,
               const std::vector<double> &mu,
               const std::vector<std::vector<double>> &sigma);

//------------------------------------------------
double dinvwish1(const std::vector<std::vector<double>> &sigma_inv,
                 const std::vector<std::vector<double>> &sigma_chol,
                 const std::vector<std::vector<double>> &psi,
                 const std::vector<std::vector<double>> &psi_chol,
                 double nu);

//------------------------------------------------
double dinvwish2(const std::vector<std::vector<double>> &sigma,
                 const std::vector<std::vector<double>> &psi,
                 double nu);

//------------------------------------------------
std::vector<double> rdirichlet1(double alpha, int n);

//------------------------------------------------
std::vector<double> rdirichlet2(std::vector<double> &alpha);


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

//------------------------------------------------
int choose(int n, int k);

//------------------------------------------------
double lchoose(int n, int k);
*/


