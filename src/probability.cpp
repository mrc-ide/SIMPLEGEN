
#include "probability.hpp"

using namespace std;

using rng_state_type = dust::random::generator<double>;
auto rng_state = dust::random::seed<rng_state_type>(12345);


// #####################################
// #        DISCRETE UNIVARIATE        #
// #####################################

//------------------------------------------------
int sample1(const std::vector<double> &p, double p_sum) {
  double rand = p_sum*runif_0_1();
  double z = 0;
  for (int i = 0; i < int(p.size()); i++) {
    z += p[i];
    if (rand < z) {
      return i;
    }
  }
  cpp11::stop("error in sample1(), ran off end of probability vector");
  return 0;
}
int sample1(const std::vector<int> &p, int p_sum) {
  int rand = sample2(1, p_sum);
  int z = 0;
  for (int i = 0; i < int(p.size()); i++) {
    z += p[i];
    if (rand <= z) {
      return i;
    }
  }
  cpp11::stop("error in sample1(), ran off end of probability vector");
  return 0;
}

//------------------------------------------------
int sample2(int a, int b) {
  return floor(runif1(a, b + 1));
}

//------------------------------------------------
bool rbernoulli1(double p) {
  return (runif_0_1() < p);
}

//------------------------------------------------
int rbinom1(int N, double p) {
  return dust::random::binomial<double>(rng_state, N, p);
}

//------------------------------------------------
int rpois1(double lambda) {
  return dust::random::poisson<double>(rng_state, lambda);
}

//------------------------------------------------
int rgeom1(double p) {
  if (p == 0) {
    cpp11::stop("cannot draw from rgeom1 with p == 0");
  } else if (p == 1) {
    return 0;
  }
  // draw uniformly from CDF and project back
  double y = runif_0_1();
  return floor(log(1 - y) / log(1 - p));
}


// #######################################
// #        DISCRETE MULTIVARIATE        #
// #######################################

//------------------------------------------------
std::vector<int> sample4(int n, int a, int b) {
  std::vector<int> ret(n);
  int t = 0, m = 0;
  int N = b - a + 1;
  if (n > N) {
    cpp11::stop("error in sample4(), attempt to sample more elements than are available");
  }
  for (int i = 0; i < N; ++i) {
    if (sample2(1, N - t) <= (n - m)) {
      ret[m] = a + i;
      m++;
      if (m == n) {
        break;
      }
    }
    t++;
  }
  return ret;
}


// #######################################
// #        CONTINUOUS UNIVARIATE        #
// #######################################

//------------------------------------------------
double runif_0_1() {
  return dust::random::uniform(rng_state, 0.0, 1.0);
}

//------------------------------------------------
double runif1(double a, double b) {
  return dust::random::uniform(rng_state, a, b);
}

//------------------------------------------------
double rbeta1(double alpha, double beta) {
  double r1 = rgamma1(alpha, 1.0);
  double r2 = rgamma1(beta, 1.0);
  return r1 / (r1 + r2);
}

//------------------------------------------------
double rgamma1(double alpha, double beta) {
  return dust::random::gamma(rng_state, alpha, 1.0 / beta);
}


