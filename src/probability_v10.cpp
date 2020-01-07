
#include "probability_v10.h"
#include "misc_v9.h"

using namespace std;

// set random seed
random_device rd;
default_random_engine generator(rd());

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
#ifdef RCPP_ACTIVE
double runif_0_1() {
  return R::runif(0,1);
}
#else
double runif_0_1() {
  uniform_real_distribution<double> uniform_0_1(0.0,1.0);
  return uniform_0_1(generator);
}
#endif

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
#ifdef RCPP_ACTIVE
double runif1(double a, double b) {
  return R::runif(a,b);
}
#else
double runif1(double a, double b) {
  uniform_real_distribution<double> uniform_a_b(a,b);
  return uniform_a_b(generator);
}
#endif

//------------------------------------------------
// draw from Bernoulli(p) distribution
#ifdef RCPP_ACTIVE
bool rbernoulli1(double p) {
  return R::rbinom(1, p);
}
#else
bool rbernoulli1(double p) {
  bernoulli_distribution dist_bernoulli(p);
  return dist_bernoulli(generator);
}
#endif

//------------------------------------------------
// draw from binomial(N,p) distribution
#ifdef RCPP_ACTIVE
int rbinom1(int N, double p) {
  return R::rbinom(N, p);
}
#else
int rbinom1(int N, double p) {
  binomial_distribution<int> dist_binom(N, p);
  return dist_binom(generator);
}
#endif

//------------------------------------------------
// draw from multinomial(N,p) distribution, where p sums to p_sum
vector<int> rmultinom1(int N, const vector<double> &p, double p_sum) {
  int k = int(p.size());
  vector<int> ret(k);
  for (int i = 0; i < (k-1); ++i) {
    ret[i] = rbinom1(N, p[i]/p_sum);
    N -= ret[i];
    if (N == 0) {
      break;
    }
    p_sum -= p[i];
  }
  ret[k-1] = N;
  return ret;
}

//------------------------------------------------
// draw from univariate normal distribution
#ifdef RCPP_ACTIVE
double rnorm1(double mean, double sd) {
  return R::rnorm(mean, sd);
}
#else
double rnorm1(double mean, double sd) {
  normal_distribution<double> dist_norm(mean,sd);
  return dist_norm(generator);
}
#endif

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b) {
  
  // draw raw value relative to a
  double ret = rnorm1(mean, sd) - a;
  
  // reflect off boundries at 0 and (b-a)
  if (ret < 0 || ret > (b-a)) {
    // use multiple reflections to bring into range [-(b-a), 2(b-a)]
    while (ret < -(b-a)) {
      ret += 2*(b-a);
    }
    while (ret > 2*(b-a)) {
      ret -= 2*(b-a);
    }
    
    // use one more reflection to bring into range [0, (b-a)]
    if (ret < 0) {
      ret = -ret;
    }
    if (ret > (b-a)) {
      ret = 2*(b-a) - ret;
    }
  }
  
  // no longer relative to a
  ret += a;
  
  // don't let ret equal exactly a or b
  if (ret == a) {
    ret += UNDERFLO_DOUBLE;
  } else if (ret == b) {
    ret -= UNDERFLO_DOUBLE;
  }
  
  return(ret);
}

//------------------------------------------------
// draw from multivariate normal distribution with mean mu and
// variance/covariance matrix sigma*scale^2. The inputs consist of mu,
// sigma_chol, and scale, where sigma_chol is the Cholesky decomposition of
// sigma. Output values are stored in x.
void rmnorm1(vector<double> &x, const vector<double> &mu,
             const vector<vector<double>> &sigma_chol, double scale) {
  
  int d = int(mu.size());
  x = mu;
  double z;
  for (int j = 0; j < d; j++) {
    z = rnorm1();
    for (int i = j; i < d; i++) {
      x[i] += sigma_chol[i][j]*scale*z;
    }
  }
}

//------------------------------------------------
// resample a vector without replacement
// reshuffle
// DEFINED IN HEADER

//------------------------------------------------
// sample single value from given probability vector (that sums to p_sum).
// Starting in probability_v10 the first value returned from this vector is 0
// rather than 1 (i.e. moving to C++-style zero-based indexing)
int sample1(const vector<double> &p, double p_sum) {
  double rand = p_sum*runif_0_1();
  double z = 0;
  for (int i=0; i<int(p.size()); i++) {
    z += p[i];
    if (rand < z) {
      return i;
    }
  }
#ifdef RCPP_ACTIVE
  Rcpp::stop("error in sample1(), ran off end of probability vector");
#else
  stop("error in sample1(), ran off end of probability vector");
#endif
  return 0;
}
int sample1(const vector<int> &p, int p_sum) {
  int rand = sample2(1,p_sum);
  int z = 0;
  for (int i=0; i<int(p.size()); i++) {
    z += p[i];
    if (rand <= z) {
      return i;
    }
  }
#ifdef RCPP_ACTIVE
  Rcpp::stop("error in sample1(), ran off end of probability vector");
#else
  stop("error in sample1(), ran off end of probability vector");
#endif
  return 0;
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal
// probability
int sample2(int a, int b) {
  return floor(runif1(a, b+1));
}

//------------------------------------------------
// sample a series of values with replacement from given probability vector that
// sums to p_sum. Results are stored in ret (passed in by reference), and the
// number of draws is dictated by the length of this vector. Option to return
// vector in shuffled order
void sample3(vector<int> &ret, const vector<double> &p, double p_sum, bool return_shuffled) {
  int n = int(ret.size());
  int j = 0;
  for (int i = 0; i < int(p.size()); ++i) {
    int n_i = rbinom1(n, p[i]/p_sum);
    if (n_i > 0) {
      fill(ret.begin()+j, ret.begin()+j+n_i, i);
      j += n_i;
      n -= n_i;
      if (n == 0) {
        break;
      }
    }
    p_sum -= p[i];
  }
  if (return_shuffled) {
    reshuffle(ret);
  }
}

//------------------------------------------------
// equivalent to sample2, but draws n values without replacement
vector<int> sample4(int n, int a, int b) {
  vector<int> ret(n);
  int t = 0, m = 0;
  int N = b - a + 1;
  if (n > N) {
    Rcpp::stop("error in sample4(), attempt to sample more elements than are available");
  }
  for (int i = 0; i < N; ++i) {
    if (sample2(1, N-t) <= (n-m)) {
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

//------------------------------------------------
// draw from gamma(shape,rate) distribution
#ifdef RCPP_ACTIVE
double rgamma1(double shape, double rate) {
  return(R::rgamma(shape, 1.0/rate));
}
#else
double rgamma1(double shape, double rate) {
  gamma_distribution<double> rgamma(shape, 1.0/rate);
  double x = rgamma(generator);
  
  // check for zero or infinite values (catches bug present in Visual Studio 2010)
  if (x == 0) {
    x = UNDERFLO_DOUBLE;
  }
  if ((1.0/x) == 0) {
    x = 1.0/UNDERFLO_DOUBLE;
  }
  
  return x;
}
#endif


//------------------------------------------------
// draw from beta(shape1,shape2) distribution
#ifdef RCPP_ACTIVE
double rbeta1(double shape1, double shape2) {
  return R::rbeta(shape1, shape2);
}
#else
double rbeta1(double shape1, double shape2) {
  double x1 = rgamma1(shape1, 1.0);
  double x2 = rgamma1(shape2, 1.0);
  return x1/double(x1+x2);
}
#endif

//------------------------------------------------
// draw from Poisson distribution with rate lambda
#ifdef RCPP_ACTIVE
int rpois1(double lambda) {
  return R::rpois(lambda);
}
#else
int rpois1(double lambda) {
  Rcpp::stop("C++ version of poisson draws not coded yet!");
}
#endif

//------------------------------------------------
// draw from zero-truncated Poisson distribution with rate lambda
// mean = lambda/(1 -exp(-lambda))
#ifdef RCPP_ACTIVE
int rztpois1(double lambda) {
  double rnd1 = runif_0_1();
  double t = -log(1 - rnd1*(1 - exp(-lambda)));
  return R::rpois(lambda - t) + 1;
}
#else
int rztpois1(double lambda) {
  Rcpp::stop("C++ version of zero-truncated poisson draws not coded yet!");
}
#endif

//------------------------------------------------
// probability mass of Poisson distribution
#ifdef RCPP_ACTIVE
double dpois1(int n, double lambda, bool return_log) {
  return R::dpois(n,lambda,return_log);
}
#else
double dpois1(int n, double lambda, bool return_log) {
  double ret = n*log(lambda) - lambda - lgamma(n+1);
  if (!return_log) {
    ret = exp(ret);
  }
  return ret;
}
#endif

//------------------------------------------------
// draw from symmetric dichlet(alpha) distribution of length n
vector<double> rdirichlet1(double alpha, int n) {
  vector<double> ret(n);
  double retSum = 0;
  for (int i=0; i<n; i++) {
    ret[i] = rgamma1(alpha,1.0);
    retSum += ret[i];
  }
  for (int i=0; i<n; i++) {
    ret[i] /= retSum;
  }
  return ret;
}

//------------------------------------------------
// draw from Geometric(p) distribution, with mean (1-p)/p
#ifdef RCPP_ACTIVE
int rgeom1(const double p) {
  return R::rgeom(p);
}
#else
int rgeom1(const double p) {
  geometric_distribution<int> dist_geom(p);
  return dist_geom(generator);
}
#endif

//------------------------------------------------
// draw from exponential(r) distribution
#ifdef RCPP_ACTIVE
double rexp1(const double r) {
  return R::rexp(1/r);
}
#else
double rexp1(const double r) {
  exponential_distribution<double> dist_exponential(r);
  return dist_exponential(generator);
}
#endif

//------------------------------------------------
// binomial coeffiient n choose k
int choose(int n, int k) {
  return int(exp(lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1)));
}

//------------------------------------------------
// binomial coeffiient n choose k, returned in log space
double lchoose(int n, int k) {
  return lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1);
}
