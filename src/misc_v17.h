
#pragma once

#include <numeric>
#include <iostream>
#include <vector>
#include <chrono>
#include <cfloat>
#include <fstream>
#include <sstream>

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
const int OVERFLO_INT = INT_MAX / 100;
const double OVERFLO_DOUBLE = DBL_MAX / 100;
const double UNDERFLO_DOUBLE = DBL_MIN / 100;

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed
// to be -inf.
double log_sum(double logA, double logB);

//------------------------------------------------
// add values in a vector in log space
double log_sum_vec(const std::vector<double> &x);

//------------------------------------------------
// basic sum over elements in a vector
template<class TYPE>
TYPE sum(const std::vector<TYPE> &x_vec) {
  TYPE ret = 0;
  for (const auto & x : x_vec) {
  	ret += x;
  }
  return ret;
}

//------------------------------------------------
// sum boolean values and return integer. Note that sum function above does not
// work for boolean values because the return type is also forced to boolean,
// and therefore returns 1 or 0
int sum_bool(const std::vector<bool> &x_vec);

//------------------------------------------------
// mean of vector
template<class TYPE>
double mean(const std::vector<TYPE> &x) {
  return sum(x) / double(x.size());
}

//------------------------------------------------
// min of vector
template<class TYPE>
TYPE min_vec(const std::vector<TYPE> &x) {
  return *min_element(x.begin(), x.end());
}

//------------------------------------------------
// max of vector
template<class TYPE>
TYPE max_vec(std::vector<TYPE> x) {
  return *max_element(x.begin(), x.end());
}

//------------------------------------------------
// square function
template<class TYPE>
TYPE sq(const TYPE x) {
  return x*x;
}

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, int to, int by = 1);

//------------------------------------------------
// Euclidian distance between points in 2 dimensions
template<class TYPE>
double dist_euclid_2d(TYPE x1, TYPE y1, TYPE x2, TYPE y2) {
  return sqrt(sq(x1 - x2) + sq(y1 - y2));
}

//------------------------------------------------
// push back multiple values to vector
template<class TYPE>
void push_back_multiple(std::vector<TYPE> &lhs, const std::vector<TYPE> &rhs) {
  lhs.insert(lhs.end(), rhs.begin(), rhs.end());
}

//------------------------------------------------
// erase particular element of a vector using efficient method. Warning - does
// not preserve original order of vector
template<class TYPE>
void quick_erase(std::vector<TYPE> &v, int index) {
  v[index] = v.back();
  v.pop_back();
}

//------------------------------------------------
// erase-remove idiom, for erasing all instances of a particular value from
// container
template<class TYPE>
void erase_remove(std::vector<TYPE> &v, TYPE x) {
  v.erase(remove(v.begin(), v.end(), x), end(v));
}

//------------------------------------------------
// test whether value can be found in vector
template<class TYPE>
bool is_in_vector(TYPE x, const std::vector<TYPE> &v) {
  return find(v.begin(), v.end(), x) != v.end();
}

//------------------------------------------------
// get the order of a vector in terms of element indices. Equivalent to R
// function order() except using zero-based indexing. Identical values are
// returned in the order they are found, for example the order of {2,1,1} is
// {1,2,0}.
template <typename TYPE>
std::vector<int> get_order(const std::vector<TYPE> &v) {
  
  // initialize original index locations
  std::vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v using std::stable_sort instead
  // of std::sort to avoid unnecessary index re-orderings when v contains
  // elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](int i1, int i2) {return v[i1] < v[i2];});
  
  return idx;
}

//------------------------------------------------
// return unique values in a vector
template<class TYPE>
std::vector<TYPE> unique(const std::vector<TYPE> &v) {
  std::vector<TYPE> ret;
  for (auto x : v) {
    if (!is_in_vector(x, ret)) {
      ret.push_back(x);
    }
  }
  return ret;
}

//------------------------------------------------
// update timer and optionally print time difference
double chrono_timer(std::chrono::high_resolution_clock::time_point &t0,
                    std::string message_before = "completed in ",
                    std::string message_after = "",
                    bool print_diff = true);

//------------------------------------------------
// helper function for printing a single value or series of values
template<typename TYPE>
void print(TYPE x) {
  std::cout << x << "\n";
}

template<typename TYPE, typename... Args>
void print(TYPE first, Args... args) {
  std::cout << first << " ";
  print(args...);
}

//------------------------------------------------
// helper function for printing contents of a vector or set
template<class TYPE>
void print_vector(const TYPE &x) {
  for (auto i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a matrix
template<class TYPE>
void print_matrix(const std::vector<std::vector<TYPE>> &x) {
  for (int i = 0; i < x.size(); ++i) {
    print_vector(x[i]);
  }
  print("");
}

//------------------------------------------------
// helper function for printing contents of a 3D array
template<class TYPE>
void print_array(const std::vector<std::vector<std::vector<TYPE>>> &x) {
  for (int i = 0; i < x.size(); ++i) {
    std::cout << "--- slice " << i + 1 << " ---\n";
    print_matrix(x[i]);
  }
  print("");
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(int n = 50, std::string title = "");

//------------------------------------------------
// print "foo", with optional number e.g. "foo2"
void foo(int n = 0);

//------------------------------------------------
// print "bar", with optional number e.g. "bar2"
void bar(int n = 0);

//------------------------------------------------
// print "foobar", with optional number e.g. "foobar2"
void foobar(int n = 0);

//------------------------------------------------
// read values from comma-separated text file to vector<int>
std::vector<int> file_to_vector_int(std::string file_path);

//------------------------------------------------
// read values from comma-separated text file to vector<double>
std::vector<double> file_to_vector_double(std::string file_path);

//------------------------------------------------
// read values from text file to vector<vector<double>>. Text file should be
// delimited by first line break, then comma. Lines do not all need to be same
// length, i.e. jagged arrays are allowed.
std::vector<std::vector<double>> file_to_matrix_double(std::string file_path);

//------------------------------------------------
// multivariate gamma function of x with dimension p, returned in log space
double lmvgamma_func(double x, double p);

//------------------------------------------------
// calculate Cholesky decomposition of positive definite matrix sigma
void cholesky(std::vector<std::vector<double>> &chol, const std::vector<std::vector<double>> &sigma);

//------------------------------------------------
// return (in log space) determinant of positive definite matrix. chol is the
// Cholesky decomposition of the target matrix.
double log_determinant(const std::vector<std::vector<double>> &chol);

//------------------------------------------------
// return diagonal matrix of dimension d with elements x
template<class TYPE>
std::vector<std::vector<TYPE>> diag(int d, TYPE x) {
  std::vector<std::vector<TYPE>> ret(d, std::vector<TYPE>(d));
  for (int i = 0; i < d; i++) {
    ret[i][i] = x;
  }
  return ret;
}

//------------------------------------------------
// find the inverse of the matrix M by Gauss-Jordan elimination. Note that M is
// modified inside the function, and therefore cannot be passed in by reference.
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> M);

//------------------------------------------------
// given [x,y] coordinates of points, and a series of values x_pred at which to
// return, calculate cubic spline interpolation between coordinates and save
// result into y_pred vector. Both x and x_pred must be increasing, and all
// values in x_pred must be inside (or equal to) x.
void cubic_spline(std::vector<double> &x, std::vector<double> &y,
                  std::vector<double> &x_pred, std::vector<double> &y_pred);

//------------------------------------------------
// works with the R function update_progress() to increment a progress bar by
// name. If display_updates is false then only update bar at 0% and 100%.
//void update_progress_cpp(Rcpp::List args_progress, Rcpp::Function update_progress_R,
//                         std::string pb_name, int i, int max_i, bool display_updates);
