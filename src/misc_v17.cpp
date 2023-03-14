
#include "misc_v17.h"

using namespace std;

// see header for function descriptions

//------------------------------------------------
double log_sum(double logA, double logB) {
  if (logA - logB > 100) {
    return(logA);
  } else if (logB - logA > 100) {
    return(logB);
  }
  return (logA < logB) ? logB + log(1 + exp(logA - logB)) : logA + log(1 + exp(logB - logA));
}

//------------------------------------------------
double log_sum_vec(const std::vector<double> &x) {
  double x_max = max_vec(x);
  if (!isfinite(x_max)) {
    return x_max;
  }
  double ret = 0.0;
  for (auto z : x) {
    ret += exp(z - x_max);
  }
  ret = x_max + log(ret);
  return ret;
}

//------------------------------------------------
int sum_bool(const vector<bool> &x_vec) {
  int ret = 0;
  for (int i = 0; i < int(x_vec.size()); ++i) {
    ret += x_vec[i];
  }
  return ret;
}

//------------------------------------------------
vector<int> seq_int(int from, int to, int by) {
  int n = floor((to - from) / double(by)) + 1;
  vector<int> ret(n, from);
  for (int i = 1; i < n; ++i) {
    from += by;
    ret[i] = from;
  }
  return ret;
}

//------------------------------------------------
double chrono_timer(chrono::high_resolution_clock::time_point &t0, string message_before,
                    string message_after, bool print_diff) {
  
  // calculate elapsed time
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t1-t0);
  double time_double = time_span.count();
  
  // print time difference
  if (print_diff) {
    message_before += to_string(time_double) + " seconds" + message_after;
    print(message_before);
  }
  
  // update timer to current time
  t0 = t1;
  
  // return time diff
  return time_double;
}

//------------------------------------------------
void print_stars(int n, string title) {
  std::cout << title << " ";
  for (int i = 0; i < n; ++i) {
    std::cout << "*";
  }
  print("");
}

//------------------------------------------------
void foo(int n) {
  if (n == 0) {
    print("foo");
  } else {
    print("foo", n);
  }
}

//------------------------------------------------
void bar(int n) {
  if (n == 0) {
    print("bar");
  } else {
    print("bar", n);
  }
}

//------------------------------------------------
void foobar(int n) {
  if (n == 0) {
    print("foobar");
  } else {
    print("foobar", n);
  }
}

//------------------------------------------------
vector<int> file_to_vector_int(string file_path) {
  
  // initialise return object
  vector<int> ret;
  
  // read in values from comma-separated file
  ifstream infile(file_path);
  std::string line1, line2;
  int x;
  while (getline(infile, line1)) {
    istringstream ss(line1);
    while (getline(ss, line2, ',')) {
      if (line2.size() > 0) {
        istringstream(line2) >> x;
        ret.push_back(x);
      }
    }
  }
  
  return ret;
}

//------------------------------------------------
vector<double> file_to_vector_double(string file_path) {
  
  // initialise return object
  vector<double> ret;
  
  // read in values from comma-separated file
  ifstream infile(file_path);
  std::string line1, line2;
  double x;
  while (getline(infile, line1)) {
    istringstream ss(line1);
    while (getline(ss, line2, ',')) {
      if (line2.size() > 0) {
        istringstream(line2) >> x;
        ret.push_back(x);
      }
    }
  }
  
  return ret;
}

//------------------------------------------------
vector<vector<double>> file_to_matrix_double(string file_path) {
  
  // initialise return object
  vector<vector<double>> ret;
  
  // read in values from comma-separated file
  ifstream infile(file_path);
  std::string line1, line2, line3;
  double x;
  vector<double> v;
  while (getline(infile, line1)) {
    v.clear();
    istringstream ss(line1);
    while (getline(ss, line2, ',')) {
      if (line2.size() > 0) {
        istringstream(line2) >> x;
        v.push_back(x);
      }
    }
    ret.push_back(v);
  }
  
  return ret;
}

//------------------------------------------------
double lmvgamma_func(double x, double p) {
  double ret = 0.25*p*(p - 1)*log(M_PI);
  for (int i = 0; i < p; ++i) {
    ret += lgamma(x - 0.5*i);
  }
  return ret;
}

//------------------------------------------------
void cholesky(vector<vector<double>> &chol, const vector<vector<double>> &sigma) {
  
  for (int i = 0; i < int(sigma.size()); ++i) {
    for (int j = 0; j < (i+1); ++j) {
      chol[i][j] = sigma[i][j];
      if (i == j) {
        if (i > 0) {
          for (int k = 0; k < i; ++k) {
            chol[i][i] -= chol[i][k]*chol[i][k];
          }
        }
        chol[i][i] = sqrt(chol[i][i]);
      } else {
        if (j > 0) {
          for (int k = 0; k < j; ++k) {
            chol[i][j] -= chol[i][k]*chol[j][k];
          }
        }
        chol[i][j] /= chol[j][j];
      }
    }
  }
  
}

//------------------------------------------------
double log_determinant(const vector<vector<double>> &chol) {
  double ret = 0;
  for (int i = 0; i < int(chol.size()); i++) {
    ret += 2*log(chol[i][i]);
  }
  return ret;
}

//------------------------------------------------
vector<vector<double>> inverse(vector<vector<double>> M) {
  
  int d = int(M.size());
  vector<vector<double>> IM = diag(d, 1.0);
  
  double temp;
  for (int j = 0; j < d; j++) {
    for (int i = j; i < d; i++) {
      if (i == j) {
        temp = M[i][i];
        for (int k = 0; k < d; k++) {
          M[i][k] /= temp;
          IM[i][k] /= temp;
        }
      } else {
        if (M[i][j]!=0) {
          temp = M[i][j];
          for (int k = 0; k < d; k++) {
            M[i][k] -= temp * M[j][k];
            IM[i][k] -= temp * IM[j][k];
          }
        }
      }
    }
  }
  for (int j = (d - 1); j > 0; j--) {
    for (int i = (j - 1); i >= 0; i--) {
      temp = M[i][j];
      for (int k = 0; k < d; k++) {
        M[i][k] -= temp * M[j][k];
        IM[i][k] -= temp * IM[j][k];
      }
    }
  }
  
  return IM;
}

//------------------------------------------------
void cubic_spline(vector<double> &x, vector<double> &y,
                  vector<double> &x_pred, vector<double> &y_pred) {
  
  // get vector sizes
  int n = x.size();
  int n_pred = x_pred.size();
  
  // define objects for storing spline coefficients
  vector<double> c(n + 1);
  vector<double> l(n + 1);
  vector<double> mu(n + 1);
  vector<double> z(n + 1);
  vector<double> h(n);
  vector<double> b(n);
  vector<double> d(n);
  vector<double> alpha(n);
  
  // compute coefficients
  for (int i = 0; i < n; ++i) {
    h[i] = x[i + 1] - x[i];
  }
  for (int i = 1; i < n; ++i) {
    alpha[i] = (3.0 / h[i])*(y[i + 1] - y[i]) - (3.0 / h[i - 1])*(y[i] - y[i - 1]);
  }
  l[0] = 1;
  for (int i = 1; i < n; ++i) {
    l[i] = 2*(x[i + 1] - x[i - 1]) - h[i - 1]*mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i] - h[i - 1]*z[i - 1]) / l[i];
  }
  l[n] = 1;
  for (int i = (n - 1); i >= 0; --i) {
    c[i] = z[i] - mu[i]*c[i + 1];
    b[i] = (y[i + 1] - y[i]) / h[i] - h[i]*(c[i + 1] + 2*c[i]) / 3.0;
    d[i] = (c[i + 1] - c[i]) / (3.0*h[i]);
  }
  
  // save spline y-values into y_pred
  int j = 0;
  for (int i = 0; i < n_pred; ++i) {
    if (x_pred[i] > x[j + 1]) {
      j++;
    }
    y_pred[i] = y[j] + b[j]*(x_pred[i] - x[j]) + c[j]*pow(x_pred[i] - x[j], 2) + d[j]*pow(x_pred[i] - x[j], 3);
  }
  
  return;
}

//------------------------------------------------
/*
void update_progress_cpp(Rcpp::List args_progress, Rcpp::Function update_progress_R,
                         string pb_name, int i, int max_i, bool display_updates) {
  if (i == 0) {
    update_progress_R(args_progress, pb_name, 0, max_i);
  } else if (i == max_i) {
    update_progress_R(args_progress, pb_name, i, max_i);
  } else {
    int remainder = i % int(ceil(double(max_i) / 100));
    if (remainder == 0 && display_updates) {
      update_progress_R(args_progress, pb_name, i, max_i);
    }
  }
}
*/
