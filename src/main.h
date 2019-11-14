
#include "misc_v8.h"

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
#endif
