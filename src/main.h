
#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
Rcpp::List indiv_sim_cpp(Rcpp::List args);
#else
int indiv_sim_cpp();
#endif
