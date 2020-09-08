
#include "misc_v9.h"

#include <vector>

//------------------------------------------------
#ifdef RCPP_ACTIVE
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// [[Rcpp::export]]
void prune_transmission_record_cpp(Rcpp::List args);
#endif

//------------------------------------------------
void add_to_pop(std::vector<int> &inoc_IDs, const std::vector<std::pair<int, int>> &first_IDs,
                std::vector<std::vector<int>> &pop, bool print_pop = false);
