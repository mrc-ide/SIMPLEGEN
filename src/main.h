
#include "misc_v11.h"

#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List prune_transmission_record_cpp(Rcpp::List args);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List sim_haplotype_tree_cpp(Rcpp::List args);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List get_haplotype_relatedness_cpp(Rcpp::List args);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List sim_block_tree_cpp(Rcpp::List args);
