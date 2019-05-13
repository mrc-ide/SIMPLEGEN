#------------------------------------------------
#' @title Simulating Plasmodium Epidemiological and Genetic Data
#'
#' @description Forwards-in time simulation of Plasmodium falciparum genetic
#'   data can be computationally intensive, as many genotypes are tracked but
#'   ultimately lost. SIMPLEGEN avoids this problem by splitting simulation into
#'   three steps - first, simulating transmission under a simple
#'   individual-based model, second, pruning the infection history to focus on
#'   nodes that contribute to the final sample, and third, simulating genetic
#'   data. A secondary advantage is that any third-party epidemiological
#'   simulator can be used for the first step as long as it outputs in the
#'   correct format. The major limitation of SIMPLEGEN is that assumes the
#'   epidemiological and genetic processes are seperable, and hence is limited
#'   to neutral variation, and cannot model selection.
#'
#' @docType package
#' @name SIMPLEGEN
NULL

#------------------------------------------------
#' @useDynLib SIMPLEGEN, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("SIMPLEGEN", libpath)  #' @skip
}
