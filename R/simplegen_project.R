
#------------------------------------------------
#' @title Define new SIMPLEGEN project
#'
#' @description Define a new SIMPLEGEN project. This project will hold all
#'   simulation inputs and outputs for a given analysis, and is initialised with
#'   the default values of all parameters.
#'
#' @export

simplegen_project <- function() {
  
  # create empty project
  project <- list(epi_parameters = NULL,
                  sampling_strategy = NULL,
                  epi_output = NULL,
                  sample_details = NULL,
                  relatedness = NULL,
                  true_genotypes = NULL,
                  observed_genotypes = NULL)
  
  class(project) <- "simplegen_project"
  
  # return
  invisible(project)
}

#------------------------------------------------
# overload print() function for simplegen_project
#' @noRd
print.simplegen_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for simplegen_project
#' @noRd
summary.simplegen_project <- function(x, ...) {
  
  
  
  message("TODO - some default print method")
  
}

