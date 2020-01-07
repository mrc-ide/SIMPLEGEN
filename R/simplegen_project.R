
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
                  sample_details = NULL)
  class(project) <- "simplegen_project"
  
  # return
  invisible(project)
}

#------------------------------------------------
# overload print() function for simplegen_project
#' @method print simplegen_project
#' @export
print.simplegen_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for simplegen_project
#' @method summary simplegen_project
#' @export
summary.simplegen_project <- function(x, ...) {
  
  message("TODO - some default print method")
  
  # return invisibly
  invisible(x)
}

