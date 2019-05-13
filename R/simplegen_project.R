
#------------------------------------------------
#' @title Define new SIMPLEGEN project
#'
#' @description Define a new SIMPLEGEN project. This project will hold all
#'   simulation inputs and outputs for a given analysis, and is initialised with
#'   the default values of all parameters.
#'
#' @export

simplegen_project <- function() {
  
  # create null project
  project <- list(sim_parameters = NULL,
                  sim_output = NULL)
  class(project) <- "simplegen_project"
  
  # use default epi parameters
  project <- define_epi_parameters(project)
  
  # use default deme parameters
  project <- define_deme_parameters(project)
  
  # use default life table
  project <- define_demograpy(project)
  
  # use default migration parameters
  project <- define_migration(project)
  
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
