
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
  
  # set all parameters to NULL
  param_names <- c("a", "p", "mu", "u", "v", "g",
                   "prob_infection", "prob_acute", "prob_AC",
                   "duration_acute", "duration_chronic",
                   "infectivity_acute", "infectivity_chronic",
                   "max_innoculations",
                   "H", "seed_infections", "M",
                   "life_table")
  project$sim_parameters <- replicate(length(param_names), NULL)
  names(project$sim_parameters) <- param_names
  
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

