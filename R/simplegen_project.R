
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
summary.simplegen_project <- function(object, ...) {
  p <- object
  
  # if empty project
  if (all(mapply(is.null, p))) {
    message("(empty project)")
    invisible(object)
  }
  
  # print epi parameters
  if (!is.null(p$epi_parameters)) {
    message("Epidemiological model:")
    n_demes <- length(p$epi_parameters$H)
    message(sprintf("  demes: %s", n_demes))
    message(sprintf("  H:\t %s", paste(p$epi_parameters$H, collapse = ", ")))
    message(sprintf("  M:\t %s", paste(p$epi_parameters$M, collapse = ", ")))
    message(sprintf("  seed infections: %s", paste(p$epi_parameters$seed_infections, collapse = ", ")))
  }
  
  # print sampling strategy
  if (!is.null(p$sampling_strategy)) {
    message("Sampling strategy:")
    n_time <- unique(p$sampling_strategy$time)
    n_samp_time <- mapply(sum, split(p$sampling_strategy$n, f = p$sampling_strategy$time))
    if (length(n_time) <= 5) {
      message(sprintf("  time: %s", paste(n_time, collapse = ", ")))
      message(sprintf("  n: %s", paste(n_samp_time, collapse = ", ")))
    } else {
      message("  time: (more than 5)")
      message("  n: (more than 5)")
    }
    
  }
  
  invisible(object)
}

