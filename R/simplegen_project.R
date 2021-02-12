
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
  project <- list(epi_model_parameters = NULL,
                  epi_sampling_parameters = NULL,
                  epi_output = NULL,
                  genetic_parameters = NULL,
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
  
  # print epi model parameters
  if (!is.null(p$epi_model_parameters)) {
    message("Epidemiological model:")
    n_demes <- length(p$epi_model_parameters$H)
    message(sprintf("  demes: %s", n_demes))
    message(sprintf("  H:\t %s", paste(p$epi_model_parameters$H, collapse = ", ")))
    message(sprintf("  M:\t %s", paste(p$epi_model_parameters$M, collapse = ", ")))
    message(sprintf("  seed infections: %s", paste(p$epi_model_parameters$seed_infections, collapse = ", ")))
  }
  
  # print epi sampling parameters
  #if (!is.null(p$epi_sampling_parameters)) {
  #  message("Sampling strategy:")
  #  n_time <- unique(p$epi_sampling_parameters$time)
  #  n_samp_time <- mapply(sum, split(p$epi_sampling_parameters$n, f = p$epi_sampling_parameters$time))
  #  if (length(n_time) <= 5) {
  #    message(sprintf("  time: %s", paste(n_time, collapse = ", ")))
  #    message(sprintf("  n: %s", paste(n_samp_time, collapse = ", ")))
  #  } else {
  #    message("  time: (more than 5)")
  #    message("  n: (more than 5)")
  #  }
  #}
  
  # print summary of epi output
  #if (!is.null(p$sample_output)) {
  #  message("Sample output:")
  #  N <- nrow(p$sample_output)
  #  n_pos <- sum(p$sample_output$positive)
  #  message(sprintf("  prevalence: %s/%s (%s%%)", n_pos, N, signif(n_pos/N*100, digits = 2)))
  #}
  
  # print genetic parameters
  if (!is.null(p$genetic_parameters)) {
    message("Genetic model:")
    
    oo_dist <- p$genetic_parameters$oocyst_distribution
    mean_oocysts <- sum(seq_along(oo_dist)*oo_dist) / sum(oo_dist)
    
    hep_dist <- p$genetic_parameters$hepatocyte_distribution
    mean_hepatocytes <- sum(seq_along(hep_dist)*hep_dist) / sum(hep_dist)
    
    message(sprintf("  recombination rate: %s", p$genetic_parameters$r))
    message(sprintf("  mean oocysts: %s", signif(mean_oocysts, digits = 2)))
    message(sprintf("  mean hepatocytes: %s", signif(mean_hepatocytes, digits = 2)))
    message(sprintf("  contigs: %s", length(p$genetic_parameters$contig_lengths)))
  }
  
  invisible(object)
}

