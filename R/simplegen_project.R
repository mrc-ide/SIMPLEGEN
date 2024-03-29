
#------------------------------------------------
#' @title Define new SIMPLEGEN project
#'
#' @description SIMPLEGEN works using projects, which are essentially just lists
#'   containing all the inputs and outputs of an analysis. Often we will pass
#'   the project in as both input and output to allow modification, as
#'   demonstrated in the examples below.
#'
#' @export
#' @examples
#' # define a new project and modify this project to hold parameters of the
#' # epidemiological model
#' myproj <- simplegen_project()

simplegen_project <- function() {
  
  # create empty project
  project <- list(epi_model_parameters = NULL,
                  epi_sampling_parameters = NULL,
                  epi_output = NULL,
                  sample_details = NULL,
                  pruned_record = NULL,
                  genetic_parameters = NULL,
                  haplotype_tree = NULL,
                  block_tree = NULL,
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
    
    message("")
  }
  
  # print epi sampling parameters
  if (!is.null(p$epi_sampling_parameters)) {
    message("Sampling strategy:")
    
    df_daily <- p$epi_sampling_parameters$daily
    n_daily_outputs <- ifelse(is.null(df_daily), 0, nrow(df_daily))
    message(sprintf("  daily outputs: %s", n_daily_outputs))
    
    df_sweeps <- p$epi_sampling_parameters$sweeps
    n_sweep_times <- ifelse(is.null(df_sweeps), 0, length(unique(df_sweeps$time)))
    message(sprintf("  sweep timepoints: %s", n_sweep_times))
    
    df_surveys <- p$epi_sampling_parameters$surveys
    n_survey_outputs <- ifelse(is.null(df_surveys), 0, nrow(df_surveys))
    message(sprintf("  surveys: %s", n_survey_outputs))
    
    message("")
  }
  
  # print summary of epi output
  if (!is.null(p$epi_output)) {
    message("Output:")
    max_time <- max(c(p$epi_output$daily$time, p$epi_output$sweeps$time))
    message(sprintf("  simulation days: %s", max_time))
    
    message("")
  }
  
  # print genetic parameters
  if (!is.null(p$genetic_parameters)) {
    message("Genetic model:")
    
    alpha <- p$genetic_parameters$alpha
    oo_dist <- p$genetic_parameters$oocyst_distribution
    mean_oocysts <- sum(seq_along(oo_dist)*oo_dist) / sum(oo_dist)
    
    hep_dist <- p$genetic_parameters$hepatocyte_distribution
    mean_hepatocytes <- sum(seq_along(hep_dist)*hep_dist) / sum(hep_dist)
    
    message(sprintf("  alpha: %s", signif(alpha, digits = 2)))
    message(sprintf("  mean oocysts: %s", signif(mean_oocysts, digits = 2)))
    message(sprintf("  mean hepatocytes: %s", signif(mean_hepatocytes, digits = 2)))
    message(sprintf("  recombination rate: %s", p$genetic_parameters$r))
    message(sprintf("  contigs: %s", length(p$genetic_parameters$contig_lengths)))
  }
  
  invisible(object)
}

