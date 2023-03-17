
#------------------------------------------------
#' @title Check that SIMPLEGEN package has loaded successfully
#'
#' @description Simple function to check that SIMPLEGEN package has loaded
#'   successfully. Prints "SIMPLEGEN loaded successfully!".
#'
#' @export

check_SIMPLEGEN_loaded <- function() {
  message("SIMPLEGEN loaded successfully!")
}

#------------------------------------------------
# perform checks on epi model parameters. NB. check_epi_model_params() is run
# before process_epi_model_params()
#' @noRd
check_epi_model_params <- function(project) {
  
  # check project class
  assert_class(project, "simplegen_project")
  
  # check that epi model parameters exist
  assert_non_null(project$epi_model_parameters, message = "no epi model parameters defined. See ?define_epi_model_parameters")
  
  # extract model parameters
  x <- project$epi_model_parameters
  
  # check basic parameter formats and ranges
  assert_single_bounded(x$a, name = "a")
  assert_single_bounded(x$p, name = "p")
  assert_single_pos(x$mu, name = "mu")
  
  # check delay parameter formats and ranges
  assert_single_pos_int(x$u, zero_allowed = FALSE, name = "u")
  assert_single_pos_int(x$v, zero_allowed = FALSE, name = "v")
  assert_single_pos_int(x$g, zero_allowed = FALSE, name = "g")
  
  # check transition probability formats and ranges
  assert_vector_bounded(x$prob_infection, name = "prob_infection")
  assert_vector_bounded(x$prob_acute, name = "prob_acute")
  assert_vector_bounded(x$prob_AC, name = "prob_AC")
  
  # check duration distribution formats and ranges. NB. these mapply() calls
  # work irrespective of whether distributions are defined as vectors or lists
  # over vectors
  mapply(assert_pos, x$duration_acute, name = "duration_acute")
  mapply(assert_pos, x$duration_chronic, name = "duration_chronic")
  mapply(assert_pos, x$time_treatment_acute, name = "time_treatment_acute")
  mapply(assert_pos, x$time_treatment_chronic, name = "time_treatment_chronic")
  mapply(assert_pos, x$duration_prophylactic, name = "duration_prophylactic")
  
  # check daily probabilities formats and ranges
  mapply(assert_bounded, x$detectability_microscopy_acute, name = "detectability_microscopy_acute")
  mapply(assert_bounded, x$detectability_microscopy_chronic, name = "detectability_microscopy_chronic")
  mapply(assert_bounded, x$detectability_PCR_acute, name = "detectability_PCR_acute")
  mapply(assert_bounded, x$detectability_PCR_chronic, name = "detectability_PCR_chronic")
  mapply(assert_bounded, x$infectivity_acute, name = "infectivity_acute")
  mapply(assert_bounded, x$infectivity_chronic, name = "infectivity_chronic")
  
  # check treatment seekeing parameters, including check that mean and standard
  # deviation are possible given constraints of Beta distribution
  assert_single_bounded(x$treatment_seeking_mean, name = "treatment_seeking_mean")
  assert_single_pos(x$treatment_seeking_sd, zero_allowed = TRUE, name = "treatment_seeking_sd")
  if ((x$treatment_seeking_mean > 0) & (x$treatment_seeking_mean < 1)) {
    sd_max <- sqrt(x$treatment_seeking_mean*(1 - x$treatment_seeking_mean))
    error_message <- sprintf(paste0("treatment_seeking_sd must be less than",
                                    " sqrt(treatment_seeking_mean*(1 - treatment_seeking_mean)),",
                                    " in this case %s, otherwise the Beta distribution is undefined"), sd_max)
    assert_le(x$treatment_seeking_sd, sd_max, message = error_message)
  }
  
  # check migration matrix formats and ranges
  assert_square_matrix(x$mig_mat, name = "mig_mat")
  assert_bounded(x$mig_mat, name = "mig_mat")
  
  # check that deme vector length equals migration matrix dimensions, or
  # alternatively are length 1 in which case the same value is applied over all
  # demes
  n_demes <- nrow(x$mig_mat)
  msg <- sprintf(paste0("H, M, and seed_infections",
                        " must be vectors with length equal to the number of rows",
                        " in the migration matrix (%s), or alternatively vectors of",
                        " length 1 in which case the same values are used over all demes"), n_demes)
  if (length(x$H) != 1) {
    assert_length(x$H, n_demes, message = msg)
  }
  if (length(x$M) != 1) {
    assert_length(x$M, n_demes, message = msg)
  }
  if (length(x$seed_infections) != 1) {
    assert_length(x$seed_infections, n_demes, message = msg)
  }
  
  # check deme properties, formats and ranges
  assert_vector_pos_int(x$H, name = "H")
  assert_vector_pos_int(x$seed_infections, name = "seed_infections")
  assert_vector_pos_int(x$M, name = "M")
  mapply(assert_leq, x$seed_infections, x$H, name_x = "seed_infections", name_y = "H")
  
  # check life table formats and ranges
  assert_bounded(x$life_table, name = "life_table")
  assert_eq(x$life_table[length(x$life_table)], 1, message = "the final value in the life table must be 1 representing a 100%% chance of dying that year to ensure a closed population")
  
  # check misc formats and ranges
  assert_single_pos_int(x$max_inoculations, zero_allowed = FALSE, name = "max_infections")
}

