
#------------------------------------------------
#' @title Check that SIMPLEGEN package has loaded successfully
#'
#' @description Simple function to check that SIMPLEGEN package has loaded
#'   successfully. Prints "SIMPLEGEN loaded successfully!" if so.
#'
#' @export

check_SIMPLEGEN_loaded <- function() {
  message("SIMPLEGEN loaded successfully!")
}

#------------------------------------------------
# perform checks on epi model parameters. NB.
# check_epi_model_params() is run before process_epi_model_params()
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
  if (x$treatment_seeking_mean > 0 & x$treatment_seeking_mean < 1) {
    sd_max <- sqrt(x$treatment_seeking_mean*(1 - x$treatment_seeking_mean))
    error_message <- sprintf(paste0("treatment_seeking_sd must be less than",
                                    " sqrt(treatment_seeking_mean*(1 - treatment_seeking_mean)),",
                                    " in this case %s, otherwise the Beta distribution is undefined"), sd_max)
    assert_le(x$treatment_seeking_sd, sd_max, message = error_message)
  }
  
  # check migration matrix formats and ranges
  assert_square_matrix(x$mig_mat, name = "mig_mat")
  assert_bounded(x$mig_mat, name = "mig_mat")

  assert_bounded(x$migration_probability, name = "migration_probability")
  
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
  
  # check deme properties formats and ranges
  assert_vector_pos_int(x$H, name = "H")
  assert_vector_pos_int(x$seed_infections, name = "seed_infections")
  assert_vector_pos_int(x$M, name = "M")
  mapply(assert_leq, x$seed_infections, x$H, name_x = "seed_infections", name_y = "H")
  
  # check life table formats and ranges
  assert_bounded(x$life_table, name = "life_table")
  assert_eq(x$life_table[length(x$life_table)], 1, message = "the final value in the life table must be 1 representing a 100%% chance of dying that year to ensure a closed population")
  
  # check misc formats and ranges
  assert_single_pos_int(x$max_infections, zero_allowed = FALSE, name = "max_infections")
}

#------------------------------------------------
# perform checks on epi sampling parameters
#' @noRd
check_epi_sampling_params <- function(project) {
  
  # check project class
  assert_class(project, "simplegen_project")
  
  # check that epi sampling parameters exist
  assert_non_null(project$epi_sampling_parameters, message = "no epi sampling parameters defined. See ?define_epi_sampling_parameters")
  
  # check individual elements of sampling output
  check_epi_sampling_params_daily(project$epi_sampling_parameters$daily)
  check_epi_sampling_params_sweeps(project$epi_sampling_parameters$sweeps)
  check_epi_sampling_params_surveys(project$epi_sampling_parameters$surveys)
  
}

#------------------------------------------------
# perform checks on daily sampling parameters
#' @noRd
check_epi_sampling_params_daily <- function(x) {
  
  # avoid "no visible binding" warning
  measure <- state <- NULL
  
  # return if null
  if (is.null(x)) {
    return()
  }
  
  # check dataframe column names
  assert_dataframe(x, message = "daily sampling parameters must be a dataframe")
  col_titles <- c("deme", "state", "measure", "diagnostic", "age_min", "age_max")
  assert_in(col_titles, names(x), message = sprintf("daily sampling parameters dataframe must contain the following columns: {%s}",
                                                    paste0(col_titles, collapse = ", ") ))
  
  # check deme format
  deme_mssg <- "deme must be a positive integer or -1"
  assert_vector_int(x$deme, message = deme_mssg)
  assert_greq(x$deme, -1, message = deme_mssg)
  assert_greq(x$deme[x$deme != -1], 1, message = deme_mssg)
  
  # check state format
  all_states <- c("A", "C", "S", "E", "P", "H", "Sv", "Ev", "Iv", "M")
  assert_in(x$state, all_states, message = sprintf("state must be one of: {%s}", paste0(all_states, collapse = ", ")))
  
  # states A and C
  if (any(x$state %in% c("A", "C"))) {
    df_sub <- subset(x, state %in% c("A", "C"))
    
    # check measure format
    measure_levels <- c("count", "prevalence", "incidence")
    assert_in(df_sub$measure, measure_levels, message = sprintf("for states A and C, measure must be one of: {%s}", paste0(measure_levels, collapse = ", ") ))
    
    # check diagnostic format
    diagnostic_levels <- c("true", "microscopy", "PCR")
    assert_in(df_sub$diagnostic, diagnostic_levels, message = sprintf("for states A and C, diagnostic must be one of: {%s}", paste0(diagnostic_levels, collapse = ", ") ))
    
    # check age format
    assert_non_NA(df_sub$age_min, message = "for all human states, age_min must be a positive integer or zero")
    assert_non_NA(df_sub$age_max, message = "for all human states, age_max must be a positive integer or zero")
    assert_pos_int(df_sub$age_min, zero_allowed = TRUE, message = "for all human states, age_min must be a positive integer or zero")
    assert_pos_int(df_sub$age_max, zero_allowed = TRUE, message = "for all human states, age_max must be a positive integer or zero")
    assert_greq(df_sub$age_max, df_sub$age_min, message = "for all human states, age_max must be greater than or equal to age_min")
  }
  
  # states S, E and P
  if (any(x$state %in% c("S", "E", "P"))) {
    df_sub <- subset(x, state %in% c("S", "E", "P"))
    
    # check measure format
    measure_levels <- c("count", "prevalence", "incidence")
    assert_in(df_sub$measure, measure_levels, message = sprintf("for states S, E and P, measure must be one of: {%s}", paste0(measure_levels, collapse = ", ") ))
    
    # check diagnostic format
    assert_NA(df_sub$diagnostic, message = "for states S, E and P, diagnostic must be NA")
    
    # check age format
    assert_non_NA(df_sub$age_min, message = "for all human states, age_min must be a positive integer or zero")
    assert_non_NA(df_sub$age_max, message = "for all human states, age_max must be a positive integer or zero")
    assert_pos_int(df_sub$age_min, zero_allowed = TRUE, message = "for all human states, age_min must be a positive integer or zero")
    assert_pos_int(df_sub$age_max, zero_allowed = TRUE, message = "for all human states, age_max must be a positive integer or zero")
    assert_greq(df_sub$age_max, df_sub$age_min, message = "for all human states, age_max must be greater than or equal to age_min")
  }
  
  # state H
  if (any(x$state == "H")) {
    df_sub <- subset(x, state == "H")
    
    # check measure format
    measure_levels <- c("count", "proportion", "EIR")
    assert_in(df_sub$measure, measure_levels, message = sprintf("for state H, measure must be one of: {%s}", paste0(measure_levels, collapse = ", ") ))
    
    # check diagnostic format
    assert_NA(df_sub$diagnostic, message = "for state H, diagnostic must be NA")
    
    # check age format
    assert_non_NA(df_sub$age_min, message = "for all human states, age_min must be a positive integer or zero")
    assert_non_NA(df_sub$age_max, message = "for all human states, age_max must be a positive integer or zero")
    assert_pos_int(df_sub$age_min, zero_allowed = TRUE, message = "for all human states, age_min must be a positive integer or zero")
    assert_pos_int(df_sub$age_max, zero_allowed = TRUE, message = "for all human states, age_max must be a positive integer or zero")
    assert_greq(df_sub$age_max, df_sub$age_min, message = "for all human states, age_max must be greater than or equal to age_min")
  }
  
  # states Sv, Ev, Iv
  if (any(x$state %in% c("Sv", "Ev", "Iv"))) {
    df_sub <- subset(x, state %in% c("Sv", "Ev", "Iv"))
    
    # check measure format
    measure_levels <- c("count", "prevalence")
    assert_in(df_sub$measure, measure_levels, message = sprintf("for states Sv, Ev and Iv, measure must be one of: {%s}", paste0(measure_levels, collapse = ", ") ))
    
    # check diagnostic format
    assert_NA(df_sub$diagnostic, message = "for states Sv, Ev and Iv, diagnostic must be NA")
    
    # check age format
    assert_NA(df_sub$age_min, message = "for all mosquito states, age_min must be NA")
    assert_NA(df_sub$age_max, message = "for all mosquito states, age_max must be NA")
    
  }
  
  # state M
  if (any(x$state == "M")) {
    df_sub <- subset(x, state == "M")
    
    # check measure format
    measure_levels <- c("count")
    assert_in(df_sub$measure, measure_levels, message = sprintf("for state M, measure must be one of: {%s}", paste0(measure_levels, collapse = ", ") ))
    
    # check diagnostic format
    assert_NA(df_sub$diagnostic, message = "for state M, diagnostic must be NA")
    
    # check age format
    assert_NA(df_sub$age_min, message = "for all mosquito states, age_min must be NA")
    assert_NA(df_sub$age_max, message = "for all mosquito states, age_max must be NA")
    
  }
  
}

#------------------------------------------------
# perform checks on population sweep sampling parameters
#' @noRd
check_epi_sampling_params_sweeps <- function(x) {
  
  # avoid "no visible binding" warning
  time <- NULL
  
  # return if null
  if (is.null(x)) {
    return()
  }
  
  # check dataframe column names
  assert_dataframe(x, message = "sweep sampling parameters must be a dataframe")
  col_titles <- c("time", "deme", "state", "measure", "diagnostic", "age_min", "age_max")
  assert_in(col_titles, names(x), message = sprintf("sweep sampling parameters dataframe must contain the following columns: {%s}",
                                                    paste0(col_titles, collapse = ", ") ))
  
  # check time format
  assert_vector_pos_int(x$time, zero_allowed = FALSE, message = "time must be a positive integer")
  
  # remaining columns should have identical format to daily dataframe
  check_epi_sampling_params_daily(subset(x, select = -time))
  
}

#------------------------------------------------
# perform checks on survey sampling parameters
#' @noRd
check_epi_sampling_params_surveys <- function(x) {
  
  # return if null
  if (is.null(x)) {
    return()
  }
  
  # check dataframe column names
  assert_dataframe(x, message = "survey sampling parameters must be a dataframe")
  col_titles <- c("t_start", "t_end", "reporting_interval", "deme", "measure", "sampling", "diagnostic", "age_min", "age_max", "sample_size")
  assert_in(col_titles, names(x), message = sprintf("survey sampling parameters dataframe must contain the following columns: {%s}",
                                                    paste0(col_titles, collapse = ", ") ))
  
  # check times
  assert_vector_pos_int(x$t_start, zero_allowed = FALSE, name = "t_start")
  assert_vector_pos_int(x$t_end, zero_allowed = FALSE, name = "t_end")
  assert_greq(x$t_end, x$t_start, name_x = "t_end", name_y = "t_start")
  
  # check reporting interval
  assert_vector_pos_int(x$reporting_interval, zero_allowed = FALSE, name = "reporting_interval")
  
  # check deme
  deme_mssg <- "deme must be a positive integer or -1"
  assert_vector_int(x$deme, message = deme_mssg)
  assert_greq(x$deme, -1, message = deme_mssg)
  assert_greq(x$deme[x$deme != -1], 1, message = deme_mssg)
  
  # check measure
  assert_in(x$measure, c("prevalence", "incidence"), name_x = "measure")
  
  # check sampling
  x %>%
    dplyr::filter(.data$measure == "prevalence") %>%
    dplyr::pull(.data$sampling) %>%
    assert_NA(message = "survey sampling method for prevalence output must be NA")
  x %>%
    dplyr::filter(.data$measure == "incidence") %>%
    dplyr::pull(.data$sampling) %>%
    assert_in(c("ACD", "PCD"), message = "survey sampling method for incidence output must be ACD or PCD")
  
  # check diagnostic
  assert_in(x$diagnostic, c("true", "microscopy", "PCR"), name_x = "diagnostic")
  
  # check age range
  assert_vector_pos_int(x$age_min, zero_allowed = TRUE, name = "age_min")
  assert_vector_pos_int(x$age_max, zero_allowed = TRUE, name = "age_max")
  assert_greq(x$age_max, x$age_min, name_x = "age_max", name_y = "age_min")
  
  # check sample size
  assert_vector_pos(x$sample_size, zero_allowed = FALSE, name = "sample_size")
}

#------------------------------------------------
# check that at least one of daily, sweeps or survey output is specified
#' @noRd
check_epi_sampling_params_present <- function(project) {
  
  # check that epi sampling parameters contain correct elements
  assert_eq(names(project$epi_sampling_parameters), c("daily", "sweeps", "surveys"))
  
  # check that at least one element is non-null
  if (is.null(project$epi_sampling_parameters$daily) &
      is.null(project$epi_sampling_parameters$sweeps) &
      is.null(project$epi_sampling_parameters$surveys)) {
    stop("at least one output from the epidemiological model must be specified. See ?define_epi_sampling_parameters")
  }
  
}

#------------------------------------------------
# perform checks on genetic parameters
#' @noRd
check_genetic_params <- function(x) {
  
  # perform checks
  assert_vector_pos(x$oocyst_distribution)
  assert_vector_pos(x$hepatocyte_distribution)
  assert_single_pos(x$alpha, zero_allowed = FALSE)
  assert_single_pos(x$r, zero_allowed = TRUE)
  assert_vector_pos_int(x$contig_lengths, zero_allowed = FALSE)
  
}
