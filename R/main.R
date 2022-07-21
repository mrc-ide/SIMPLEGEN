
#------------------------------------------------
#' @title Load system file
#'
#' @description Load a file from within the inst/extdata folder of the SIMPLEGEN
#'   package. File extension must be one of .csv, .txt, or .rds.
#'
#' @param name the name of a file within the inst/extdata folder.
#'
#' @importFrom data.table fread
#' @export

simplegen_file <- function(name) {
  
  # check that valid file extension
  ext <- strsplit(name, "\\.")[[1]]
  ext <- ext[length(ext)]
  assert_in(ext, c("txt", "csv", "rds"), message = "file extension not valid")
  
  # get full file path
  name_full <- system.file("extdata/", name, package = 'SIMPLEGEN', mustWork = TRUE)
  
  # read in file
  if (ext == "rds") {
    ret <- readRDS(name_full)
  } else {
    ret <- data.table::fread(name_full, data.table = FALSE)
  }
  
  return(ret)
}

#------------------------------------------------
#' @title A life table taken from Mali
#'
#' @description The default life table used within the inbuilt SIMPLEGEN
#'   transmission model. Values represent probabilities of death in one-year age
#'   bands.
#'
#' @export

life_table_Mali <- function() {
  life_table_raw <- simplegen_file("Mali_life_table.csv")
  ret <- life_table_raw$prop_death
  return(ret)
}

#------------------------------------------------
#' @title Define parameters of transmission simulation model
#'
#' @description Define the parameters that will be used when simulating data
#'  from the inbuilt SIMPLEGEN transmission model. Any parameters that are
#'  defined as \code{NULL} will be set to default values.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. \code{mu = -log(p)} unless
#'   otherwise specified.
#' @param u intrinsic incubation period. The number of days from infection to
#'   blood-stage infection in a human host.
#' @param v extrinsic incubation period. The number of days from infection to
#'   becoming infectious in a mosquito.
#' @param g lag time between human blood-stage infection and production of
#'   gametocytes.
#' @param prob_infection probability a human becomes infected after being bitten
#'   by an infected mosquito. If a vector then each subsequent value applies to
#'   a subsequent infection.
#' @param prob_acute probability an infection passes to the acute stage, rather
#'   than going directly to the chronic stage. If a vector then each subsequent
#'   value applies to a subsequent infection.
#' @param prob_AC probability an infection in the acute stage passes to the
#'   chronic stage prior to recovery, rather than recovering directly. If a
#'   vector then each subsequent value applies to a subsequent infection.
#' @param duration_acute,duration_chronic specifies the probability distribution
#'   of duration (in days) of acute or chronic disease. Can be a vector or a
#'   list. If a list then the first element specifies the probability
#'   distribution for the first incident of disease, the second element
#'   for the second incident and so on (the final element is used for all
#'   remaining incidents). If a vector then the same probability distribution is
#'   used for all incidents of disease. Values are automatically
#'   normalised to a proper probability mass distribution internally (i.e. a
#'   distribution that sums to 1).
#' @param detectability_microscopy_acute,detectability_microscopy_chronic
#'   specifies the probability of parasite detection via microscopy at each day
#'   since entering the acute or chronic stages. Can be a vector or a list. If a
#'   list then the first element specifies the probability distribution for the
#'   first incident of disease, the second element for the second incident and
#'   so on (the final element is used for all remaining incidents). If a vector
#'   then the same probability distribution is used for all incidents of
#'   disease. Values are automatically normalised to a proper probability mass
#'   distribution internally (i.e. a distribution that sums to 1).
#' @param detectability_PCR_acute,detectability_PCR_chronic equivalent to
#'   `detectability_microscopy_acute` and `detectability_microscopy_chronic`,
#'   but for detection via polymerase chain reaction (PCR).
#' @param time_treatment_acute,time_treatment_chronic the probability
#'   distribution of time (in days since entering the stage) until considering
#'   seeking treatment for acute or chronic infection. If the infection clears
#'   naturally prior to this time then the host does not seek treatment,
#'   otherwise the host seeks treatment according to their host-specific
#'   treatment seeking parameter (see \code{treatment_seeking_alpha} and
#'   \code{treatment_seeking_beta}). Can be a vector or a list. If a list then
#'   the first element specifies the probability distribution for the first
#'   incident of disease, the second element for the second incident and so on
#'   (the final element is used for all remaining incidents). If a vector then
#'   the same probability distribution is used for all incidents of disease.
#'   Values are automatically normalised to a proper probability mass
#'   distribution internally (i.e. a distribution that sums to 1).
#' @param treatment_seeking_mean,treatment_seeking_sd treatment seeking is
#'   modelled in two stages; first, the time at which hosts \emph{consider}
#'   seeking treatment is drawn from the probability distribution in
#'   \code{time_treatment_acute} or \code{time_treatment_chronic}, and second,
#'   whether they go ahead with treatment-seeking at this point in time depends
#'   on their personal host-specific treatment-seeking probability. These latter
#'   probabilities are drawn independently for each host at birth from a Beta
#'   distribution with mean `treatment_seeking_mean` and standard deviation
#'   `treatment_seeking_sd`. Hence, these parameters describe the average
#'   propensity to seek treatment in the population, and the degree of
#'   inequality in treatment seeking. Note that `treatment_seeking_sd`` must be
#'   less than `sqrt(treatment_seeking_mean*(1 - treatment_seeking_mean))``,
#'   otherwise the Beta distribution is undefined.
#' @param duration_prophylactic vector defining the probability distribution of
#'   duration (in days) of prophylaxis following treatment.
#' @param infectivity_acute,infectivity_chronic the probability that a mosquito
#'   becomes infected upon biting an acutely or chronically infective human host
#'   at each day since entering the infective state. Can be a vector or a list.
#'   If a list then the first element specifies the infectivity for the first
#'   incident of acute disease, the second element for the second incident and
#'   so on (the final element is used for all remaining incidents). If a vector
#'   then the same infectivity is used for all episodes of acute disease.
#' @param max_infections maximum number of infections that an individual
#'   can carry simultaneously.
#' @param H vector specifying human population size in each deme.
#' @param seed_infections vector specifying the initial number of infected human
#'   hosts in each deme. Infected hosts are assumed to have just been bitten,
#'   meaning they have just entered the latent phase.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#' @param mig_mat migration matrix specifying the daily probability of human
#'   migration.
#' @param life_table vector specifying the probability of death of a human host
#'   in each one-year age group. The final value must be 1 to ensure a closed
#'   population. Defaults to a life table taken from Mali - see
#'   \code{?life_table_Mali} for details.
#'
#' @importFrom stats dgeom
#' @export

define_epi_model_parameters <- function(project,
                                        a = 0.3,
                                        p = 0.85,
                                        mu = -log(p),
                                        u = 12,
                                        v = 10,
                                        g = 12,
                                        prob_infection = 0.6,
                                        prob_acute = 1.0,
                                        prob_AC = 1.0,
                                        duration_acute = dgeom(1:25, 1/5),
                                        duration_chronic = dgeom(1:250, 1/50),
                                        detectability_microscopy_acute = 1,
                                        detectability_microscopy_chronic = 0.1,
                                        detectability_PCR_acute = 1,
                                        detectability_PCR_chronic = 1,
                                        time_treatment_acute = dgeom(1:100, 1/20),
                                        time_treatment_chronic = dgeom(1:100, 1/20),
                                        treatment_seeking_mean = 0.5,
                                        treatment_seeking_sd = 0.1,
                                        duration_prophylactic = dgeom(1:25, 1/5),
                                        infectivity_acute = 0.07,
                                        infectivity_chronic = 0.07,
                                        max_infections = 5,
                                        H = 1000,
                                        seed_infections = 100,
                                        M = 1000,
                                        mig_mat = diag(1),
                                        migration_probability = c(0.1),
                                        life_table = life_table_Mali()) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # get list of all input values, including those set by default
  all_args <- within(as.list(environment()), rm(project))
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values
  if (is.null(project$epi_model_parameters)) {
    project$epi_model_parameters <- all_args
  }
  
  # get list of only input values defined by user
  user_arg_names <- names(as.list(match.call()))
  user_arg_names <- setdiff(user_arg_names, c("", "project"))
  user_args <- all_args[user_arg_names]
  
  # overwrite parameters defined by user
  project$epi_model_parameters[user_arg_names] <- user_args
  
  # perform checks on parameters
  check_epi_model_params(project)
  
  # standardise and process parameters
  project <- process_epi_model_params(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
# convert epi model parameters to standardised types. NB.
# check_epi_model_params() is run before process_epi_model_params()
#' @noRd
process_epi_model_params <- function(project) {
  
  # function that forces objects that can be defined as list or vector into a
  # list
  force_list <- function(x) {
    if (!is.list(x)) {
      x <- list(x)
    }
    return(x)
  }
  
  # force duration distributions and daily probabilities to list
  name_vec <- c("duration_acute", "duration_chronic",
                "time_treatment_acute", "time_treatment_chronic",
                "duration_prophylactic",
                "detectability_microscopy_acute", "detectability_microscopy_chronic",
                "detectability_PCR_acute", "detectability_PCR_chronic",
                "infectivity_acute", "infectivity_chronic")
  project$epi_model_parameters[name_vec] <- mapply(force_list, project$epi_model_parameters[name_vec], SIMPLIFY = FALSE)
  
  # function that forces vectors to a given length by replicating scalar values
  force_veclength <- function(x, n) {
    if (length(x) == 1) {
      x <- rep(x, n)
    }
    return(x)
  }
  
  # force deme properties to vectors over demes
  name_vec <- c("H", "M", "seed_infections")
  n_demes <- nrow(project$epi_model_parameters$mig_mat)
  project$epi_model_parameters[name_vec] <- mapply(force_veclength, project$epi_model_parameters[name_vec],
                                                   MoreArgs = list(n = n_demes), SIMPLIFY = FALSE)
  
  # return
  invisible(project)
}

#------------------------------------------------
# get stable demography distribution from life table
#' @noRd
get_demography <- function(life_table) {
  
  # check that life_table values are in the range [0,1]
  assert_bounded(life_table)
  
  # compute distribution of age of death
  n <- length(life_table)
  age_death <- rep(0, n)
  remaining <- 1
  for (i in 1:n) {
    age_death[i] <- remaining*life_table[i]
    remaining <- remaining*(1 - life_table[i])
  }
  
  # check that age_death is a proper probability distribution with no
  # probability mass escaping
  assert_eq(sum(age_death), 1, message = "life table does not result in a 100% probability of eventual death")
  
  # convert life table to transition matrix
  m <- matrix(0, n, n)
  m[col(m) == (row(m)+1)] <- 1 - life_table[1:(n-1)]
  m[,1] <- 1 - rowSums(m)
  
  # convert to rates
  r = m - diag(n)
  
  # compute Eigenvalues of the rate matrix
  E = eigen(t(r))
  
  # there should be one Eigenvalue that is zero (up to limit of computational
  # precision). Find which Eigenvalue this is
  w <- which.min(abs(E$values))
  
  # the stable solution is the corresponding Eigenvector, suitably normalised
  age_stable <- Re(E$vectors[,w]/sum(E$vectors[,w]))
  
  # return list of distributions
  ret <- list(life_table = life_table,
              age_death = age_death,
              age_stable = age_stable)
  return(ret)
}

#------------------------------------------------
#' @title Define outputs produced from the transmission model
#'
#' @description Loads one or more dataframes into the SIMPLEGEN project that
#'   specify which outputs will be returned from the epidemiological model.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param daily a dataframe of daily outputs.
#' @param sweeps a dataframe of outputs at specific time points.
#' @param surveys a dataframe specifying surveys to be conducted on the
#'   population.
#'
#' @export

define_epi_sampling_parameters <- function(project,
                                           daily = NULL,
                                           sweeps = NULL,
                                           surveys = NULL) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # at least one input must be non-NULL
  if (is.null(daily) & is.null(sweeps) & is.null(surveys)) {
    stop("must define at least one output type")
  }
  
  # get list of all input values, including those set by default
  all_args <- within(as.list(environment()), rm(project))
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values
  if (is.null(project$epi_sampling_parameters)) {
    project$epi_sampling_parameters <- all_args
  }
  
  # get list of only input values defined by user
  user_arg_names <- names(as.list(match.call()))
  user_arg_names <- setdiff(user_arg_names, c("", "project"))
  user_args <- all_args[user_arg_names]
  
  # overwrite parameters defined by user
  project$epi_sampling_parameters[user_arg_names] <- user_args
  
  # perform checks on final parameters
  check_epi_sampling_params(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
# expand surveys dataframe to give precise number of samples needed on each day
# of study
#' @importFrom utils tail
#' @noRd
expand_surveys <- function(x) {
  
  # return if null
  if (is.null(x)) {
    return()
  }
  
  # apply over all rows of survey dataframe
  ret <- mapply(function(i) {
    
    # get sequence of breaks corresponding to reporting times
    break_vec <- seq(x$t_start[i] - 1, x$t_end[i], x$reporting_interval[i])
    if (tail(break_vec, 1) != x$t_end[i]) {
      break_vec <- c(break_vec, x$t_end[i])
    }
    
    # make dataframe with sample times and reporting times
    ret <- data.frame(study_ID = i,
                      sampling_time = x$t_start[i]:x$t_end[i]) %>%
      dplyr::mutate(reporting_time = cut(.data$sampling_time, breaks = break_vec)) %>%
      dplyr::mutate(reporting_time = break_vec[as.numeric(.data$reporting_time) + 1])
    
    return(ret)
  }, seq_len(nrow(x)), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
  
  # sort by increasing sampling time
  ret <- dplyr::arrange(ret, .data$sampling_time)
  
  return(ret)
}

#------------------------------------------------
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from the inbuilt SIMPLEGEN transmission model. Model
#'   parameters are taken from the \code{epi_model_parameters} slot of the
#'   project, and the sampling strategy is taken from the
#'   \code{epi_sampling_parameters} slot. Outputs are written to the
#'   \code{epi_output} slot.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param max_time run simulation for this many days.
#' @param save_transmission_record whether to write the transmission record to
#'   file.
#' @param transmission_record_location the file path that the transmission
#'   record will be written to.
#' @param overwrite_transmission_record if \code{TRUE} the transmission record
#'   will overwrite any existing file by the same name. \code{FALSE} by default.
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100% to avoid large amounts of output
#'   being printed to markdown files.
#' @param silent whether to suppress written messages to the console.
#'
#' @importFrom utils txtProgressBar
#' @export

sim_epi <- function(project,
                    max_time = 365,
                    save_transmission_record = FALSE,
                    transmission_record_location = "",
                    overwrite_transmission_record = FALSE,
                    pb_markdown = FALSE,
                    silent = FALSE) {
  
  
  # avoid "no visible binding" warning
  numer <- denom <- measure <- NULL
  
  # ---------- check inputs ----------
  
  assert_class(project, "simplegen_project")
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_single_logical(save_transmission_record)
  assert_string(transmission_record_location)
  assert_single_logical(overwrite_transmission_record)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # optionally return warning if will overwrite transmission record file
  if (save_transmission_record & !overwrite_transmission_record) {
    if (file.exists(transmission_record_location)) {
      stop(sprintf("file already exists at %s. Change target location, or use argument `overwrite_transmission_record = TRUE` to manually override this warning", transmission_record_location))
    }
  }
  
  # ensure that parameters are loaded and pass checks
  check_epi_model_params(project)
  check_epi_sampling_params(project)
  
  # ---------- define arguments  ----------
  
  # create argument list by concatenating project parameters
  args <- c(project$epi_model_parameters,
            project$epi_sampling_parameters)
  
  # append function arguments
  args <- c(args,
            list(max_time = max_time,
                 save_transmission_record = save_transmission_record,
                 transmission_record_location = transmission_record_location,
                 pb_markdown = pb_markdown,
                 silent = silent))
  
  # get migration matrix into list
  args$mig_mat <- matrix_to_rcpp(args$mig_mat)
  
  # get complete demography from life table
  demog <- get_demography(project$epi_model_parameters$life_table)
  args <- c(args, list(age_death = demog$age_death,
                       age_stable = demog$age_stable))
  
  # establish which outputs are required
  args$any_daily_outputs <- !is.null(args$daily)
  args$any_sweep_outputs <- !is.null(args$sweeps)
  args$any_survey_outputs <- !is.null(args$surveys)
  
  # replace "proportion" with "prevalence", as this uses the same calculation
  if (args$any_daily_outputs) {
    args$daily$measure <- replace(args$daily$measure, args$daily$measure == "proportion", "prevalence")
  }
  
  # get sampling strategy indices into 0-indexed numerical format
  sampling_to_cpp_format <- function(x) {
    if (is.null(x)) {
      return(x)
    }
    w <- which(x$deme != -1)
    x$deme[w] <- x$deme[w] - 1
    x$measure <- match(x$measure, c("count", "prevalence", "incidence", "EIR")) - 1
    if ("state" %in% names(x)) {
      x$state <- match(x$state, c("S", "E", "A", "C", "P", "H", "Sv", "Ev", "Iv", "M")) - 1
    }
    if ("diagnostic" %in% names(x)) {
      x$diagnostic <- match(x$diagnostic, c("true", "microscopy", "PCR")) - 1
    }
    if ("sampling" %in% names(x)) {
      x$sampling <- match(x$sampling, c(NA, "ACD", "PCD")) - 1
    }
    return(x)
  }
  args$daily <- sampling_to_cpp_format(args$daily)
  args$sweeps <- sampling_to_cpp_format(args$sweeps)
  surveys_raw <- args$surveys
  args$surveys <- sampling_to_cpp_format(args$surveys)
  
  # get unique times at which sweeps happen
  args$sweep_time_ordered <- NULL
  if (args$any_sweep_outputs) {
    args$sweep_time_ordered <- sort(unique(args$sweeps$time))
  }
  
  # add columns to surveys
  if (args$any_survey_outputs) {
    args$surveys$n_days <- args$surveys$t_end - args$surveys$t_start + 1
  }
  
  # get expanded version of surveys dataframe
  args$surveys_expanded <- expand_surveys(args$surveys)
  
  # functions
  args_functions <- list(update_progress = update_progress)
  
  # make progress bars
  pb_sim <- txtProgressBar(min = 0, max = max_time, initial = NA, style = 3)
  args_progress <- list(pb_sim = pb_sim)
  
  
  # ---------- run simulation ----------
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args, args_functions, args_progress)
  
  #return(output_raw)
  
  # ---------- process output ----------
  
  # wrangle daily output
  daily_output <- NULL
  if (args$any_daily_outputs) {
    
    # make output dataframe with raw values
    daily_df <- project$epi_sampling_parameters$daily
    daily_output <- cbind(time = rep(seq_len(max_time), each = nrow(daily_df)),
                          daily_df,
                          value = unlist(output_raw$daily_output),
                          row.names = NULL)
  }
  
  # wrangle sweep output
  sweeps_output <- NULL
  if (args$any_sweep_outputs) {
    
    # make output dataframe with raw values
    sweeps_output <- cbind(project$epi_sampling_parameters$sweeps,
                           value = unlist(output_raw$sweep_output),
                           row.names = NULL)
  }
  
  # wrangle surveys output
  sample_details <- surveys_output <- NULL
  if (args$any_survey_outputs) {
    
    # make individual-level output dataframe
    sample_details <- mapply(as.data.frame, output_raw$survey_output, SIMPLIFY = FALSE) %>%
      dplyr::bind_rows()
    sample_details$infection_IDs <- output_raw$survey_output_infection_IDs
    sample_details <- dplyr::arrange(sample_details, .data$study_ID)
    
    # get summary output from individual-level
    surveys_output <- get_survey_summary(sample_details,
                                         surveys_raw,
                                         args$surveys_expanded)
    
  }
  
  # append to project
  project$epi_output <- list(daily = daily_output,
                             sweeps = sweeps_output,
                             surveys = surveys_output)
  
  if (!is.null(sample_details)) {
    project$sample_details <- sample_details
  }
  
  # return
  invisible(project)
}

#------------------------------------------------
# get summary survey output from individual-level
#' @noRd
get_survey_summary <- function(surveys_indlevel,
                               surveys_raw,
                               surveys_expanded) {
  
  # get skeleton dataframe of all possible reporting days for each study
  df_skeleton <- surveys_expanded %>%
    dplyr::group_by(.data$study_ID, .data$reporting_time) %>%
    dplyr::summarise(reporting_interval = dplyr::n()) %>%
    dplyr::ungroup()
  
  # get counts over individual-level output
  df_counts <- surveys_indlevel %>%
    dplyr::group_by(.data$study_ID, .data$reporting_time) %>%
    dplyr::summarise(n_sampled = dplyr::n(),
                     n_true_pos = sum(.data$true_positive),
                     n_micro_pos = sum(.data$microscopy_positive),
                     n_PCR_pos = sum(.data$PCR_positive)) %>%
    dplyr::ungroup()
  
  # merge back with skeleton dataframe, get into long format over diagnostics and replace
  # NAs with 0s
  df_long <- df_counts %>%
    dplyr::full_join(df_skeleton, by = c("study_ID", "reporting_time")) %>%
    dplyr::arrange(.data$study_ID, .data$reporting_time) %>%
    tidyr::pivot_longer(cols = c(.data$n_true_pos, .data$n_micro_pos, .data$n_PCR_pos), names_to = "observed_diagnostic", values_to = "n_pos") %>%
    dplyr::mutate(n_sampled = ifelse(is.na(.data$n_sampled), 0, .data$n_sampled),
                  n_pos = ifelse(is.na(.data$n_pos), 0, .data$n_pos)) %>%
    dplyr::mutate(observed_diagnostic = match(.data$observed_diagnostic, c("n_true_pos", "n_micro_pos", "n_PCR_pos"))) %>%
    dplyr::mutate(observed_diagnostic = c("true", "microscopy", "PCR")[.data$observed_diagnostic])
  
  # merge back with original survey dataframe and subset to chosen diagnostic
  df_filtered <- surveys_raw %>%
    dplyr::mutate(study_ID = seq_along(.data$t_start)) %>%
    dplyr::select(.data$study_ID, .data$deme, .data$measure, .data$sampling, .data$diagnostic, .data$age_min, .data$age_max) %>%
    dplyr::left_join(df_long, by = "study_ID") %>%
    dplyr::filter(.data$observed_diagnostic == .data$diagnostic) %>%
    dplyr::select(-.data$observed_diagnostic) %>%
    dplyr::relocate(c("study_ID", "reporting_time", "reporting_interval", "deme", "measure", "sampling", "diagnostic", "age_min", "age_max", 
                      "n_pos", "n_sampled"))
  
  # finalise incidence and prevalence units
  df_final <- df_filtered %>%
    dplyr::rename(numerator = .data$n_pos,
                  denominator = .data$n_sampled) %>%
    dplyr::mutate(denominator = ifelse(.data$measure == "incidence", NA, .data$denominator)) %>%
    dplyr::mutate(value = ifelse(.data$measure == "prevalence",
                                 .data$numerator / .data$denominator * 100,
                                 .data$numerator / .data$reporting_interval * 365))
  
  return(df_final)
}

#------------------------------------------------
#' @title Prune the transmission record
#'
#' @description Reads in a saved transmission record from file. Combines this
#'   with sampling information in the project to produce a pruned version of the
#'   transmission record that contains only the events relevant to the final
#'   sample. This pruned record is saved to the project.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param transmission_record_location the file path to a transmission record
#'   already written to file.
#' @param silent whether to suppress written messages to the console.
#' 
#' @importFrom  utils read.csv
#' @export

prune_transmission_record <- function(project,
                                      transmission_record_location = "",
                                      silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_string(transmission_record_location)
  assert_single_logical(silent)
  
  # check that project contains individual-level info
  sample_details <- project$sample_details
  indlevel_error <- "project must contain a dataframe in project$sample_details slot, and this dataframe must contain the columns 'study_ID' and 'infection_IDs'"
  assert_dataframe(sample_details, message = indlevel_error)
  assert_in(c("study_ID", "infection_IDs"), names(sample_details), message = indlevel_error)
  assert_vector_pos_int(sample_details$study_ID, name = "study_ID")
  assert_vector_pos_int(unlist(sample_details$infection_IDs), name = "infection_IDs")
  
  # check transmission record exists
  if (!file.exists(transmission_record_location)) {
    stop(sprintf("could not find file at %s", transmission_record_location))
  }
  
  # check that headers formatted correctly
  first_row <- read.csv(transmission_record_location, nrows = 1)
  required_names <- c("time", "event", "human_ID", "mosquito_ID", "child_infection_ID", "parent_infection_ID")
  assert_in(required_names, names(first_row), message = sprintf("transmission record must contain the following column headers: {%s}", paste(required_names, collapse = ", ")))
  
  # extract starting infection IDs from sample info
  infection_IDs <- unlist(project$sample_details$infection_IDs)
  if (length(infection_IDs) == 0) {
    stop("no infection_IDs found in sample details. Check that there is at least one malaria positive host")
  }
  
  # define argument list
  args <- list(transmission_record_location = transmission_record_location,
               infection_IDs = infection_IDs,
               silent = silent)
  
  # run efficient C++ code
  output_raw <- prune_transmission_record_cpp(args)
  
  # process output into dataframe. parent_infection_ID column needs to be dealt
  # with separately because each element can be a vector
  output_processed <- output_raw[-which(names(output_raw) == "parent_infection_ID")] %>%
    as.data.frame() %>%
    dplyr::mutate(parent_infection_ID = output_raw$parent_infection_ID, .after = "child_infection_ID")
  
  # reverse order so same as original transmission record
  output_processed <- output_processed[rev(seq_len(nrow(output_processed))),]
  row.names(output_processed) <- NULL
  
  # add to project
  project$pruned_record <- output_processed
  
  # return project
  invisible(project)
}

#------------------------------------------------
#' @title Define genetic parameters
#'
#' @description Defines all parameters relating to the genetic model. This
#'   includes parameters related to recombination etc. that are used when
#'   simulating relatedness between lineages, as well as parameters related to
#'   mutation and sequencing errors etc. that are used at a later stage when
#'   converting relatedness into actual genotypes.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param oocyst_distribution vector specifying the probability distribution of
#'   each number of oocysts within the mosquito midgut.
#' @param hepatocyte_distribution vector specifying the probability distribution
#'   of the number of infected hepatocytes in a human host. More broadly, this
#'   defines the number of independent draws from the oocyst products that make
#'   it into the host bloodstream upon a bite from an infectious mosquito.
#' @param alpha parameter dictating the skew of lineage densities. Small
#'   values of \code{alpha} create a large skew, and hence make it likely that
#'   an oocyst will be produced from the same parents.
#' @param r the rate of recombination. The expected number of base pairs in a
#'   single recombinant block is 1/r.
#' @param contig_lengths vector of lengths (in bp) of each contig.
#' 
#' @importFrom stats dpois
#' @export

define_genetic_params <- function(project,
                                  r = 1e-6,
                                  alpha = 1.0,
                                  oocyst_distribution = dpois(1:10, lambda = 2),
                                  hepatocyte_distribution = dpois(1:10, lambda = 5),
                                  contig_lengths = c(643292, 947102, 1060087, 1204112, 1343552,
                                                     1418244, 1501717, 1419563, 1541723, 1687655,
                                                     2038337, 2271478, 2895605, 3291871)) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # get list of all input values, including those set by default
  all_args <- within(as.list(environment()), rm(project))
  
  # get list of only input values defined by user
  user_arg_names <- names(as.list(match.call()))
  user_arg_names <- setdiff(user_arg_names, c("", "project"))
  user_args <- all_args[user_arg_names]
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values
  if (is.null(project$genetic_parameters)) {
    project$genetic_parameters <- all_args
  }
  
  # otherwise overwrite parameters defined by user
  project$genetic_parameters[user_arg_names] <- user_args
  
  # perform checks on parameters
  check_genetic_params(project$genetic_parameters)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Simulate haplotype tree
#'
#' @description Reads in the pruned transmission record from file and uses this,
#'   along with a specified genetic model, to simulate a tree relating
#'   haplotypes to one another.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

sim_haplotype_tree <- function(project,
                               silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_single_logical(silent)
  
  # check pruned record exists
  assert_non_null(project$pruned_record,
                  message = "no pruned transmission record found in project")
  
  # check sample details exist
  assert_non_null(project$sample_details,
                  message = "no sample details dataframe found in project")
  
  # check for defined genetic params
  assert_non_null(project$genetic_parameters,
                  message = "no genetic parameters defined. See ?define_genetic_params")
  
  # subset sample to those with infections
  sample_positive <- project$sample_details %>%
    dplyr::filter(.data$true_positive)
  
  # determine whether optional extra columns have been used
  defined_densities <- ("parent_infection_density" %in% names(project$pruned_record))
  defined_deme <- ("deme" %in% names(project$pruned_record))
  
  # define arguments
  args <- append(project$genetic_parameters,
                 list(sample_human_IDs = sample_positive$sample_ID,
                      sample_infection_IDs = sample_positive$infection_IDs,
                      pruned_record = project$pruned_record,
                      defined_densities = defined_densities,
                      defined_deme = defined_deme,
                      silent = silent))
  
  # run efficient C++ code
  output_raw <- sim_haplotype_tree_cpp(args)
  
  # get haplotype tree into dataframe by ading columns to pruned transmission
  # record
  haplotype_tree <- project$pruned_record[output_raw$record_row + 1, ] %>%
    dplyr::mutate(child_haplo_ID = output_raw$child_haplo_ID,
                  parent_haplo_ID = output_raw$parent_haplo_ID)
  row.names(haplotype_tree) <- NULL
  
  # append haplotype IDs and densities to sample_details dataframe
  sample_details <- project$sample_details
  w <- which(sample_details$true_positive)
  haplo_ID_list <- replicate(nrow(sample_details), integer())
  haplo_ID_list[w] <- output_raw$sample_haplo_IDs
  sample_details$haplo_IDs <- haplo_ID_list
  
  # save objects back into project
  project$haplotype_tree <- haplotype_tree
  project$sample_details <- sample_details
  
  # return project
  invisible(project)
}

#------------------------------------------------
#' @title Get the relatedness between a series of haplotypes
#'
#' @description Calculates the relatedness up to a defined number of generations
#'   between all pairs of input haplotypes.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param haplo_IDs a vector of haplotype IDs on which to calculate pairwise
#'   relatedness.
#' @param generations number of generations back to search.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

get_haplotype_relatedness <- function(project,
                                      haplo_IDs,
                                      generations = 3,
                                      silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_vector_pos_int(haplo_IDs, zero_allowed = FALSE)
  assert_greq(length(haplo_IDs), 2, message = "must input at least 2 haplo IDs")
  assert_single_pos_int(generations, zero_allowed = FALSE)
  assert_single_logical(silent)
  
  # check that project has a haplotype tree computed
  assert_non_null(project$haplotype_tree, message = "project has no haplotype_tree object")
  
  # check that all haplo_IDs are within the haplotype tree
  if (!all(haplo_IDs %in% project$haplotype_tree$child_haplo_ID)) {
    stop("haplo_IDs could not be found within haplotype_tree child_haplo_ID column")
  }
  
  # define arguments
  args <- list(time = project$haplotype_tree$time,
               child_haplo_ID = project$haplotype_tree$child_haplo_ID,
               parent_haplo_ID = project$haplotype_tree$parent_haplo_ID,
               target_haplo_IDs = haplo_IDs,
               generations = generations,
               silent = silent)
  
  # run efficient C++ function
  output_raw <- get_haplotype_relatedness_cpp(args)
  
  # process output
  n_haplo_IDs <- length(haplo_IDs)
  pair_df <- pairwise_long(n_haplo_IDs)
  ret <- data.frame(haplo_ID_1 = haplo_IDs[pair_df$x],
                    haplo_ID_2 = haplo_IDs[pair_df$y],
                    relatedness = output_raw$pairwise_relatedness)
  
  return(ret)
}

#------------------------------------------------
#' @title Get the relatedness between a series of samples
#'
#' @description TODO
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param sample_IDs a vector of sample IDs on which to calculate pairwise
#'   relatedness.
#' @param generations number of generations back to search.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

get_sample_relatedness <- function(project,
                                   sample_IDs,
                                   generations = 3,
                                   silent = FALSE) {

  # check inputs
  assert_class(project, "simplegen_project")
  assert_vector_pos_int(sample_IDs, zero_allowed = FALSE)
  assert_greq(length(sample_IDs), 2, message = "must input at least 2 sample IDs")
  assert_single_pos_int(generations, zero_allowed = FALSE)
  assert_single_logical(silent)
  
  # check that project has sample details loaded
  assert_dataframe(project$sample_details, message = "project has no sample_details object")
  
  # check that sample_IDs are within sample dataframe
  if (!all(sample_IDs %in% c(project$sample_details$sample_ID))) {
    stop("not all sample_IDs could be found within sample_details dataframe")
  }
  
  # subset to target samples and columns
  sample_target <- project$sample_details %>%
    dplyr::filter(.data$true_positive == TRUE) %>%
    dplyr::filter(.data$sample_ID %in% sample_IDs) %>%
    dplyr::select(.data$sample_ID, .data$haplo_IDs)
  
  # remove all samples that are not positive (contain at least one haplotype)
  n_haplos <- mapply(function(x) {
    length(unlist(x))
  }, sample_target$haplo_IDs)
  if (!all(n_haplos > 0)) {
    sample_target<-sample_target[n_haplos > 0,]
    warning("Warning: Not all sample_IDs were positive (contained at least one haplotype), non-positives removed")
  }
  
  
  if (!silent) {
    message("Calculating relatedness between samples")
  }
  t0 <- Sys.time()
  
  # get dataframe of all pairwise comparisons between rows
  df_pairwise <- pairwise_long(nrow(sample_target), include_diagonal = TRUE)
  
  # get dataframe of all pairwise comparisons between samples
  tmp <- sample_target %>%
    dplyr::rename(sample_ID_1 = .data$sample_ID,
                  haplo_ID_1 = .data$haplo_IDs)
  sample_pairs <- tmp[df_pairwise$x,] %>%
    dplyr::bind_cols(sample_target[df_pairwise$y,]) %>%
    dplyr::rename(sample_ID_2 = .data$sample_ID,
                  haplo_ID_2 = .data$haplo_IDs)
  
  # get similar dataframe but with haplotypes expanded over rows
  sample_haplo_pairs <- apply(sample_pairs, 1, function(x) {
    ret <- tidyr::expand_grid(haplo_ID_1 = unlist(x$haplo_ID_1),
                              haplo_ID_2 = unlist(x$haplo_ID_2)) %>%
      dplyr::mutate(sample_ID_1 = x$sample_ID_1,
                    sample_ID_2 = x$sample_ID_2, .before = 1)
    ret
  }) %>%
    dplyr::bind_rows()
  
  # calculate relatedness between all pairs of haplotypes. NB, this will not
  # compare haplotypes against themselves
  haplo_relatedness <- get_haplotype_relatedness(project = project,
                                                 haplo_IDs = unlist(sample_target$haplo_IDs),
                                                 generations = generations,
                                                 silent = TRUE)
  
  # get similar dataframe but at the sample level, with relatedness averaged
  # over all haplos
  sample_relatedness <- sample_haplo_pairs %>%
    dplyr::left_join(haplo_relatedness, by = c("haplo_ID_1", "haplo_ID_2")) %>%
    dplyr::filter(!is.na(.data$relatedness)) %>%
    dplyr::group_by(.data$sample_ID_1, .data$sample_ID_2) %>%
    dplyr::summarise(relatedness = mean(.data$relatedness), .groups = "keep") %>%
    dplyr::ungroup()
  
  # merge back with complete list of pairwise samples. For comparisons where a
  # monogenomic sample is compaired against itself, set relatedness to 1 by
  # convension
  ret <- sample_pairs %>%
    dplyr::select(-c(.data$haplo_ID_1, .data$haplo_ID_2)) %>%
    dplyr::full_join(sample_relatedness, by = c("sample_ID_1", "sample_ID_2")) %>%
    dplyr::mutate(relatedness = ifelse(is.na(.data$relatedness), 1.0, .data$relatedness))
  
  if (!silent) {
    t1 <- Sys.time()
    message(sprintf("completed in %s seconds", signif(t1 - t0, 3)))
  }
  
  return(ret)
}

#------------------------------------------------
#' @title Simulate recombinant block tree
#'
#' @description Starting with the haplotype tree stored within the project,
#'   simulates the relatedness of lengths of the genome from a recombination
#'   model.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

sim_block_tree <- function(project,
                           silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_single_logical(silent)
  
  # check sample exists and contains at least two haplotypes
  assert_non_null(project$sample_details, message = "no sample details found in project")
  target_haplos <- unlist(project$sample_details$haplo_IDs)
  assert_greq(length(target_haplos), 2, message = "sample must contain at least two haplotype IDs within malaria-positive individuals")
  
  # check haplotype tree exists
  assert_non_null(project$haplotype_tree, message = "no haplotype tree found in project")
  
  # check for defined genetic params
  assert_non_null(project$genetic_parameters,
                  message = "no genetic parameters defined. See ?define_genetic_params")
  
  # define genetic model arguments
  args <- list(r = project$genetic_parameters$r,
               contig_lengths = project$genetic_parameters$contig_lengths,
               silent = silent)
  
  # append haplotype tree to arguments
  args <- append(args,
                 list(time = project$haplotype_tree$time,
                      child_haplo_ID = project$haplotype_tree$child_haplo_ID,
                      parent_haplo_ID = project$haplotype_tree$parent_haplo_ID))
  
  # run efficient C++ code
  output_raw <- sim_block_tree_cpp(args)
  
  # process output
  output_processed <- as.data.frame(output_raw)
  
  # add to project
  project$block_tree <- output_processed
  
  # get coalescence between blocks
  haplo_coalescence <- get_haplotype_coalescence(project, haplo_IDs = target_haplos, generations = )
  
  # return project
  invisible(project)
}

#------------------------------------------------
#' @title Get genome-wide coalescence times between haplotypes
#'
#' @description TODO
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param haplo_IDs a vector of haplotype IDs on which to calculate coalescence
#'   times.
#' @param generations number of generations back to search.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

get_haplotype_coalescence <- function(project,
                                      haplo_IDs,
                                      generations = 3,
                                      silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_vector_pos_int(haplo_IDs, zero_allowed = FALSE)
  assert_greq(length(unique(haplo_IDs)), 2, message = "must input at least 2 different haplo IDs")
  assert_single_pos_int(generations, zero_allowed = FALSE)
  assert_single_logical(silent)
  
  # check that project has contig lengths defined
  assert_non_null(project$genetic_parameters$contig_lengths)
  
  # check that project has a block tree computed
  assert_non_null(project$block_tree, message = "project has no block_tree object")
  
  # check that all haplo_IDs are within the block tree
  if (!all(haplo_IDs %in% unique(project$block_tree$child_haplo_ID))) {
    stop("haplo_IDs could not be found within block_tree child_haplo_ID column")
  }
  
  # check that contig values are within range specified by genetic parameters
  if (!all(project$block_tree$contig <= length(project$genetic_parameters$contig_lengths))) {
    stop("block_tree contains contig values that exceed the number of contigs found in genetic_parameters")
  }
  
  # define arguments
  args <- list(time = project$block_tree$time,
               child_haplo_ID = project$block_tree$child_haplo_ID,
               contig = project$block_tree$contig - 1,
               contig_start = project$block_tree$start,
               contig_end = project$block_tree$end,
               parent_haplo_ID = project$block_tree$parent_haplo_ID,
               target_haplo_IDs = haplo_IDs,
               contig_lengths = project$genetic_parameters$contig_lengths,
               generations = generations,
               silent = silent)
  
  # run efficient C++ function
  output_raw <- get_haplotype_coalescence_cpp(args)
  
  # process output
  output_processed <- output_raw$coalescence_events %>%
    dplyr::bind_rows()
  
  if (nrow(output_processed) > 0) {
    output_processed <- output_processed %>%
      dplyr::mutate(l1 = ifelse(.data$lineage1 < .data$lineage2, .data$lineage1, .data$lineage2),
                    l2 = ifelse(.data$lineage1 < .data$lineage2, .data$lineage2, .data$lineage1),
                    .before = 1) %>%
      dplyr::select(-c(.data$lineage1, .data$lineage2)) %>%
      dplyr::rename(lineage1 = .data$l1, lineage2 = .data$l2)
  }
  
  return(output_processed)
}

#------------------------------------------------
#' @title Get genome-wide coalescence times between haplotypes
#'
#' @description TODO
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param haplo_IDs a vector of haplotype IDs on which to calculate coalescence
#'   times.
#' @param sampling_time time at which sample taken.
#' @param generations number of generations back to search.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

get_haplotype_coalescence2 <- function(project,
                                      haplo_IDs,
                                      sampling_time,
                                      generations = 3,
                                      silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_vector_pos_int(haplo_IDs, zero_allowed = FALSE)
  assert_greq(length(unique(haplo_IDs)), 2, message = "must input at least 2 different haplo IDs")
  assert_single_pos_int(sampling_time, zero_allowed = TRUE)
  assert_single_pos_int(generations, zero_allowed = FALSE)
  assert_single_logical(silent)
  
  # check that project has contig lengths defined
  assert_non_null(project$genetic_parameters$contig_lengths)
  
  # check that project has a block tree computed
  assert_non_null(project$block_tree, message = "project has no block_tree object")
  
  # check that all haplo_IDs are within the block tree
  if (!all(haplo_IDs %in% unique(project$block_tree$child_haplo_ID))) {
    stop("haplo_IDs could not be found within block_tree child_haplo_ID column")
  }
  
  # check that contig values are within range specified by genetic parameters
  if (!all(project$block_tree$contig <= length(project$genetic_parameters$contig_lengths))) {
    stop("block_tree contains contig values that exceed the number of contigs found in genetic_parameters")
  }
  
  # define arguments
  args <- list(time = project$block_tree$time,
               child_haplo_ID = project$block_tree$child_haplo_ID,
               contig = project$block_tree$contig - 1,
               contig_start = project$block_tree$start,
               contig_end = project$block_tree$end,
               parent_haplo_ID = project$block_tree$parent_haplo_ID,
               target_haplo_IDs = haplo_IDs,
               sampling_time = sampling_time,
               contig_lengths = project$genetic_parameters$contig_lengths,
               generations = generations,
               silent = silent)
  
  # run efficient C++ function
  output_raw <- get_haplotype_coalescence2_cpp(args)
  
  # process output
  output_processed <- NULL
  if (length(output_raw) > 0) {
    output_processed <- mapply(function(x) {
      as.data.frame(x[1:5])
    }, output_raw$coalescence_events, SIMPLIFY = FALSE) %>%
      dplyr::bind_rows()
    output_processed$descendants <- mapply(function(x) x$descendants, output_raw$coalescence_events)
  }
  
  return(output_processed)
}

#------------------------------------------------
#' @title Simulate mutations and write vcf to file
#'
#' @description TODO
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param output_location where to write vcfs to file.
#' @param sampling_time time of sample.
#' @param generations how many generations to look back when computing ancestral
#'   recombination graph.
#' @param mu mutation rate per base per unit time.
#' @param silent whether to suppress written messages to the console.
#'
#' @importFrom  stats rpois runif
#' @export

sim_vcf <- function(project,
                    output_location,
                    sampling_time,
                    generations = 3,
                    mu = 1e-7,
                    silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_single_string(output_location)
  assert_single_pos_int(sampling_time, zero_allowed = TRUE)
  assert_single_pos_int(generations, zero_allowed = FALSE)
  assert_single_pos(mu, zero_allowed = FALSE)
  assert_single_logical(silent)
  
  # check sample exists and contains at least two haplotypes
  assert_non_null(project$sample_details, message = "no sample details found in project")
  target_haplos <- unlist(project$sample_details$haplo_IDs)
  assert_greq(length(target_haplos), 2, message = "sample must contain at least two haplotype IDs within malaria-positive individuals")
  
  # get coalescence tree between blocks in format2
  message("Getting ancestral recombination graph")
  haplo_coalescence <- get_haplotype_coalescence2(project, haplo_IDs = target_haplos,
                                                  sampling_time = sampling_time,
                                                  generations = generations,
                                                  silent = silent)
  
  # simulate mutations
  n_haplos <- length(target_haplos)
  n_contigs <- length(project$genetic_parameters$contig_lengths)
  mut_map <- replicate(n_haplos, replicate(n_contigs, NULL), simplify = FALSE)
  
  message("Simulating mutations")
  for (i in 1:nrow(haplo_coalescence)) {
    n_muts <- rpois(1, mu*(haplo_coalescence$time_new[i] - haplo_coalescence$time_old[i])*(haplo_coalescence$right[i] - haplo_coalescence$left[i]))
    if (n_muts > 0) {
      mut_pos <- floor(runif(n_muts, haplo_coalescence$left[i], haplo_coalescence$right[i] + 1))
      w <- match(haplo_coalescence$descendants[i][[1]], target_haplos)
      for (j in seq_along(w)) {
        mut_map[[w[j]]][[haplo_coalescence$contig[i]]] <- c(mut_map[[w[j]]][[haplo_coalescence$contig[i]]], mut_pos)
      }
    }
  }
  
  # replace NULL with -1 before C++
  for (i in seq_along(mut_map)) {
    for (j in seq_along(mut_map[[i]])) {
      if (is.null(mut_map[[i]][[j]])) {
        mut_map[[i]][[j]] <- -1
      } else {
        mut_map[[i]][[j]] <- sort(unique(mut_map[[i]][[j]]))
      }
    }
  }
  
  # write each sequence to a different vcf
  for (i in seq_along(mut_map)) {
    
    args <- list(mut_map = mut_map[[i]],
                 output_location = sprintf("%s%s.vcf", output_location, i),
                 silent = FALSE)
    
    write_vcf_cpp(args)
  }
  
  
}
