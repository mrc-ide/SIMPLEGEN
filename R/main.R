
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
#' @param max_inoculations maximum number of inoculations that an individual
#'   can hold simultaneously.
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
                                        max_inoculations = 5,
                                        H = 1000,
                                        seed_infections = 100,
                                        M = 1000,
                                        mig_mat = diag(1),
                                        life_table = life_table_Mali()) {
  
  # NB. This function is written so that only parameters specified by the user
  # are updated. Any parameters that already have values within the project are
  # left alone
  
  # basic checks on inputs (more thorough checks on parameter values will be
  # carried out later)
  assert_class(project, "simplegen_project")
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values where not specified by user
  if (is.null(project$epi_model_parameters)) {
    project$epi_model_parameters <- within(as.list(environment()), rm(project))
  }
  
  # find which parameters are user-defined
  userlist <- as.list(match.call())
  userlist <- userlist[!(names(userlist) %in% c("", "project"))]
  
  # replace project parameters with user-defined
  project$epi_model_parameters[names(userlist)] <- mapply(eval, userlist, SIMPLIFY = FALSE)
  
  # standardise and process parameters
  project <- process_epi_model_params(project)
  
  # perform checks on final parameters
  check_epi_model_params(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
# convert epi model parameters to standardised types
#' @noRd
process_epi_model_params <- function(project) {
  
  # function that forces objects that can be defined as list or vector to list
  force_list <- function(x) {
    if (!is.list(x)) {
      x <- list(x)
    }
    return(x)
  }
  
  # force parameters to list
  name_vec <- c("duration_acute", "duration_chronic",
                "detectability_microscopy_acute", "detectability_microscopy_chronic",
                "detectability_PCR_acute", "detectability_PCR_chronic",
                "time_treatment_acute", "time_treatment_chronic",
                "infectivity_acute", "infectivity_chronic")
  project$epi_model_parameters[name_vec] <- mapply(force_list, project$epi_model_parameters[name_vec], SIMPLIFY = FALSE)
  
  # return
  invisible(project)
}

#------------------------------------------------
# perform checks on epi model parameters
#' @noRd
check_epi_model_params <- function(project) {
  
  # check project class
  assert_class(project, "simplegen_project")
  
  # check that epi model parameters exist
  assert_non_null(project$epi_model_parameters, message = "no epi model parameters defined. See ?define_epi_model_parameters")
  
  # extract model parameters
  x <- project$epi_model_parameters
  
  # perform checks
  assert_single_bounded(x$a, name = "a")
  assert_single_bounded(x$p, name = "p")
  assert_single_pos(x$mu, name = "mu")
  
  assert_single_pos_int(x$u, zero_allowed = FALSE, name = "u")
  assert_single_pos_int(x$v, zero_allowed = FALSE, name = "v")
  assert_single_pos_int(x$g, zero_allowed = FALSE, name = "g")
  
  assert_bounded(x$prob_infection, name = "prob_infection")
  assert_bounded(x$prob_acute, name = "prob_acute")
  assert_bounded(x$prob_AC, name = "prob_AC")
  
  mapply(assert_pos, x$duration_acute, name = "duration_acute")
  mapply(assert_pos, x$duration_chronic, name = "duration_chronic")
  mapply(assert_pos, x$detectability_microscopy_acute, name = "detectability_microscopy_acute")
  mapply(assert_pos, x$detectability_microscopy_chronic, name = "detectability_microscopy_chronic")
  mapply(assert_pos, x$detectability_PCR_acute, name = "detectability_PCR_acute")
  mapply(assert_pos, x$detectability_PCR_chronic, name = "detectability_PCR_chronic")
  mapply(assert_bounded, x$time_treatment_acute, name = "time_treatment_acute")
  mapply(assert_bounded, x$time_treatment_chronic, name = "time_treatment_chronic")
  
  assert_single_bounded(x$treatment_seeking_mean, name = "treatment_seeking_mean")
  assert_single_pos(x$treatment_seeking_sd, zero_allowed = TRUE, name = "treatment_seeking_sd")
  if (x$treatment_seeking_mean > 0 & x$treatment_seeking_mean < 1) {
    sd_max <- sqrt(x$treatment_seeking_mean*(1 - x$treatment_seeking_mean))
    error_message <- sprintf(paste0("treatment_seeking_sd must be less than",
                                    " sqrt(treatment_seeking_mean*(1 - treatment_seeking_mean)),",
                                    " in this case %s, otherwise the Beta distribution is undefined"), sd_max)
    assert_le(x$treatment_seeking_sd, sd_max, message = error_message)
  }
  
  assert_pos(x$duration_prophylactic, name = "duration_prophylactic")
  mapply(assert_bounded, x$infectivity_acute, name = "infectivity_acute")
  mapply(assert_bounded, x$infectivity_chronic, name = "infectivity_chronic")
  
  assert_single_pos_int(x$max_inoculations, zero_allowed = FALSE, name = "max_inoculations")
  assert_pos_int(x$H, name = "H")
  assert_pos_int(x$seed_infections, name = "seed_infections")
  assert_pos_int(x$M, name = "M")
  assert_same_length_multiple(x$H, x$seed_infections, x$M)
  assert_leq(x$seed_infections, x$H, name_x = "seed_infections", name_y = "H")
  
  assert_square_matrix(x$mig_mat, name = "mig_mat")
  assert_nrow(x$mig_mat, length(x$H), name = "mig_mat")
  assert_bounded(x$mig_mat, name = "mig_mat")
  
  assert_bounded(x$life_table, name = "life_table")
  assert_eq(x$life_table[length(x$life_table)], 1, message = "the final value in the life table must be 1, representing a 100%% chance of dying, to ensure a closed population")
  
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
#' @title Define the outputs that are produced from the transmission model
#'
#' @description Loads a dataframe into the SIMPLEGEN project that specifies how
#'   one or more samples are taken from the population. This is an important
#'   step in the SIMPLEGEN pipeline, as ultimately genotypes are only
#'   constructed for these sampled individuals. The user can specify details of
#'   the sampling scheme, including the type of sampling and the
#'   locations/timepoints. Probabilities of detection are taken into account
#'   when producing the final sample.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param df_sample a dataframe containing all of the following columns:
#'   \itemize{
#'     \item time: the time (in days) at which samples are taken.
#'     \item deme: the deme from which samples are taken.
#'     \item case_detection: the method by which cases are identified. Either
#'     "active" or "passive".
#'     \item diagnosis: the method by which infected individuals are diagnosed.
#'     Either "microscopy" or "PCR".
#'     \item n: the number of individuals screened. Note that the actual number
#'     of infected individuals (and hence the number of genotypes) may be lower
#'     than this number.
#'   }
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
  
  # if there are no defined parameters then create all parameters from
  # scratch using default values where not specified by user
  if (is.null(project$epi_sampling_parameters)) {
    project$epi_sampling_parameters <- within(as.list(environment()), rm(project))
  }
  
  # find which parameters are user-defined
  userlist <- as.list(match.call())
  userlist <- userlist[!(names(userlist) %in% c("", "project"))]
  
  # replace project parameters with user-defined
  project$epi_sampling_parameters[names(userlist)] <- mapply(eval, userlist, SIMPLIFY = FALSE)
  
  # perform checks on final parameters
  check_epi_sampling_params(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
# perform checks on epi sampling parameters
#' @noRd
check_epi_sampling_params <- function(project) {
  
  # check project class
  assert_class(project, "simplegen_project")
  
  # model parameters must be defined before sampling parameters
  assert_non_null(project$epi_model_parameters, message = "model parameters must be defined before sampling parameters. See ?define_epi_model_parameters")
  
  # check that epi sampling parameters exist
  assert_non_null(project$epi_sampling_parameters, message = "no epi sampling parameters defined. See ?define_epi_sampling_parameters")
  
  # check individual elements
  check_epi_sampling_params_daily(project$epi_sampling_parameters$daily)
  check_epi_sampling_params_sweeps(project$epi_sampling_parameters$sweeps)
}

#------------------------------------------------
# perform checks on daily sampling parameters
#' @noRd
check_epi_sampling_params_daily <- function(x) {
  
  # return if null
  if (is.null(x)) {
    return()
  }
  
  # check dataframe column names
  assert_dataframe(x, message = "daily sampling parameters must be a dataframe")
  col_titles <- c("deme", "measure", "state", "diagnostic", "age_min", "age_max", "inoculations")
  assert_in(col_titles, names(x), message = sprintf("daily sampling parameters dataframe must contain the following columns: {%s}",
                                                    paste0(col_titles, collapse = ", ") ))
  
  # check deme and measure formats
  deme_mssg <- "deme must be a positive integer or -1"
  assert_vector_int(x$deme, message = deme_mssg)
  assert_greq(x$deme, -1, message = deme_mssg)
  assert_greq(x$deme[x$deme != -1], 1, message = deme_mssg)
  
  x$measure <- as.character(x$measure)
  measure_levels <- c("count", "prevalence", "incidence", "EIR")
  assert_in(x$measure, measure_levels, message = sprintf("measure must be one of: {%s}", paste0(measure_levels, collapse = ", ") ))
  
  # split into sub-dataframes based on measure, and check state and diagnostic columns
  if (any(x$measure == "EIR")) {
    
    df_EIR <- subset(x, measure == "EIR")
    
    # check state and detection columns
    assert_NA(df_EIR$state, message = "state must be NA when measure is EIR")
    assert_NA(df_EIR$diagnostic, message = "diagnostic must be NA when measure is EIR")
    
  }
  if (any(x$measure != "EIR")) {
    
    df_main = subset(x, measure != "EIR")
    
    # check state and diagnostic columns
    state_levels <- c("S", "E", "A", "C", "P", "H", "Sv", "Ev", "Iv", "M")
    assert_in(df_main$state, state_levels, message = sprintf("state must be one of {%s}", paste0(state_levels, collapse = ", ")))
    diagnostic_levels <- c("true", "microscopy", "PCR")
    assert_in(df_main$diagnostic, diagnostic_levels, message = sprintf("diagnostic must be one of {%s}", paste0(diagnostic_levels, collapse = ", ")))
    
  }
  
  # check age_min, age_max and inoculations columns
  assert_pos_int(x$age_min, zero_allowed = TRUE, message = "age_min must be a positive integer or zero")
  assert_pos_int(x$age_max, zero_allowed = TRUE, message = "age_max must be a positive integer or zero")
  assert_greq(x$age_max, x$age_min, message = "age_max must be greater than or equal to age_min")
  
  inoc_mssg <- "inoculations must be a positive integer or -1 to indicate any number of inoculations"
  assert_vector_int(x$inoculations, message = inoc_mssg)
  assert_greq(x$inoculations, -1, message = inoc_mssg)
  
}

#------------------------------------------------
# perform checks on population sweep sampling parameters
#' @noRd
check_epi_sampling_params_sweeps <- function(x) {
  
  # return if null
  if (is.null(x)) {
    return()
  }
  
  # check dataframe column names
  assert_dataframe(x, message = "sweep sampling parameters must be a dataframe")
  col_titles <- c("time", "deme", "measure", "state", "diagnostic", "age_min", "age_max", "inoculations")
  assert_in(col_titles, names(x), message = sprintf("sweep sampling parameters dataframe must contain the following columns: {%s}",
                                                    paste0(col_titles, collapse = ", ") ))
  
  # check time format
  assert_vector_pos_int(x$time, zero_allowed = FALSE, message = "time must be a positive integer")
  
  # remaining columns should have identical format to daily dataframe
  check_epi_sampling_params_daily(subset(x, select = -time))
  
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
#' @param output_format several options exist for the output format:
#'   \itemize{
#'     \item 1 (default) = return final values only
#'     \item 2 = also return numerator and denominator of prevalence and
#'     incidence calculations.
#'   }
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
                    output_format = 1,
                    save_transmission_record = FALSE,
                    transmission_record_location = "",
                    overwrite_transmission_record = FALSE,
                    pb_markdown = FALSE,
                    silent = FALSE) {
  
  
  # ---------- check inputs ----------
  
  assert_class(project, "simplegen_project")
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_in(output_format, c(1,2))
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
                 output_format = output_format,
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
  args$any_sweep_outputs <- !is.null(args$sweep)
  
  # get sampling strategy indices into 0-indexed (C++) format
  sampling_to_cpp_format <- function(x) {
    w <- which(x$deme != -1)
    x$deme[w] <- x$deme[w] - 1
    x$measure <- match(x$measure, c("count", "prevalence", "incidence", "EIR")) - 1
    x$state <- match(x$state, c("S", "E", "A", "C", "P", "H", "Sv", "Ev", "Iv", "M")) - 1
    x$diagnostic <- match(x$diagnostic, c("true", "microscopy", "PCR")) - 1
    return(x)
  }
  if (args$any_daily_outputs) {
    args$daily <- sampling_to_cpp_format(args$daily)
  }
  if (args$any_sweep_outputs) {
    args$sweeps <- sampling_to_cpp_format(args$sweeps)
  }
  
  # unique times at which sweeps happen
  args$sweep_time_ordered <- sort(unique(args$sweeps$time))
  
  # functions
  args_functions <- list(update_progress = update_progress)
  
  # make progress bars
  pb_sim <- txtProgressBar(min = 0, max = max_time, initial = NA, style = 3)
  args_progress <- list(pb_sim = pb_sim)
  
  
  # ---------- run simulation ----------
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args, args_functions, args_progress)
  
  
  # ---------- process output ----------
  
  # wrangle daily output
  daily_output <- NULL
  daily_sampling <- project$epi_sampling_parameters$daily
  if (!is.null(daily_sampling)) {
    
    # calculate final values from numerator and denominator
    daily_numer <- unlist(output_raw$daily_numer)
    daily_denom <- unlist(output_raw$daily_denom)
    daily_values <- daily_numer / daily_denom
    
    # make output dataframe with raw values
    daily_output <- cbind(time = rep(seq_len(max_time), each = nrow(daily_sampling)),
                          daily_sampling,
                          value = daily_values,
                          numer = daily_numer,
                          denom = daily_denom,
                          row.names = NULL)
    
    # insert NAs as needed
    w <- which(!(daily_output$measure %in% c("prevalence", "incidence")))
    daily_output$value[w] <- daily_output$numer[w]
    daily_output$numer[w] <- daily_output$denom[w] <- NA
    
    # format dataframe
    if (output_format == 1) {
      daily_output <- subset(daily_output, select = -c(numer, denom))
    }
  }
  
  # wrangle sweep output
  sweeps_output <- NULL
  sweeps_sampling <- project$epi_sampling_parameters$sweeps
  if (!is.null(sweeps_sampling)) {
    
    # calculate final values from numerator and denominator
    sweep_numer <- unlist(output_raw$sweep_numer)
    sweep_denom <- unlist(output_raw$sweep_denom)
    sweep_values <- sweep_numer / sweep_denom
    
    # make output dataframe with raw values
    sweeps_output <- cbind(sweeps_sampling,
                           value = sweep_values,
                           numer = sweep_numer,
                           denom = sweep_denom,
                           row.names = NULL)
    
    # insert NAs as needed
    w <- which(!(sweeps_output$measure %in% c("prevalence", "incidence")))
    sweeps_output$value[w] <- sweeps_output$numer[w]
    sweeps_output$numer[w] <- sweeps_output$denom[w] <- NA
    
    # format dataframe
    if (output_format == 1) {
      sweeps_output <- subset(sweeps_output, select = -c(numer, denom))
    }
  }
  
  # wrangle surveys output
  surveys_output <- NULL
  
  # append to project
  project$epi_output <- list(daily = daily_output,
                             sweeps = sweeps_output,
                             surveys = surveys_output)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Prune the transmission record
#'
#' @description Reads in a saved transmission record from file. Combines this
#'   with sampling information in the project to produce a pruned version of the
#'   transmission record containing only the events relevant to the final
#'   sample. This is also saved to file.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param transmission_record_location the file path to a transmission record
#'   already written to file.
#' @param pruned_record_location the file path that the pruned transmission
#'   record will be written to.
#' @param overwrite_pruned_record if \code{TRUE} the pruned transmission record
#'   will overwrite any existing file by the same name. \code{FALSE} by default.
#' @param silent whether to suppress written messages to the console.
#' 
#' @export

prune_transmission_record <- function(project,
                                      transmission_record_location = "",
                                      pruned_record_location = "",
                                      overwrite_pruned_record = FALSE,
                                      silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_string(transmission_record_location)
  assert_string(pruned_record_location)
  assert_single_logical(overwrite_pruned_record)
  assert_single_logical(silent)
  
  # check transmission record exists
  if (!file.exists(transmission_record_location)) {
    stop(sprintf("could not find file at %s", transmission_record_location))
  }
  
  # optionally return warning if will overwrite pruned transmission record file
  if (!overwrite_pruned_record) {
    if (file.exists(pruned_record_location)) {
      stop(sprintf("file already exists at %s. Change target location, or use argument `overwrite_pruned_record = TRUE` to manually override this warning", pruned_record_location))
    }
  }
  
  # subset sample details to a vector of inoc_IDs
  inoc_IDs <- unlist(project$sample_output$inoc_IDs)
  if (length(inoc_IDs) == 0) {
    stop("no malaria positive hosts in sample")
  }
  
  # define argument list
  args <- list(transmission_record_location = transmission_record_location,
               pruned_record_location = pruned_record_location,
               inoc_IDs = inoc_IDs,
               silent = silent)
  
  # run efficient C++ code
  prune_transmission_record_cpp(args)
  
  # return project unchanged
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
#' @param r the rate of recombination. The expected number of base pairs in a
#'   single recombinant block is 1/r.
#' @param alpha parameter dictating the skew of lineage densities. Small
#'   values of \code{alpha} create a large skew, and hence make it likely that
#'   an oocyst will be produced from the same parents. Large values of
#'   \code{alpha} tend toward more even densities.
#' @param oocyst_distribution vector specifying the probability distribution of
#'   each number of oocysts within the mosquito midgut.
#' @param hepatocyte_distribution vector specifying the probability distribution
#'   of the number of infected hepatocytes in a human host. More broadly, this
#'   defines the number of independent draws from the oocyst products that make
#'   it into the host bloodstream upon a bite from an infectious mosquito.
#' @param contig_lengths vector of lengths (in bp) of each contig.
#' 
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
  
  # if there are no defined genetic parameters then create all parameters from
  # scratch using default values where not specified by user
  if (is.null(project$genetic_parameters)) {
    project$genetic_parameters <- as.list(environment())
    invisible(project)
  }
  
  # find which parameters are user-defined
  userlist <- as.list(match.call())
  userlist <- userlist[!(names(userlist) %in% c("", "project"))]
  
  # replace project parameters with user-defined
  project$genetic_parameters[names(userlist)] <- mapply(eval, userlist, SIMPLIFY = FALSE)
  
  # perform checks on parameters
  check_genetic_params(project$genetic_parameters)
  
  # return
  invisible(project)
}

#------------------------------------------------
# perform checks on genetic parameters
#' @noRd
check_genetic_params <- function(x) {
  
  # perform checks
  assert_single_pos(x$r, zero_allowed = TRUE)
  assert_single_pos(x$alpha, zero_allowed = FALSE)
  assert_vector_pos(x$oocyst_distribution)
  assert_vector_pos(x$hepatocyte_distribution)
  assert_vector_pos_int(x$contig_lengths, zero_allowed = FALSE)
  
}

#------------------------------------------------
#' @title Draw tree of genome-wide relatedness from pruned transmission record
#'
#' @description Reads in the pruned transmission record from file and uses this
#'   information to simulate relatedness between parasite lineages from a
#'   genetic model. The map of relatedness output by this step can be used to
#'   generate genotypic information of various types.
#'   
#'
#' @details Reads in the pruned transmission record and creates a new node for
#'   each inoculation ID. Nodes at time zero are initialised with a single
#'   unique lineage created de novo. In subsequent generations the children of
#'   any given node are known from the pruned transmission record. The lineages
#'   from parental nodes are sampled at random, and brought together in pairs to
#'   produce oocysts. The recombinant products of these oocysts are then sampled
#'   down to produce a new generation of lineage IDs for this node, along with
#'   the relative densities of each lineage.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param pruned_record_location the file path from which the pruned
#'   transmission record will be read.
#' @param silent whether to suppress written messages to the console.
#'
#' @importFrom stats dpois
#' @export

sim_relatedness <- function(project,
                            pruned_record_location = "",
                            silent = FALSE) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_string(pruned_record_location)
  assert_neq(pruned_record_location, "", message = "pruned_record_location cannot be empty")
  assert_single_logical(silent)
  
  # check pruned record exists
  if (!file.exists(pruned_record_location)) {
    stop(sprintf("could not find file at %s", pruned_record_location))
  }
  
  # check for defined genetic params
  assert_non_null(project$genetic_parameters, message = "no genetic parameters defined. See ?define_genetic_params")
  
  # define arguments
  args <- append(project$genetic_parameters,
                 list(pruned_record_location = pruned_record_location,
                      silent = silent))
  
  # run efficient C++ code
  output_raw <- sim_relatedness_cpp(args)
  
  # get the inoculation IDs of every node in the pruned transmission record
  inoc_IDs <- mapply(function(x) x$details$inoc_ID, output_raw)
  
  # get the lineage IDs of the sampled hosts
  sample_lineage_IDs <- list()
  for (i in seq_len(nrow(project$sample_output))) {
    
    # get lineage IDs from inoculation IDs
    node_inoc_IDs <- project$sample_output$inoc_IDs[[i]]
    if (length(node_inoc_IDs) == 0) {
      sample_lineage_IDs[[i]] <- integer()
    } else {
      # which output elements are we interested in
      w <- match(node_inoc_IDs, inoc_IDs)
      
      # concatenate all lineage IDs over these elements
      sample_lineage_IDs[[i]] <- unlist(mapply(function(x) x$details$lineage_IDs, output_raw[w], SIMPLIFY = FALSE))
    }
  }
  
  # add sampled lineage IDs back into project
  project$sample_output$lineage_IDs <- sample_lineage_IDs
  
  # get list over lineages, ignoring hosts
  lineage_list <- unlist(mapply(function(k) {
    ret <- output_raw[[k]]$lineages
  }, seq_along(output_raw), SIMPLIFY = FALSE), recursive = FALSE)
  
  # add lineage list back into project
  project$relatedness <- lineage_list
  
  # return project
  invisible(project)
}

#------------------------------------------------
#' @title Get coalescent time of two lineages
#'
#' @description TODO.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param lineage_IDs a vector of two lineage IDs with which to calculate coalescent times.
#' @param max_reps maximum number of steps when searching for coalescent times.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

get_coalescent_times <- function(project,
                                 lineage_IDs,
                                 max_reps = 1e3,
                                 silent = FALSE) {
  
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_vector_int(lineage_IDs)
  assert_length(lineage_IDs, 2)
  assert_single_pos_int(max_reps, zero_allowed = FALSE)
  assert_single_logical(silent)
  
  # check lineage_IDs can be found in project output
  assert_in("lineage_IDs", names(project$sample_output),
            message = paste0("column 'lineage_IDs' not found in project$sample_output.",
                             "Check you have run all necessary steps in the pipeline up to this point."))
  assert_in(lineage_IDs, unlist(project$sample_output$lineage_IDs),
            message = "lineage_IDs not found in project$sample_output$lineage_IDs")
  
  # check project contains relatedness info
  assert_non_null(project$relatedness, message = "project contains no relatedness information")
  
  #------------------------------------------------
  
  # function replaces element x[[i]] with ancestral intervals
  # defined inside function as hard-coded for project$relatedness
  get_ancestor <- function(x, i, contig) {
    parent <- x[[i]][length(x[[i]])]
    if (parent == -1) {
      return(x)
    }
    x_ancestor <- project$relatedness[[parent]][[contig]]
    x_replace <- list()
    for (j in seq_along(x_ancestor)) {
      if (interval_intersect(x_ancestor[[j]], x[[i]])) {
        left <- max(x[[i]][1], x_ancestor[[j]][1])
        right <- min(x[[i]][2], x_ancestor[[j]][2])
        x_replace <- append(x_replace, list(c(left, right, x[[i]][-(1:2)], x_ancestor[[j]][3])))
      }
    }
    ret <- append(x[-i], x_replace, i-1)
    return(ret)
  }
  
  # loop through contigs
  contigs <- length(project$relatedness[[lineage_IDs[1]]])
  ret_list <- list()
  for (contig in seq_len(contigs)) {
    message(sprintf("contig %s of %s", contig, contigs))
    
    # get starting relatedness info for both lineages
    x1 <- project$relatedness[[lineage_IDs[1]]][[contig]]
    x2 <- project$relatedness[[lineage_IDs[2]]][[contig]]
    
    # loop back through ancestry, splitting intervals and checking for coalescence
    for (rep in 1:max_reps) {
      
      # get complete list of breakpoints over both x1 and x2
      all_breaks <- c(mapply(function(x) x[1], x1),
                      mapply(function(x) x[1], x2))
      all_breaks <- sort(unique(all_breaks))
      all_breaks <- c(all_breaks, x1[[length(x1)]][2] + 1)
      
      # make intervals the same between x1 and x2
      for (i in 1:(length(all_breaks) - 1)) {
        if (x1[[i]][2] != all_breaks[i+1] - 1) {
          x1 <- append(x1, list(c(x1[[i]][1], all_breaks[i+1] - 1, x1[[i]][-(1:2)])), i - 1)
          x1[[i+1]][1] <- all_breaks[i+1]
        }
        if (x2[[i]][2] != all_breaks[i+1] - 1) {
          x2 <- append(x2, list(c(x2[[i]][1], all_breaks[i+1] - 1, x2[[i]][-(1:2)])), i - 1)
          x2[[i+1]][1] <- all_breaks[i+1]
        }
      }
      
      # find next element of x1/x2 that has not already coalesced
      all_coal <- TRUE
      for (i in seq_along(x1)) {
        if (length(intersect(x1[[i]][-(1:2)], x2[[i]][-(1:2)])) == 0) {
          all_coal <- FALSE
          break
        }
      }
      if (all_coal) {
        break
      }
      
      # split uncoalesced elements of x1 and x2 based on ancestry 
      x1 <- get_ancestor(x1, i, contig)
      x2 <- get_ancestor(x2, i, contig)
      
    }  # end rep loop
    
    # check that didn't time out
    if (rep == max_reps) {
      stop("max_reps reached")
    }
    
    # make dataframe of coalescent intervals and timings
    x_coal <- x1
    for (i in seq_along(x_coal)) {
      coal_lineage <- intersect(x1[[i]][-(1:2)], x2[[i]][-(1:2)])
      gen1 <- which(x1[[i]][-(1:2)] == coal_lineage)
      gen2 <- which(x2[[i]][-(1:2)] == coal_lineage)
      x_coal[[i]] <- c(x_coal[[i]][1:2], coal_lineage, gen1, gen2)
    }
    x_coal <- as.data.frame(do.call(rbind, x_coal))
    names(x_coal) <- c("start", "end", "ancestor", "generations1", "generations2")
    
    # simplify by merging adjacent intervals with same ancestor
    i <- 2
    while (i <= nrow(x_coal)) {
      if (x_coal$ancestor[i] == x_coal$ancestor[i-1]) {
        x_coal$end[i-1] <- x_coal$end[i]
        x_coal <- x_coal[-i,]
      } else {
        i <- i + 1
      }
    }
    
    # add to list over contigs
    ret_list[[contig]] <- x_coal
  }
  
  # convert to dataframe
  ret_df <- do.call(rbind, mapply(function(i) {
    cbind(contig = i, ret_list[[i]])
  }, seq_along(ret_list), SIMPLIFY = FALSE))
  
  # add segment length column
  ret_df <- as.data.frame(append(as.list(ret_df), list(length = ret_df$end - ret_df$start + 1), 3))
  
  return(ret_df)
}
