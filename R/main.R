
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

define_epi_params <- function(project,
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
  assert_custom_class(project, "simplegen_project")
  
  # if the project is completely new then create all parameters de novo, using
  # default values where not specified by user
  if (is.null(project$epi_parameters)) {
    project$epi_parameters <- as.list(environment())
    invisible(project)
  }
  
  # find which parameters are user-defined
  userlist <- as.list(match.call())
  userlist <- userlist[!(names(userlist) %in% c("", "project"))]
  
  # replace project parameters with user-defined
  project$epi_parameters[names(userlist)] <- mapply(eval, userlist, SIMPLIFY = FALSE)
  
  # standardise parameters (e.g. normalise distributions)
  params_processed <- process_epi_params(project$epi_parameters)
  
  # perform checks on parameters
  check_epi_params(params_processed)
  
  # return
  invisible(project)
}

#------------------------------------------------
# convert epi parameters to standardised types
#' @noRd
process_epi_params <- function(x) {
  
  # for objects that can be defined as list or vector, force to list
  if (!is.list(x$duration_acute)) {
    x$duration_acute <- list(x$duration_acute)
  }
  if (!is.list(x$duration_chronic)) {
    x$duration_chronic <- list(x$duration_chronic)
  }
  if (!is.list(x$detectability_microscopy_acute)) {
    x$detectability_microscopy_acute <- list(x$detectability_microscopy_acute)
  }
  if (!is.list(x$detectability_microscopy_chronic)) {
    x$detectability_microscopy_chronic <- list(x$detectability_microscopy_chronic)
  }
  if (!is.list(x$detectability_PCR_acute)) {
    x$detectability_PCR_acute <- list(x$detectability_PCR_acute)
  }
  if (!is.list(x$detectability_PCR_chronic)) {
    x$detectability_PCR_chronic <- list(x$detectability_PCR_chronic)
  }
  if (!is.list(x$time_treatment_acute)) {
    x$time_treatment_acute <- list(x$time_treatment_acute)
  }
  if (!is.list(x$time_treatment_chronic)) {
    x$time_treatment_chronic <- list(x$time_treatment_chronic)
  }
  if (!is.list(x$infectivity_acute)) {
    x$infectivity_acute <- list(x$infectivity_acute)
  }
  if (!is.list(x$infectivity_chronic)) {
    x$infectivity_chronic <- list(x$infectivity_chronic)
  }
  
  return(x)
}

#------------------------------------------------
# perform checks on parameters
#' @noRd
check_epi_params <- function(x) {
  
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
#' @title Define how samples are taken from transmission model
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
#' @param x a dataframe containing all of the following columns:
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

define_sampling_strategy <- function(project, x) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  assert_dataframe(x)
  assert_eq(names(x), c("time", "deme", "case_detection", "diagnosis", "n"))
  assert_pos_int(x$time, zero_allowed = FALSE)
  assert_pos_int(x$deme, zero_allowed = FALSE)
  assert_in(x$case_detection, c("active", "passive"))
  assert_in(x$diagnosis, c("microscopy", "PCR"))
  assert_pos_int(x$n, zero_allowed = FALSE)
  
  # specify formats
  x$case_detection <- as.character(x$case_detection)
  x$diagnosis <- as.character(x$diagnosis)
  
  # load into project
  project$sampling_strategy <- x
  
  invisible(project)
}

#------------------------------------------------
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from the inbuilt SIMPLEGEN transmission model.
#'   Parameters are taken from the \code{epi_parameters} slot of the project,
#'   and basic outputs are written to the \code{epi_output} slot. If a sampling
#'   strategy has been defined then samples will also be obtained and saved in
#'   the \code{sample_details} slot (see \code{?define_sampling_strategy()}).
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
#' @param output_daily_counts whether to output daily counts of key quantities,
#'   such as the number of infected hosts and the EIR.
#' @param output_age_distributions whether to output complete age distributions
#'   at certain times.
#' @param output_age_times a vector of times at which complete age distributions
#'   are output.
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
                    output_daily_counts = TRUE,
                    output_age_distributions = TRUE,
                    output_age_times = max_time,
                    pb_markdown = FALSE,
                    silent = FALSE) {
  
  
  # ---------- check inputs ----------
  
  assert_custom_class(project, "simplegen_project")
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_single_logical(save_transmission_record)
  assert_string(transmission_record_location)
  assert_single_logical(overwrite_transmission_record)
  assert_single_logical(output_daily_counts)
  assert_single_logical(output_age_distributions)
  assert_vector_pos_int(output_age_times, zero_allowed = FALSE)
  assert_leq(output_age_times, max_time)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # optionally return warning if will overwrite transmission record file
  if (save_transmission_record & !overwrite_transmission_record) {
    if (file.exists(transmission_record_location)) {
      stop(sprintf("file already exists at %s. Change target location, or use argument `overwrite_transmission_record = TRUE` to manually override this warning", transmission_record_location))
    }
  }
  
  # check for defined epi params
  if (is.null(project$epi_parameters)) {
    stop("no epi parameters defined. See ?define_epi_params")
  }
  
  
  # ---------- define argument lists ----------
  
  # get project params into standardised format and perform checks
  args <- process_epi_params(project$epi_parameters)
  check_epi_params(args)
  
  # get migration matrix into list
  args$mig_mat <- matrix_to_rcpp(args$mig_mat)
  
  # get complete demography from life table
  demog <- get_demography(project$epi_parameters$life_table)
  args <- c(args,
            list(age_death = demog$age_death,
                 age_stable = demog$age_stable))
  
  # add sampling strategy info
  if (is.null(project$sampling_strategy)) {
    args <- c(args, obtain_samples = FALSE)
  } else {
    ss <- project$sampling_strategy
    args <- c(args,
              list(obtain_samples = TRUE,
                   ss_time = ss$time,
                   ss_deme = ss$deme - 1,  # NB, subtract 1 to go from R to C++ indexing
                   ss_case_detection = ss$case_detection,
                   ss_diagnosis = ss$diagnosis,
                   ss_n = ss$n))
  }
  
  # append arguments to this function to the list of args that are passed forward
  args <- c(args,
            list(max_time = max_time,
                 save_transmission_record = save_transmission_record,
                 transmission_record_location = transmission_record_location,
                 output_daily_counts = output_daily_counts,
                 output_age_distributions = output_age_distributions,
                 output_age_times = output_age_times,
                 pb_markdown = pb_markdown,
                 silent = silent))
  
  # functions
  args_functions <- list(update_progress = update_progress)
  
  # make progress bars
  pb_sim <- txtProgressBar(min = 0, max = max_time, initial = NA, style = 3)
  args_progress <- list(pb_sim = pb_sim)
  
  
  # ---------- run simulation ----------
  
  # internal flag, not visible to user. If TRUE then this function writes
  # parameter lists to file and returns without running simulation. Parameters
  # can then be read directly from file into Xcode.
  #xcode_on <- FALSE
  #if (xcode_on) {
  #  write_xcode_params(args)
  #  return()
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args, args_functions, args_progress)
  
  
  # ---------- process output ----------
  
  # wrangle daily values into dataframe
  daily_values <- NULL
  if (output_daily_counts) {
    daily_values <- do.call(rbind, mapply(function(i) {
      ret <- rcpp_to_matrix(output_raw$daily_values[[i]])
      ret <- as.data.frame(cbind(seq_len(nrow(ret)), i, ret))
      names(ret) <- c("time", "deme", "H", "S", "E", "A", "C", "P",
                      "Sv", "Ev", "Iv",
                      "EIR", "inc_infection", "inc_acute", "inc_chronic",
                      "A_detectable_microscopy", "C_detectable_microscopy",
                      "A_detectable_PCR", "C_detectable_PCR","n_inoc")
      return(ret)
    }, seq_along(output_raw$daily_values), SIMPLIFY = FALSE))
  }
  
  # wrangle age distributions into dataframe
  age_distributions <- NULL
  if (output_age_distributions) {
    age_distributions <- do.call(rbind, mapply(function(j) {
      ret <- do.call(rbind, mapply(function(i) {
        ret <- do.call(rbind, output_raw$age_distributions[[j]][[i]])
        colnames(ret) <- c("S", "E", "A", "C", "P", "inc_infection", "inc_acute", "inc_chronic")
        data.frame(cbind(deme = i, age = seq_len(nrow(ret)) - 1, ret))
      }, seq_along(output_raw$age_distributions[[j]]), SIMPLIFY = FALSE))
      cbind(sample_time = j, ret)
    }, seq_along(output_raw$age_distributions), SIMPLIFY = FALSE))
  }
  
  # wrangle sample details into dataframe
  sample_details_list <- mapply(function(x) {
    ret <- data.frame(time = x[1], deme = x[2], host_ID = x[3], positive = x[4])
    ret$inoc_IDs <- list(x[-(1:4)])
    return(ret)
  }, output_raw$sample_details, SIMPLIFY = FALSE)
  sample_details <- do.call(rbind, sample_details_list)
  
  # append to project
  project$epi_output <- list(daily_values = daily_values,
                             age_distributions = age_distributions)
  if (!is.null(sample_details)) {
    project$sample_details <- sample_details
  }
  
  invisible(project)
}

