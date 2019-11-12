
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
#' @description The default life table used within the SIMPLEGEN epidemiological
#'   model. Values represent probabilities of death in one-year age bands, and
#'   are taken from (TODO - insert reference for this table).
#'
#' @export

life_table_Mali <- function() {
  life_table_raw <- simplegen_file("Mali_life_table.csv")
  ret <- life_table_raw$prop_death
  return(ret)
}

#------------------------------------------------
#' @title Define parameters of epidemiological simulation model
#'
#' @description Define the parameters that will be used when simulating data
#'   from the inbuilt SIMPLEGEN epidemiological model. Any parameters that are
#'   left \code{NULL} will be set to default values.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. \code{mu = -log(p)} unless
#'   otherwise specified.
#' @param u intrinsic incubation period. The number of days from infection to
#'   blood-stage in a human host.
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
#' @param treatment_seeking_alpha,treatment_seeking_beta treatment seeking is
#'   modelled in two stages; first, the time at which hosts \emph{consider}
#'   seeking treatment is drawn from the probability distribution in
#'   \code{time_treatment_acute} or \code{time_treatment_chronic}, and at this
#'   point whether they go ahead with treatment-seeking depends on their
#'   personal host-specific treatment-seeking probability. These latter
#'   probabilities are drawn independently for each host at birth from a Beta
#'   distribution with shape parameters `treatment_seeking_alpha` and
#'   `treatment_seeking_beta`. The mean of this distribution defines the average
#'   treatment probability in the population, and is equal to
#'   `treatment_seeking_alpha/(treatment_seeking_alpha +
#'   treatment_seeking_beta)`.
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
                              duration_chronic = dgeom(1:100, 1/20),
                              detectability_microscopy_acute = 1,
                              detectability_microscopy_chronic = 0.1,
                              detectability_PCR_acute = 1,
                              detectability_PCR_chronic = 1,
                              time_treatment_acute = dgeom(1:25, 1/5),
                              time_treatment_chronic = dgeom(1:100, 1/20),
                              treatment_seeking_alpha = 1,
                              treatment_seeking_beta = 1,
                              duration_prophylactic = dgeom(1:25, 1/5),
                              infectivity_acute = 0.07,
                              infectivity_chronic = 0.07,
                              max_inoculations = 5,
                              H = 1000,
                              seed_infections = 100,
                              M = 1000,
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
    project$epi_parameters <- list(a = a,
                                   p = p,
                                   mu = mu,
                                   u = u,
                                   v = v,
                                   g = g,
                                   prob_infection = prob_infection,
                                   prob_acute = prob_acute,
                                   prob_AC = prob_AC,
                                   duration_acute = duration_acute,
                                   duration_chronic = duration_chronic,
                                   detectability_microscopy_acute = detectability_microscopy_acute,
                                   detectability_microscopy_chronic = detectability_microscopy_chronic,
                                   detectability_PCR_acute = detectability_PCR_acute,
                                   detectability_PCR_chronic = detectability_PCR_chronic,
                                   time_treatment_acute = time_treatment_acute,
                                   time_treatment_chronic = time_treatment_chronic,
                                   treatment_seeking_alpha = treatment_seeking_alpha,
                                   treatment_seeking_beta = treatment_seeking_beta,
                                   duration_prophylactic = duration_prophylactic,
                                   infectivity_acute = infectivity_acute,
                                   infectivity_chronic = infectivity_chronic,
                                   max_inoculations = max_inoculations,
                                   H = H,
                                   seed_infections = seed_infections,
                                   M = M,
                                   life_table = life_table)
    
    return(project)
  }
  
  # find which parameters are user-defined
  userlist <- as.list(match.call()) %>%
    (function(x) x[!(names(x) %in% c("", "project"))])
  
  # replace project parameters with user-defined
  project$epi_parameters[names(userlist)] <- mapply(eval, userlist, SIMPLIFY = FALSE)
  
  # standardise parameters (e.g. normalise distributions) and perform checks
  params_processed <- process_params(project$epi_parameters)
  check_params(params_processed)
  
  # return
  invisible(project)
}

#------------------------------------------------
# convert parameters to standardised types
#' @noRd
process_params <- function(x) {
  
  # standardise parameters
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
check_params <- function(x) {
  
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
  assert_single_pos(x$treatment_seeking_alpha, name = "treatment_seeking_alpha")
  assert_single_pos(x$treatment_seeking_beta, name = "treatment_seeking_beta")
  assert_pos(x$duration_prophylactic, name = "duration_prophylactic")
  mapply(assert_bounded, x$infectivity_acute, name = "infectivity_acute")
  mapply(assert_bounded, x$infectivity_chronic, name = "infectivity_chronic")
  assert_single_pos_int(x$max_inoculations, zero_allowed = FALSE, name = "max_inoculations")
  assert_pos_int(x$H, name = "H")
  assert_pos_int(x$seed_infections, name = "seed_infections")
  assert_pos_int(x$M, name = "M")
  assert_same_length_multiple(x$H, x$seed_infections, x$M)
  assert_leq(x$seed_infections, x$H, name_x = "seed_infections", name_y = "H")
  assert_bounded(x$life_table, name = "life_table")
  assert_eq(x$life_table[length(x$life_table)], 1, message = "the final value in the life table must be 1, representing a 100%% chance of dying, to ensure a closed population")
  
}

#------------------------------------------------
# get stable demography distribution from life table
#' @noRd
get_demography <- function(life_table) {
  
  # compute distribution of age of death
  n <- length(life_table)
  age_death <- rep(0,n)
  remaining <- 1
  for (i in 1:n) {
    age_death[i] <- remaining*life_table[i]
    remaining <- remaining*(1 - life_table[i])
  }
  
  # convert life table to transition matrix
  m <- matrix(0,n,n)
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
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from the inbuilt epidemiological model. Parameters are
#'   taken from the \code{epi_parameters} slot of the project, and basic outputs are
#'   written to the \code{epi_output} slot.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param max_time run simulation for this many days.
#' @param save_transmission_record whether to write the transmission record to
#'   file.
#' @param transmission_record_location the file path that the transmission
#'   record will be written to.
#' @param output_daily_counts whether to output daily counts of key quantities,
#'   such as the number of infected hosts and the EIR.
#' @param output_age_distributions whether to output complete age distributions
#' @param output_age_times a vector of times at which complete age distributions
#'   are output.
#' @param silent whether to suppress written messages to the console.
#'
#' @export

sim_epi <- function(project,
                    max_time = 365,
                    save_transmission_record = FALSE,
                    transmission_record_location = "",
                    output_daily_counts = TRUE,
                    output_age_distributions = TRUE,
                    output_age_times = max_time,
                    silent = FALSE) {
  
  
  # ---------- check inputs ----------
  
  assert_custom_class(project, "simplegen_project")
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_single_logical(save_transmission_record)
  assert_string(transmission_record_location)
  assert_single_logical(output_daily_counts)
  assert_single_logical(output_age_distributions)
  assert_vector(output_age_times)
  assert_pos_int(output_age_times)
  assert_leq(output_age_times, max_time)
  assert_single_logical(silent)
  
  
  # ---------- define argument lists ----------
  
  # check for defined epi params
  if (is.null(project$epi_parameters)) {
    stop("no epi parameters defined. See ?define_epi_params")
  }
  
  # get project params into standardised format and perform checks
  args <- process_params(project$epi_parameters)
  check_params(args)
  
  # get complete demog from life table
  demog <- get_demography(project$epi_parameters$life_table)
  args <- c(args,
            list(age_death = demog$age_death,
                 age_stable = demog$age_stable))
  
  # append arguments
  args <- c(args,
            list(max_time = max_time,
                 save_transmission_record = save_transmission_record,
                 transmission_record_location = transmission_record_location,
                 output_daily_counts = output_daily_counts,
                 output_age_distributions = output_age_distributions,
                 output_age_times = output_age_times,
                 silent = silent))
  
  
  # ---------- run simulation ----------
  
  # internal flag, not visible to user. If TRUE then write parameter lists to
  # file and return without running simulation. Parameters can then be read
  # directly from file into Xcode.
  xcode_on <- FALSE
  if (xcode_on) {
    write_xcode_params(args)
    return()
  }
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # ---------- process output ----------
  
  return(output_raw)
}


#------------------------------------------------
# print parameters to file, to be read in directly from Xcode
#' @noRd
write_xcode_params <- function(args) {
  print(args)
  message("writing arguments to file")
  
  # functions for writing vectors and matrices to file
  vector_to_file <- function(x, file_path) {
    writeLines(paste(x, collapse = ","), con = file_path)
  }
  matrix_to_file <- function(x, file_path) {
    writeLines(mapply(function(y) paste(y, collapse = ","), x), con = file_path)
  }
  
  # write scalar epi parameters to file
  arg_file_path <- "R_ignore/SIMPLEGEN_Xcode/args/"
  scalar_epi_params <- c("a", "p", "mu", "u", "v", "g", "treatment_seeking_alpha",
                         "treatment_seeking_beta", "max_inoculations")
  vector_to_file(unlist(args[scalar_epi_params]), paste0(arg_file_path, "scalar_epi_params.txt"))
  
  print(args[scalar_epi_params])
  
  # write epi vectors to file
  vector_to_file(args$prob_infection, paste0(arg_file_path, "prob_infection.txt"))
  vector_to_file(args$prob_acute, paste0(arg_file_path, "prob_acute.txt"))
  vector_to_file(args$prob_AC, paste0(arg_file_path, "prob_AC.txt"))
  vector_to_file(args$duration_prophylactic, paste0(arg_file_path, "duration_prophylactic.txt"))
  
  # write epi matrices (or lists of vectors) to file
  matrix_to_file(args$duration_acute, paste0(arg_file_path, "duration_acute.txt"))
  matrix_to_file(args$duration_chronic, paste0(arg_file_path, "duration_chronic.txt"))
  matrix_to_file(args$detectability_microscopy_acute, paste0(arg_file_path, "detectability_microscopy_acute.txt"))
  matrix_to_file(args$detectability_microscopy_chronic, paste0(arg_file_path, "detectability_microscopy_chronic.txt"))
  matrix_to_file(args$detectability_PCR_acute, paste0(arg_file_path, "detectability_PCR_acute.txt"))
  matrix_to_file(args$detectability_PCR_chronic, paste0(arg_file_path, "detectability_PCR_chronic.txt"))
  matrix_to_file(args$time_treatment_acute, paste0(arg_file_path, "time_treatment_acute.txt"))
  matrix_to_file(args$time_treatment_chronic, paste0(arg_file_path, "time_treatment_chronic.txt"))
  vector_to_file(args$infectivity_acute, paste0(arg_file_path, "infectivity_acute.txt"))
  vector_to_file(args$infectivity_chronic, paste0(arg_file_path, "infectivity_chronic.txt"))
  
  # write deme vectors to file
  vector_to_file(args$H, paste0(arg_file_path, "H.txt"))
  vector_to_file(args$seed_infections, paste0(arg_file_path, "seed_infections.txt"))
  vector_to_file(args$M, paste0(arg_file_path, "M.txt"))
  
  # write demog vectors to file
  vector_to_file(args$life_table, paste0(arg_file_path, "life_table.txt"))
  vector_to_file(args$age_death, paste0(arg_file_path, "age_death.txt"))
  vector_to_file(args$age_stable, paste0(arg_file_path, "age_stable.txt"))
  
  # write scalar run parameters to file
  scalar_run_params <- c("max_time", "save_transmission_record", "transmission_record_location",
                         "output_daily_counts", "output_age_distributions", "silent")
  vector_to_file(args[scalar_run_params], paste0(arg_file_path, "scalar_run_params.txt"))
  
  # write vector run parameters to file
  vector_to_file(args$output_age_times, paste0(arg_file_path, "output_age_times.txt"))
  
}
