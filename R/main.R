
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
#' @title Define SIMPLEGEN parameters
#'
#' @description Define the parameters of a SIMPLEGEN project. These parameters
#'   will be used when simulating epidemiolocal data.
#'
#' @param project a SIMPLEGEN project.
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
#'   by an infected mosquito. If a vector then each value applies to each
#'   subsequent infection.
#' @param prob_acute probability an infection passes to the acute stage, rather
#'   than going directly to the chronic stage. If a vector then each value
#'   applies to each subsequent infection.
#' @param prob_AC probability an infection in the acute stage passes to the
#'   chronic stage prior to recovery, rather than recovering directly. If a
#'   vector then each value applies to each subsequent infection.
#' @param duration_acute specifies the probability distribution of duration (in
#'   days) of acute disease. Can be vector or list: if a list then the first
#'   element specifies the probability distribution for the first incident of
#'   acute disease, the second element for the second incident and so on (the
#'   final element is used for all remaining incidents); if a vector then the
#'   same probability distribution is used for all incidents of acute disease.
#'   Values are automatically normalised to a proper probability mass
#'   distribution internally (i.e. a distribution that sums to 1).
#' @param duration_chronic equivalent to \code{duration_acute} but for chronic
#'   disease.
#' @param time_treatment_acute specifies the probability distribution of time
#'   (in days since entering the acute stage) until considering seeking
#'   treatment. If the infection clears naturally prior to this time then the
#'   host does not seek treatment. Otherwise the host seeks treatment according
#'   to their host-specific treatment seeking parameter (see
#'   \code{treatment_seeking_alpha} and \code{treatment_seeking_beta}). Can be
#'   vector or list: if a list then the first element specifies the probability
#'   distribution for the first incident of acute disease, the second element
#'   for the second incident and so on (the final element is used for all
#'   remaining incidents); if a vector then the same probability distribution is
#'   used for all incidents of acute disease. Values are automatically
#'   normalised to a proper probability mass distribution internally (i.e. a
#'   distribution that sums to 1).
#' @param time_treatment_chronic equivalent to \code{time_treatment_acute} but
#'   for chronic disease.
#' @param infectivity_acute the probability that a mosquito becomes infected
#'   upon biting an acutely infective human host, at each day since entering the
#'   infective state. Can be vector or list: if a list then the first element
#'   specifies the infectivity for the first incident of acute disease, the
#'   second element for the second incident and so on (the final element is used
#'   for all remaining incidents); if a vector then the same infectivity is used
#'   for all incidents of acute disease.
#' @param infectivity_chronic equivalent to \code{infectivity_acute} but for
#'   chronic disease.
#' @param max_inoculations maximum number of inoculations that an individual
#'   can hold simultaneously.
#' @param H vector specifying human population size in each deme.
#' @param seed_infections vector specifying the initial number of infected
#'   humans in each deme. Infected humans are assumed to have just been bitten,
#'   meaning they have just entered the latent phase.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#' @param life_table vector specifying probability of death in each one-year age
#'   group. Final value must be 1 to ensure a closed population. If \code{NULL}
#'   then an age distribution from Mali is used by default.
#'
#' @importFrom stats dgeom
#' @export

define_params <- function(project,
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
                          time_treatment_acute = dgeom(1:25, 1/5),
                          time_treatment_chronic = dgeom(1:100, 1/20),
                          infectivity_acute = 0.07,
                          infectivity_chronic = 0.07,
                          max_inoculations = 5,
                          H = 1000,
                          seed_infections = 100,
                          M = 1000,
                          life_table = NULL) {
  
  # get current life table
  life_table_before <- project$life_table
  
  # use Mali life table by default
  if (is.null(life_table)) {
    life_table_raw <- simplegen_file("Mali_life_table.csv")
    life_table <- life_table_raw[,2]
  }
  
  # replace project parameters with user-defined
  userlist <- as.list(match.call())[-(1:2)]
  project$sim_parameters[names(userlist)] <- mapply(eval, userlist, SIMPLIFY = FALSE)
  
  # give all NULL parameters default values
  arglist <- list(a = a,
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
                  time_treatment_acute = time_treatment_acute,
                  time_treatment_chronic = time_treatment_chronic,
                  infectivity_acute = infectivity_acute,
                  infectivity_chronic = infectivity_chronic,
                  max_inoculations = max_inoculations,
                  H = H,
                  seed_infections = seed_infections,
                  M = M,
                  life_table = life_table)
  argnull <- mapply(is.null, project$sim_parameters[names(arglist)])
  project$sim_parameters[which(argnull)] <- arglist[which(argnull)]
  
  # standardise and perform checks
  params_processed <- process_params(project$sim_parameters)
  check_params(params_processed)
  
  # update demography distributions as needed
  if (!isTRUE(all.equal(life_table_before, project$sim_parameters$life_table))) {
    demog <- get_demography(life_table)
    project$sim_parameters$age_death <- demog$age_death
    project$sim_parameters$age_stable <- demog$age_stable
  }
  
  # return
  invisible(project)
}

#------------------------------------------------
# convert parameters to standardised type
#' @noRd
process_params <- function(x) {
  
  # standardise parameters
  if (!is.list(x$duration_acute)) {
    x$duration_acute <- list(x$duration_acute)
  }
  if (!is.list(x$duration_chronic)) {
    x$duration_chronic <- list(x$duration_chronic)
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
  mapply(assert_bounded, x$time_treatment_acute, name = "time_treatment_acute")
  mapply(assert_bounded, x$time_treatment_chronic, name = "time_treatment_chronic")
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
#' @description Simulate data from a simple individual-based model of humans and
#'   mosquitoes. Parameters are taken from the values stored within the
#'   \code{sim_parameters} slot of the project, and outputs are written to the
#'   \code{sim_output} slot.
#'
#' @param project a SIMPLEGEN project.
#' @param max_time run simulation for this many days.
#' @param output_daily_counts whether to output daily counts of key quantities,
#'   such as the number of infected hosts and the EIR.
#' @param output_age_distributions whether to output complete age distributions
#' @param output_age_times a vector of times at which complete age distributions
#'   are output.
#' @param output_infection_history whether to output complete infection history.
#' @param silent whether to write messages to console.
#'
#' @export

sim_epi <- function(project,
                    max_time = 365,
                    output_daily_counts = TRUE,
                    output_age_distributions = TRUE,
                    output_age_times = max_time,
                    output_infection_history = FALSE,
                    silent = FALSE) {
  
  # ---------- check inputs ----------
  
  assert_custom_class(project, "simplegen_project")
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_single_logical(output_daily_counts)
  assert_single_logical(output_age_distributions)
  assert_vector(output_age_times)
  assert_pos_int(output_age_times)
  assert_leq(output_age_times, max_time)
  assert_single_logical(output_infection_history)
  assert_single_logical(silent)
  
  # ---------- define argument lists ----------
  
  # get project params into standardised format and perform checks
  args <- process_params(project$sim_parameters)
  check_params(args)
  
  # append arguments
  args <- c(args,
            list(max_time = max_time,
                 output_daily_counts = output_daily_counts,
                 output_age_distributions = output_age_distributions,
                 output_age_times = output_age_times,
                 output_infection_history = output_infection_history,
                 silent = silent))
  
  # ---------- run simulation ----------
  
  # internal flag, not visible to user. If TRUE then write parameter lists to
  # file and return without running simulation. Parameters can then be read
  # directly from file into Xcode.
  xcode_on <- FALSE
  if (xcode_on) {
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
    scalar_epi_params <- c("a", "p", "mu", "u", "v", "g", "max_inoculations")
    vector_to_file(unlist(args[scalar_epi_params]), paste0(arg_file_path, "epi_scalar.txt"))
    
    # write epi vectors to file
    vector_to_file(args$prob_infection, paste0(arg_file_path, "prob_infection.txt"))
    vector_to_file(args$infectivity, paste0(arg_file_path, "infectivity.txt"))
    
    # write epi matrices (or lists of vectors) to file
    matrix_to_file(args$duration_disease, paste0(arg_file_path, "duration_disease.txt"))
    
    # write deme vectors to file
    vector_to_file(args$deme_parameters$H, paste0(arg_file_path, "H.txt"))
    vector_to_file(args$deme_parameters$seed_infections, paste0(arg_file_path, "seed_infections.txt"))
    vector_to_file(args$deme_parameters$M, paste0(arg_file_path, "M.txt"))
    
    # write demog vectors to file
    vector_to_file(args$demography$life_table, paste0(arg_file_path, "life_table.txt"))
    vector_to_file(args$demography$age_death, paste0(arg_file_path, "age_death.txt"))
    vector_to_file(args$demography$age_stable, paste0(arg_file_path, "age_stable.txt"))
    
    # write scalar run parameters to file
    scalar_run_params <- c("max_time", "output_daily_counts", "output_age_distributions",
                           "output_infection_history", "silent")
    vector_to_file(args$run_parameters[scalar_run_params], paste0(arg_file_path, "run_scalar.txt"))
    
    # write run vectors to file
    vector_to_file(args$run_parameters$output_age_times, paste0(arg_file_path, "output_age_times.txt"))
    
    return()
  }
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # ---------- process output ----------
  
  return(output_raw)
}
