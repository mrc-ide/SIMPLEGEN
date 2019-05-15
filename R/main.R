
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
#' @title Define epidemiological parameters
#'
#' @description Define the epidemiological parameters of a SIMPLEGEN project.
#'   These parameters will be used when simulating epidemiolocal data.
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
#'   by an infected mosquito.
#' @param prob_acute probability an infection goes through an acute phase.
#' @param prob_AC probability of acute infection transitioning to chronic before
#'   clearing, as opposed to clearing directly.
#' @param duration_acute vector or list specifying probability distribution of
#'   time (in days) of acute phase of disease. If a list then the first element
#'   specifies the distribution for the first incident of acute disease, the
#'   second element for the second incident and so on (the final distribution is
#'   used for all remaining incidents). If a vector then the same distribution
#'   is used for all incidents of acute disease.
#' @param duration_chronic equivalent to \code{duration_acute} but for chronic
#'   phase of disease.
#' @param infectivity_acute probability a mosquito becomes infected after biting
#'   a human host in the acute phase.
#' @param infectivity_chronic probability a mosquito becomes infected after
#'   biting a human host in the chronic phase.
#' @param max_innoculations maximum number of innoculations that an individual
#'   can hold simultaneously.
#'
#' @importFrom stats dgeom
#' @export

define_epi_parameters <- function(project,
                                  a = 0.3,
                                  p = 0.85,
                                  mu = -log(p),
                                  u = 12,
                                  v = 10,
                                  g = 12,
                                  prob_infection = 0.6,
                                  prob_acute = c(1,0.5,0),
                                  prob_AC = 0.2,
                                  duration_acute = dgeom(1:25, 1/5),
                                  duration_chronic = dgeom(1:1000, 1/200),
                                  infectivity_acute = 0.07,
                                  infectivity_chronic = 0.07,
                                  max_innoculations = 5) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  assert_single_bounded(a)
  assert_single_bounded(p)
  assert_single_pos(mu)
  assert_single_pos_int(u, zero_allowed = FALSE)
  assert_single_pos_int(v, zero_allowed = FALSE)
  assert_single_pos_int(g, zero_allowed = FALSE)
  assert_bounded(prob_infection)
  assert_bounded(prob_acute)
  assert_single_bounded(prob_AC)
  if (!is.list(duration_acute)) {  # force to list
    duration_acute <- list(duration_acute)
  }
  mapply(assert_pos, duration_acute)
  if (!is.list(duration_chronic)) {  # force to list
    duration_chronic <- list(duration_chronic)
  }
  mapply(assert_pos, duration_chronic)
  assert_bounded(infectivity_acute)
  assert_bounded(infectivity_chronic)
  assert_single_pos_int(max_innoculations, zero_allowed = FALSE)
  
  # normalise input distributions to sum to 1
  for (i in 1:length(duration_acute)) {
    duration_acute[[i]] <- duration_acute[[i]]/sum(duration_acute[[i]])
  }
  for (i in 1:length(duration_chronic)) {
    duration_chronic[[i]] <- duration_chronic[[i]]/sum(duration_chronic[[i]])
  }
  
  # modify project to store parameters
  project$sim_parameters$epi_parameters <- list(a = a,
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
                                                infectivity_acute = infectivity_acute,
                                                infectivity_chronic = infectivity_chronic,
                                                max_innoculations = max_innoculations)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define deme parameters
#'
#' @description Define the parameters of each deme in a SIMPLEGEN simulation.
#'   These parameters will be used when simulating epidemiolocal data.
#'
#' @param project a SIMPLEGEN project.
#' @param H vector specifying human population size in each deme.
#' @param seed_infections vector specifying the initial number of infected
#'   humans in each deme. Infected humans are assumed to have just been bitten,
#'   meaning they have just entered the latent phase.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#'
#' @export

define_deme_parameters <- function(project,
                                   H = 1000,
                                   seed_infections = 100,
                                   M = 1000) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  assert_pos_int(H)
  assert_pos_int(seed_infections)
  assert_pos_int(M)
  assert_same_length_multiple(H, seed_infections, M)
  assert_leq(seed_infections, H)
  
  # modify project to store parameters
  project$sim_parameters$deme_parameters <- list(H = H,
                                                 seed_infections = seed_infections,
                                                 M = M)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define demographic parameters
#'
#' @description Define demography used in SIMPLEGEN simulation. The raw input is
#'   a life table giving the probability of death in each one-year age group.
#'   This is used to derive the distribution of age of death, and the stable age
#'   distribution.
#'
#' @param project a SIMPLEGEN project.
#' @param life_table vector specifying probability of death in each one-year age
#'   group. Final value must be 1 to ensure a closed population. If \code{NULL}
#'   then an age-distribution from Mali is used by default.
#'
#' @export

define_demograpy <- function(project,
                             life_table = NULL) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  if (is.null(life_table)) {
    life_table_raw <- simplegen_file("Mali_life_table.csv")
    life_table <- life_table_raw[,2]
  }
  assert_bounded(life_table)
  assert_eq(life_table[length(life_table)], 1, message = "the final value in the life table must be 1, representing a 100%% chance of dying, to ensure a closed population")
  
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
  
  # modify project to store distributions
  project$sim_parameters$demography <- list(life_table = life_table,
                                            age_death = age_death,
                                            age_stable = age_stable)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define migration parameters
#'
#' @description Define default migration parameters for SIMPLEGEN simulation.
#'
#' @param project a SIMPLEGEN project.
#' @param migration_matrix TODO - rethink this input type to allow "trip"
#'   migration.
#'
#' @export

define_migration <- function(project,
                             migration_matrix = matrix(1)) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  
  # modify project
  project$sim_parameters$migration <- migration_matrix
  
  # return
  invisible(project)
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
  
  # parameters from project
  args <- project$sim_parameters
  
  # simulation parameters from function arguments
  args$run_parameters <- list(max_time = max_time,
                              output_daily_counts = output_daily_counts,
                              output_age_distributions = output_age_distributions,
                              output_age_times = output_age_times,
                              output_infection_history = output_infection_history,
                              silent = silent)
  
  # ---------- run simulation ----------
  
  # internal flag, not visible to user. If TRUE then write parameter lists to
  # file and return without running simulation. Parameters can then be read
  # directly from file into Xcode.
  xcode_on <- FALSE
  if (xcode_on) {
    message("writing arguments to file")
    
    # write scalar epi parameters to file
    scalar_params <- c("a", "p", "mu", "u", "v", "g", "prob_AC", "max_innoculations")
    writeLines(paste(unlist(args$epi_parameters[scalar_params]), collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/scalar.txt")
    
    # write epi vectors to file
    writeLines(paste(args$epi_parameters$prob_infection, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/prob_infection.txt")
    writeLines(paste(args$epi_parameters$prob_acute, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/prob_acute.txt")
    writeLines(paste(args$epi_parameters$infectivity_acute, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/infectivity_acute.txt")
    writeLines(paste(args$epi_parameters$infectivity_chronic, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/infectivity_chronic.txt")
    
    # write epi matrices (or lists of vectors) to file
    writeLines(mapply(function(x) paste(x, collapse = ","), args$epi_parameters$duration_acute),
               con = "R_ignore/SIMPLEGEN_Xcode/args/duration_acute.txt")
    writeLines(mapply(function(x) paste(x, collapse = ","), args$epi_parameters$duration_chronic),
               con = "R_ignore/SIMPLEGEN_Xcode/args/duration_chronic.txt")
    
    # write deme vectors to file
    writeLines(paste(args$deme_parameters$H, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/H.txt")
    writeLines(paste(args$deme_parameters$seed_infections, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/seed_infections.txt")
    writeLines(paste(args$deme_parameters$M, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/M.txt")
    
    # write demog vectors to file
    writeLines(paste(args$demography$life_table, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/life_table.txt")
    writeLines(paste(args$demography$age_death, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/age_death.txt")
    writeLines(paste(args$demography$age_stable, collapse = ","),
               con = "R_ignore/SIMPLEGEN_Xcode/args/age_stable.txt")
    
    return()
  }
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # ---------- process output ----------
  
  return(output_raw)
}
