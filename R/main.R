
#------------------------------------------------
#' @title Dummy function
#'
#' @description Simple test function that demonstrates some of the features of
#'   this package.
#'
#' @details Takes a vector of values, returns the square.
#'
#' @param x vector of values
#'
#' @export
#' @examples
#' # Find square of first 100 values
#' dummy1(1:100)

dummy1 <- function(x = 1:5) {
  
  # print message to console
  message("running R dummy1 function")
  
  # get arguments in list form
  args <- list(x = x)
  
  # run C++ function with these arguments
  output_raw <- dummy1_cpp(args)
  
  # some optional processing of output
  message("processing output")
  ret <- output_raw$x_squared
  
  # return
  return(ret)
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
