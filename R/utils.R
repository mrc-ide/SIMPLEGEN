#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#------------------------------------------------
# force objects that can be defined as list or vector into a list
force_list <- function(x) {
  if (!is.list(x)) {
    x <- list(x)
  }
  return(x)
}

#------------------------------------------------
# force vectors to a given length by replicating scalar values if needed
force_veclength <- function(x, n) {
  if (length(x) == 1) {
    x <- rep(x, n)
  }
  return(x)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if ((i == max_i) & close) {
    close(pb_list[[name]])
  }
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
#'   bands. TODO - replace with new demography table that I can properly
#'   reference!
#'
#' @export

life_table_Mali <- function() {
  life_table_raw <- simplegen_file("Mali_life_table.csv")
  ret <- life_table_raw$prop_death
  return(ret)
}

#------------------------------------------------
# get stable demography distribution from life table via Eigenvalues
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
  m[col(m) == (row(m) + 1)] <- 1 - life_table[1:(n - 1)]
  m[, 1] <- 1 - rowSums(m)
  
  # convert to rates
  r = m - diag(n)
  
  # compute Eigenvalues of the rate matrix
  E = eigen(t(r))
  
  # there should be one Eigenvalue that is zero (up to limit of computational
  # precision). Find which Eigenvalue this is
  w <- which.min(abs(E$values))
  
  # the stable solution is the corresponding Eigenvector, suitably normalised
  age_stable <- Re(E$vectors[, w] / sum(E$vectors[, w]))
  
  # return list of distributions
  ret <- list(life_table = life_table,
              age_death = age_death,
              age_stable = age_stable)
  return(ret)
}

#------------------------------------------------
# convert epi model parameters to standardised types - for example, deals with
# arguments that can be defined as either list or vector. Ensures that values
# passed to downstream functions are unambiguous. NB. check_epi_model_params()
# is run before process_epi_model_params()
#' @noRd
process_epi_model_params <- function(project) {
  
  # cpp11 requires integers to be defined as distinct from doubles
  name_vec <- c("H", "seed_infections", "M")
  project$epi_model_parameters[name_vec] <- mapply(as.integer, project$epi_model_parameters[name_vec], SIMPLIFY = FALSE)
  
  # force duration distributions and daily probabilities to lists
  name_vec <- c("duration_acute", "duration_chronic",
                "time_treatment_acute", "time_treatment_chronic",
                "duration_prophylactic",
                "detectability_microscopy_acute", "detectability_microscopy_chronic",
                "detectability_PCR_acute", "detectability_PCR_chronic",
                "infectivity_acute", "infectivity_chronic")
  project$epi_model_parameters[name_vec] <- mapply(force_list, project$epi_model_parameters[name_vec], SIMPLIFY = FALSE)
  
  # force deme properties to vectors over demes
  name_vec <- c("H", "M", "seed_infections")
  n_demes <- nrow(project$epi_model_parameters$mig_mat)
  project$epi_model_parameters[name_vec] <- mapply(force_veclength, project$epi_model_parameters[name_vec],
                                                   MoreArgs = list(n = n_demes), SIMPLIFY = FALSE)
  
  # return
  invisible(project)
}

#------------------------------------------------
# process sampling data.frames to make compatible with cpp11 format. For daily
# and sweep outputs this involves subsetting columns, converting everything to
# numeric, and adding a time column if not present
#' @noRd
process_sweep <- function(sweep_df) {
  
  # subset columns
  ret <- sweep_df %>%
    dplyr::select(measure, state, diagnostic, deme, age_min, age_max)
  
  # add dummy time column if not present
  if (!("time" %in% names(ret))) {
    ret$time <- -1
  }
  
  # convert all strings to numeric starting at 0
  ret <- ret %>%
    dplyr::mutate(measure = match(measure, c("count", "prevalence", "incidence_active", "incidence_passive", "EIR")) - 1,
                  state = match(state, c("S", "E", "A", "C", "P", "H")) - 1,
                  diagnostic = match(diagnostic, c("True", "Microscopy", "PCR")) - 1)
  
  # start deme index at 0
  ret$deme[ret$deme != -1] <- ret$deme[ret$deme != -1] - 1
  
  # convert all NAs to -1
  ret$state[is.na(ret$state)] <- -1
  ret$diagnostic[is.na(ret$diagnostic)] <- -1
  
  # convert numeric values to integer type (cpp11 makes this distinction)
  ret <- mapply(as.integer, ret) %>%
    as.data.frame()
  
  return(ret)
}
