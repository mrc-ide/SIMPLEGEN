
# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

# -----------------------------------
# takes list format returned from Rcpp and converts to three-dimensional array.
# Array indexing is in the same order as the underlying list, for example
# x[i,j,k] is equivalent to l[[i]][[j]][[k]]
#' @noRd
rcpp_to_array <- function(x) {
  ret <- array(unlist(x), dim = c(length(x[[1]][[1]]), length(x[[1]]), length(x)))
  ret <- aperm(ret, perm = c(3,2,1))
  return(ret)
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
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
# does interval x intersect interval y
# both integer vectors of length 2
#' @noRd
interval_intersect <- function(x, y) {
  (x[1] < y[2]) & (x[2] > y[1])
}

#-----------------------------
#' @title Retrieve prevalence
#' 
#' @Description Function to add detection layer to estimated prevalence
#'
#' @param data Model output from SIMPLEGEN sim_epi() function epi_output$daily_values
#' @param case_detection Method of case detection, "Active" or "Passive"
#' @param diagnosis Method of diagnosis, "Microscopy" or "PCR"
#' @param sampling_time Sampling time, currently input should be a numeric day e.g. 100 would represent the 100th day of simulation.
#' @param sampled TO DO - sample from bernoulli distribution to generate simulated numbers of cases detected on each day - done in SIMPLEGEN internally already
#' @param deme deme of interest - currently assumes just a single deme analysed - i.e. can't pool results of multiple demes
#'
#' @importFrom stats aggregate
#' @return returns named vector :  c("annual_EIR", "prevalence_true", "prevalence_detected")
#' @export


retrieve_prev <- function(data,
                          case_detection,
                          diagnosis,
                          sampling_time,
                          sampled = FALSE,
                          deme = 1) {
  
  results <- vector(length=3)
  names(results) <- c("annual_EIR", "prevalence_true", "prevalence_detected")
  
  # Calculating Annual EIR --------------------------------------------------
  
  # subset to desired rows and columns
  df_wide <- data[, c("time", "deme", "EIR")]
  df_wide <- df_wide[df_wide$deme == deme,]
  
  df_wide$year <- lubridate::year(lubridate::as_date(df_wide$time, origin = lubridate::origin))
  
  annual_EIR <- stats::aggregate(EIR ~ year + deme, data = df_wide, FUN = "sum")
  
  sample_year <- lubridate::year(lubridate::as_date(sampling_time, origin = lubridate::origin))
  
  results["annual_EIR"] <- annual_EIR[which(annual_EIR$year == sample_year),"EIR"]
  
  
  # Calculating true prevalence ---------------------------------------------
  data <- data[data$time == sampling_time,]
  
  results["prevalence_true"] <- sum(data$C,data$A) / data$H
  
  # Calculating observed prevalence -----------------------------------------
  
  if (case_detection == "Active") {
    if (diagnosis == "PCR") {
      results["prevalence_detected"] <- sum(data$A_detectable_PCR, data$C_detectable_PCR) / data$H
    } else if (diagnosis == "Microscopy") {
      results["prevalence_detected"] <-sum(data$A_detectable_microscopy,data$C_detectable_microscopy) / data$H
    } else {
      warning("diagnosis type not recognized, please enter Microscopy or PCR")
    }
  } else if (case_detection == "Passive") {
    if (diagnosis == "PCR") {
      results["prevalence_detected"] <- data$A_detectable_PCR / data$H
    } else if (diagnosis == "Microscopy") {
      results["prevalence_detected"] <- data$A_detectable_microscopy / data$H
    } else {
      warning("diagnosis type not recognized, please enter Microscopy or PCR")
    }
  } else {
    warning("case_detection type not recognized, please enter Active or Passive")
  }
  
  return(results)
}
