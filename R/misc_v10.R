
#------------------------------------------------
# if NULL then replace with chosen value, otherwise keep original value
#' @noRd
define_default <- function(x, default) {
  if (is.null(x)) {
    return(default)
  } else {
    return(x)
  }
}

#------------------------------------------------
# if a single value is provided then expand to a vector of length n
#' @noRd
force_vector <- function(x, n) {
  if (length(x) == 1) {
    return(rep(x,n))
  } else {
    return(x)
  }
}

#------------------------------------------------
# calculate midpoints of a vector
#' @noRd
midpoints <- function(x) {
  return((x[-1] + x[-length(x)])/2)
}

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
# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x - x_max)))
  return(ret)
}

#------------------------------------------------
# geweke_pvalue
# return p-value of Geweke's diagnostic convergence statistic, estimated from
# package coda
#' @importFrom stats pnorm
#' @importFrom coda geweke.diag
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(coda::geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check that geweke p-value non-significant at alpha significance level on
# values x[1:n]
#' @importFrom coda mcmc
#' @noRd
test_convergence <- function(x, n, alpha = 0.01) {
  # fail if n = 1
  if (n == 1) {
    return(FALSE)
  }
  
  # fail if ESS too small
  ESS <- try(coda::effectiveSize(x[1:n]), silent = TRUE)
  if (class(ESS) == "try-error") {
    return(FALSE)
  }
  if (ESS < 10) {
    return(FALSE)
  }
  
  # fail if geweke p-value < threshold
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g > alpha)
  if (is.na(ret)) {
    ret <- FALSE;
  }
  
  # return
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

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE.
#' @noRd
user_yes_no <- function(x = "continue? (Y/N): ") {
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}

#------------------------------------------------
# recursive function for converting nested list (of any depth) to single string.
# Each list is enclosed in "{}" each element within a list is separated with
# ",", and values in a vector are separated with " ".
#' @noRd
list_to_text <- function(x, s = NULL) {
  s <- paste0(s, "{")
  for (i in 1:length(x)) {
    if (i > 1) {
      s <- paste0(s, ",")
    }
    if (is.list(x[[i]])) {
      s <- list_to_text(x[[i]], s)
    } else {
      s <- paste0(s, paste(x[[i]], collapse = " "))
    }
  }
  s <- paste0(s, "}")
  return(s)
}

#------------------------------------------------
# write a list x to file. See list_to_text() for string format
#' @noRd
write_text_list <- function(x, file_path) {
  
  # convert list to single string
  s <- list_to_text(x)
  
  # write to file
  writeLines(s, file_path)
  
}

#-----------------------------
#' @title Retrieve prevalence
#' 
#' @Description Function to add detection layer to estimated prevalence
#'
#' @param data  Model output from SIMPLEGEN sim_epi() function epi_output$daily_values
#' @param case_detection Method of case detection, "Active" or "Passive"
#' @param diagnosis Method of diagnosis, "Microscopy" or "PCR"
#' @param sampling_time Sampling time, currently input should be a numeric day e.g. 100 would represent the 100th day of simulation.
#' @param sampled TO DO - sample from bernoulli distribution to generate simulated numbers of cases detected on each day - done in SIMPLEGEN internally already
#' @param deme deme of interest - currently assumes just a single deme analysed - i.e. can't pool results of multiple demes
#'
#' @importFrom stats aggregate
#' @return returns named vector :  c("annual_EIR", "prevalence_true", "prevalence_detected")
#' @export


retrieve_prev <-
  function(data,
           case_detection,
           diagnosis,
           sampling_time,
           sampled = FALSE,
           deme = 1) {
    results<-vector(length=3)
    names(results) <-
      c("annual_EIR", "prevalence_true", "prevalence_detected")
    
    # Calculating Annual EIR --------------------------------------------------
    
    # subset to desired rows and columns
    df_wide <-
      data[, c("time", "deme", "EIR")]
    df_wide <- df_wide[df_wide$deme == deme,]
    
    df_wide$year <-
      lubridate::year(lubridate::as_date(df_wide$time, origin = lubridate::origin))
    
    annual_EIR <-
      stats::aggregate(EIR ~ year + deme, data = df_wide, FUN = "sum")
    #annual_A<-aggregate(A~year + deme,data=df_wide, FUN= "mean" )
    
    sample_year<-lubridate::year(lubridate::as_date(sampling_time, origin = lubridate::origin))
    
    results["annual_EIR"] <- annual_EIR[which(annual_EIR$year==sample_year),"EIR"]
    
    
    # Calculating true prevalence ---------------------------------------------
    data<-data[data$time==sampling_time,]
    
    
    results["prevalence_true"]<- sum(data$C,data$A)/data$H
    
    # Calculating observed prevalence -----------------------------------------
    
    if(case_detection == "Active") {
      if (diagnosis == "PCR") {
        results["prevalence_detected"] <-
          sum(data$A_detectable_PCR, data$C_detectable_PCR) / data$H
      }
      else if (diagnosis == "Microscopy") {
        results["prevalence_detected"] <-
          sum(data$A_detectable_microscopy,data$C_detectable_microscopy) / data$H
        
      } else{
        warning("diagnosis type not recognized, please enter Microscopy or PCR")
      }
    } else if (case_detection == "Passive") {
      if (diagnosis == "PCR") {
        results["prevalence_detected"] <- data$A_detectable_PCR/data$H
      }
      else if (diagnosis == "Microscopy") {
        results["prevalence_detected"] <- data$A_detectable_microscopy/data$H
        
      } else {
        warning("diagnosis type not recognized, please enter Microscopy or PCR")
      }
      
    } else {
      warning("case_detection type not recognized, please enter Active or Passive")
    }
    
    return(results)
  }
