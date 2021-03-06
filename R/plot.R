
#------------------------------------------------
# a series of internally-used colours
#' @noRd
simplegen_cols <- function() {
  c("firebrick1", "chartreuse3", "dodgerblue", "dodgerblue4", "purple", "darkorange", "firebrick4")
}

#------------------------------------------------
# change the alpha value (transparency) of any named colour
#' @importFrom grDevices col2rgb rgb
#' @noRd
set_col_alpha <- function(col, alpha) {
  apply(col2rgb(col), 2, function(x) rgb(x[1], x[2], x[3], alpha*255, maxColorValue = 255))
}

#------------------------------------------------
#' @title Generic function for plotting a SIMPLEGEN flexible distribution
#'
#' @description Many distributions within the inbuilt SIMPLEGEN transmission
#'   model can be defined flexibly, as either a vector or a list of vectors.
#'   This function facilitates quick visualisation of these distributions.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param name the name of the distribution to be plotted, as it is found within
#'   the \code{project$epi_model_parameters}
#'
#' @importFrom stats dbeta
#' @export

plot_epi_distribution <- function(project, name = "duration_acute") {
  
  # define types of allowed distribution
  duration_dist <- c("duration_acute", "duration_chronic", "time_treatment_acute",
                     "time_treatment_chronic", "duration_prophylactic")
  transition_prob <- c("prob_infection", "prob_acute", "prob_AC")
  daily_prob <- c("detectability_microscopy_acute", "detectability_microscopy_chronic",
                  "detectability_PCR_acute", "detectability_PCR_chronic", "infectivity_acute",
                  "infectivity_chronic")
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_in(name, c(duration_dist, transition_prob, daily_prob))
  
  # switch between plotting functions for durations vs. transition probabilities
  # vs. daily probabilities
  if (name %in% duration_dist) {
    ret <- plot_epi_distribution_duration(project, name, normalise = TRUE)
  } else if (name %in% transition_prob) {
    ret <- plot_epi_distribution_transition(project, name)
  } else {
    ret <- plot_epi_distribution_duration(project, name)
  }
  
  return(ret)
}

#------------------------------------------------
# plot duration distribution
#' @noRd
plot_epi_distribution_duration <- function(project, name, normalise = FALSE) {
  
  # required to remove notes
  x <- NULL
  
  # check epi model parameters
  check_epi_model_params(project)
  
  # make plotting dataframe
  y <- project$epi_model_parameters[[name]]
  plot_df <- do.call(rbind, mapply(function(i) {
    ret <- data.frame(x = seq_along(y[[i]]), y = y[[i]], rep = i)
    if (normalise) {
      ret$y <- ret$y/sum(ret$y)
    }
    return(ret)
  }, seq_along(y), SIMPLIFY = FALSE))
  
  # produce plot
  ret <- ggplot2::ggplot(data = plot_df) + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank()) +
    ggplot2::geom_bar(ggplot2::aes(x = x, y = y), color = NA, fill = simplegen_cols()[4],
                      stat = "identity", width = 0.9) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05*max(plot_df$y)), expand = c(0, 0)) +
    ggplot2::xlab("day") + ggplot2::ylab("probability") +
    ggplot2::ggtitle(name)
  
  # facet if multiple distributions
  if (length(y) > 1) {
    ret <- ret + ggplot2::facet_wrap(. ~rep) +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = grey(1.0)))
  }
  
  return(ret)
}

#------------------------------------------------
# plot transition probabilities
#' @noRd
plot_epi_distribution_transition <- function(project, name) {
  
  # required to remove notes
  x <- NULL
  
  # check epi model parameters
  check_epi_model_params(project)
  
  # make plotting dataframe
  y <- project$epi_model_parameters[[name]]
  plot_df <- data.frame(x = seq_along(y), y = y)
  
  # produce plot
  ret <- ggplot2::ggplot(data = plot_df) + ggplot2::theme_bw() +
    ggplot2::geom_bar(ggplot2::aes(x = as.factor(x), y = y), color = NA, fill = simplegen_cols()[4],
                      stat = "identity", width = 0.9) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05*max(plot_df$y)), expand = c(0, 0)) +
    ggplot2::xlab("exposure") + ggplot2::ylab("probability") +
    ggplot2::ggtitle(name)
  
  return(ret)
}

#------------------------------------------------
#' @title Plot distribution of treatment seeking in the population
#'
#' @description Within the inbuilt SIMPLEGEN transmission model, each host has
#'   its own treatment seeking parameter defining the probability of
#'   actively seeking treatment once disease reaches a point that this becomes
#'   necessary. These parameter values are drawn independently for each host
#'   from a Beta distribution with user-defined mean and standard deviation.
#'   This function plots the Beta distribution for the current project values.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#'
#' @importFrom stats dbeta
#' @export

plot_treatment_seeking <- function(project) {
  
  # check inputs
  assert_class(project, "simplegen_project")
  
  # check epi model parameters
  check_epi_model_params(project)
  
  # extract treatment seeking parameters
  mu <- project$epi_model_parameters$treatment_seeking_mean
  sigma <- project$epi_model_parameters$treatment_seeking_sd
  alpha <- mu^2*(1 - mu)/sigma^2 - mu
  beta <- mu*(1 - mu)^2/sigma^2 - (1 - mu)
  
  # create plotting dataframe
  df_plot <- data.frame(x = seq(0, 1, 0.001))
  
  # define either Beta or Dirac delta distribution
  if  (mu == 0 | mu == 1 | sigma == 0) {
    df_plot$y <- 0
    w <- which.min(abs(df_plot$x - mu))
    df_plot$y[w] <- 1
  } else {
    df_plot$y <- dbeta(df_plot$x, shape1 = alpha, shape2 = beta)
  }
  
  # produce plot
  ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_area(ggplot2::aes_(x = ~x, y = ~y), colour = "black", fill = simplegen_cols()[4]) +
    ggplot2::geom_vline(xintercept = mu, linetype = "dashed") +
    ggplot2::xlab("Pr(seek treatment)") + ggplot2::ylab("probability density") +
    ggplot2::ggtitle("Treatment seeking distribution")
  
}

#------------------------------------------------
#' @title Plot daily prevalence
#'
#' @description Plot the daily prevalence in a series of defined model states.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param state,diagnostic,age_min,age_max,inoculations which rows of daily output to subset to.
#'
#' @importFrom grDevices grey
#' @import tidyr
#' @export

plot_daily_prevalence <- function(project,
                                  state = c("S", "E", "A", "C", "P"),
                                  diagnostic = "true",
                                  age_min = 0,
                                  age_max = 100,
                                  inoculations = -1) {
  
  # needed to avoid error "no visible global variable"
  #state <- count <- NULL
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_in(state, c("S", "E", "A", "C", "P"))
  
  # check that daily output exists
  assert_non_null(project$epi_output$daily)
  
  # subset to prevalence in specified states
  df_plot <- project$epi_output$daily
  df_plot <- df_plot[(df_plot$measure == "prevalence") &
                       (df_plot$state %in% state) &
                       (df_plot$diagnostic == diagnostic) &
                       (df_plot$age_min == age_min) &
                       (df_plot$age_max == age_max) &
                       (df_plot$inoculations == inoculations),]
  
  # error if nothing to plot
  if (nrow(df_plot) == 0) {
    stop("no daily prevalence values matching this combination")
  }
  
  # choose x-axis scale
  max_time <- max(df_plot$time)
  if (max_time <= 30) {
    x_breaks <- seq(0, max(df_plot$time))
    x_labels <-seq_along(x_breaks) - 1
    x_lab <- "time (days)"
  } else if (max_time <= 365) {
    x_breaks <- seq(0, max(df_plot$time), 365/12)
    x_labels <-seq_along(x_breaks) - 1
    x_lab <- "time (months)"
  } else {
    x_breaks <- seq(0, max(df_plot$time), 365)
    x_labels <-seq_along(x_breaks) - 1
    x_lab <- "time (years)"
  }
  
  # produce plot
  ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes_(x = ~time, y = ~value, color = ~state)) +
    ggplot2::scale_x_continuous(x_lab, breaks = x_breaks, labels = x_labels) +
    ggplot2::facet_wrap(~deme)
  
}

#------------------------------------------------
#' @title Plot age distribution of states
#'
#' @description If project contains age distributions in the output from running
#'   the transmission model, then this function can be used to produce a simple
#'   barplot for a specified state.
#'   
#'   TODO - this function is now depricated.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param sample_time which time-orderd sample to plot.
#' @param deme which deme to plot.
#' @param state which host state to plot on the y-axis.
#'
#' @export

plot_age_states <- function(project, sample_time = NULL, deme = 1, state = "S") {
  
  # check inputs
  assert_non_null(project$epi_output$age_distributions, message = "no age distribution output found")
  if (is.null(sample_time)) {
    sample_time <- project$epi_output$age_distributions$sample_time[1]
  }
  assert_single_pos_int(sample_time)
  assert_in(sample_time, project$epi_output$age_distributions$sample_time,
            message = "sample_time values not found in age distribution output")
  assert_single_pos_int(deme)
  assert_in(deme, project$epi_output$age_distributions$deme,
            message = "deme values not found in age distribution output")
  assert_single_string(state)
  assert_in(state, names(project$epi_output$age_distributions),
            message = "state not found in age distribution output")
  
  # subset data
  dat <- project$epi_output$age_distributions
  dat <- dat[(dat$sample_time == sample_time) & (dat$deme == deme), c("age", state)]
  names(dat) <- c("age", "y")
  
  # produce plot
  ggplot2::ggplot(dat) + ggplot2::theme_bw() +
    ggplot2::geom_bar(ggplot2::aes_(x = ~age, y = ~y), stat = "identity", width = 1) +
    ggplot2::xlab("age (years)") + ggplot2::ylab("prevalence") +
    ggplot2::ggtitle(sprintf("Age-distribution of state %s", state))
  
}


#------------------------------------------------
#' @title Plot EIR versus prevalence
#'
#' @description Plots relationship between annual EIR and prevalence
#'
#' @param data dataset produced by the retrieve_prev function
#' @param plot_studies logical argument whether to plot data from previous
#'   empirical studies (currently Hay et al., 2005)
#' @param scale_x linear or log scale for x axis
#'
#' @importFrom ggsci scale_color_lancet
#' @importFrom rlang .data
#' @export
#' 
plot_EIR_prevalence <- function(data, plot_studies = TRUE, scale_x = "linear"){
  
  # load data (this is apparently not good practice!)
  #data("EIRprev_hay2005")
  #data("EIRprev_beier1999")
  
  EIRprev <- tidyr::gather(as.data.frame(EIRprev), key = "detection_type", value = "prevalence", -1)
  
  # remove zeros to avoid issues when logging
  EIRprev <- EIRprev[is.finite(log(EIRprev[,"annual_EIR"])),]
  
  if (plot_studies) {
    EIRprev_hay2005$detection_type <- rep("Hay et al., 2005", length.out = nrow( EIRprev_hay2005))
    EIRprev_beier1999$detection_type <- rep("Beier et al., 1999", length.out = nrow( EIRprev_beier1999))
    EIRprev <- rbind( EIRprev, EIRprev_hay2005)
    EIRprev <- rbind(EIRprev,EIRprev_beier1999)
  }
  if (scale_x == "linear") {
    p <- ggplot2::ggplot(data = as.data.frame(EIRprev), ggplot2::aes(.data$annual_EIR, .data$prevalence, colour = .data$detection_type)) +
      ggplot2::theme_bw() + ggsci::scale_color_lancet() +
      ggplot2::geom_point() +
      ggplot2::labs(title = "Annual EIR and Prevalence", x = "Annual EIR", y = "Prevalence") +
      ggplot2::coord_cartesian(ylim = c(0,1), xlim = c(0,400))
  } else if (scale_x == "log") {
    p <- ggplot2::ggplot(as.data.frame(EIRprev), ggplot2::aes(.data$annual_EIR, .data$prevalence, colour = .data$detection_type)) +
      ggplot2::theme_bw() +
      ggplot2::geom_point() + ggplot2::scale_x_log10() + ggsci::scale_color_lancet() +
      ggplot2::labs(title = "Annual EIR (log scale) and Prevalence", x = "Annual EIR", y = "Prevalence") +
      ggplot2::coord_cartesian(ylim = c(0,1))
  } else {
    warning("scale_x must be log or linear")
  }
  print(p)
}

