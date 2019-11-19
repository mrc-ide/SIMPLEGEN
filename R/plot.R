
#------------------------------------------------
#' @title Plot daily counts of each host state
#'
#' @description TODO.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param deme which deme to plot.
#' @param states which states to plot. Can be any subset of \code{c("S", "E",
#'   "A", "C", "P", "Sv", "Ev", "Iv")}.
#'
#' @importFrom grDevices grey
#' @export

plot_daily_states <- function(project, deme = 1, states = c("S", "E", "A", "C", "P")) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  assert_dataframe(project$epi_output$daily_values)
  assert_in(deme, project$epi_output$daily_values$deme)
  assert_in(states, c("S", "E", "A", "C", "P", "Sv", "Ev", "Iv"))
  
  # subset to desired rows and columns
  df_wide <- project$epi_output$daily_values[, c("time", "deme", states)]
  df_wide <- df_wide[df_wide$deme == deme,]
  
  # get to long format
  state <- count <- NULL  # needed to avoid error "no visible global variable"
  df_long <- df_wide
  df_long <- tidyr::gather(df_wide, state, count, states, factor_key = TRUE)
  
  # choose plotting colours
  plot_cols <- c(S = grey(0.5), E = "chartreuse3", A = "firebrick1", C = "dodgerblue", P = "purple",
                 Sv = grey(0.8), Ev = "darkorange", Iv = "firebrick4",
                 EIR = "black")
  
  # choose x-axis scale
  max_time <- max(df_wide$time)
  if (max_time <= 30) {
    x_breaks <- seq(0, max(df_wide$time))
    x_labels <- 0:(length(x_breaks)-1)
    x_lab <- "time (days)"
  } else if (max_time <= 365) {
    x_breaks <- seq(0, max(df_wide$time), 365/12)
    x_labels <- 0:(length(x_breaks)-1)
    x_lab <- "time (months)"
  } else {
    x_breaks <- seq(0, max(df_wide$time), 365)
    x_labels <- 0:(length(x_breaks)-1)
    x_lab <- "time (years)"
  }
  
  
  # produce plot
  ggplot2::ggplot(df_long) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes_(x = ~time, y = ~count, color = ~state)) +
    ggplot2::scale_color_manual(values = plot_cols) +
    ggplot2::scale_x_continuous(x_lab, breaks = x_breaks, labels = x_labels)
  
}


#------------------------------------------------
#' @title Plot distribution of treatment seeking in the population
#'
#' @description Within the SIMPLEGEN transmission model each host has its own
#'   treatment seeking parameter, which defines the probability of actively
#'   seeking treatment once disease reaches a point that this becomes necessary.
#'   These parameter values are drawn independently for each host from a Beta
#'   distribution with user-defined mean and standard deviation. This function
#'   plots the Beta distribution for the current project values.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#'
#' @importFrom stats dbeta
#' @export

plot_treatment_seeking <- function(project) {
  
  # check inputs
  assert_custom_class(project, "simplegen_project")
  
  # check epi parameters
  check_epi_params(project$epi_parameters)
  
  # extract treatment seeking parameters
  mu <- project$epi_parameters$treatment_seeking_mean
  sigma <- project$epi_parameters$treatment_seeking_sd
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
    ggplot2::geom_area(ggplot2::aes_(x = ~x, y = ~y), colour = "black", fill = "#4575B499") +
    ggplot2::geom_vline(xintercept = mu, linetype="dashed") +
    ggplot2::xlab("Pr(seek treatment)") + ggplot2::ylab("probability density") +
    ggplot2::ggtitle("Treatment seeking distribution")
  
}

#------------------------------------------------
#' @title Plot age distribution of states
#'
#' @description TODO.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param sample_time which time-orderd sample to plot.
#' @param deme which deme to plot.
#' @param state which host state to plot on the y-axis.
#'
#' @export

plot_age_states <- function(project, sample_time = 1, deme = 1, state = "S") {
  
  # TODO - check inputs
  
  # subset data
  dat <- project$epi_output$age_distributions
  dat <- dat[dat$sample_time == sample_time & dat$deme == deme, c("age", state)]
  names(dat) <- c("age", "y")
  
  # produce plot
  ggplot2::ggplot(dat) + ggplot2::theme_bw() +
    ggplot2::geom_bar(ggplot2::aes_(x = ~age, y = ~y), stat = "identity", width = 1) +
    ggplot2::xlab("age (years)") + ggplot2::ylab("prevalence") +
    ggplot2::ggtitle(sprintf("Age-distribution of state %s", state))
  
}


