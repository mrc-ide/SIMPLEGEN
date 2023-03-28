
#------------------------------------------------
#' @title A series of internally-used colours
#'
#' @description These colours are purely used to make plots look consistent.
#' 
#' @export

simplegen_cols <- function() {
  c("firebrick1", "chartreuse3", "dodgerblue", "dodgerblue4", "purple", "darkorange", "firebrick4")
}

#------------------------------------------------
#' @title Generic function for plotting flexible distributions
#'
#' @description Many distributions within the inbuilt SIMPLEGEN transmission
#'   model can be defined flexibly as either a vector, or a list of vectors.
#'   This function plots these distributions.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param name the name of the distribution to be plotted, which must be one of:
#'   \itemize{
#'     \item \code{"duration_acute"}
#'     \item \code{"duration_chronic"}
#'     \item \code{"time_treatment"}
#'     \item \code{"duration_prophylactic"}
#'     \item \code{"prob_infection"}
#'     \item \code{"prob_acute"}
#'     \item \code{"prob_AC"}
#'     \item \code{"detectability_microscopy_acute"}
#'     \item \code{"detectability_microscopy_chronic"}
#'     \item \code{"detectability_PCR_acute"}
#'     \item \code{"detectability_PCR_chronic"}
#'     \item \code{"infectivity_acute"}
#'     \item \code{"infectivity_chronic"}
#'   }
#'
#' @importFrom stats dbeta
#' @export

plot_flex_distribution <- function(project, name = "duration_acute") {
  
  # define types of allowed distribution
  duration_dist <- c("duration_acute", "duration_chronic", "time_treatment", "duration_prophylactic")
  transition_prob <- c("prob_infection", "prob_acute", "prob_AC")
  daily_prob <- c("detectability_microscopy_acute", "detectability_microscopy_chronic",
                  "detectability_PCR_acute", "detectability_PCR_chronic",
                  "infectivity_acute", "infectivity_chronic")
  
  # check inputs
  assert_class(project, "simplegen_project")
  assert_in(name, c(duration_dist, transition_prob, daily_prob))
  
  # switch between plotting functions for durations vs. transition probabilities
  # vs. daily probabilities
  if (name %in% duration_dist) {
    ret <- plot_flex_distribution_duration(project, name, normalise = TRUE)
  } else if (name %in% transition_prob) {
    ret <- plot_flex_distribution_transition(project, name)
  } else if (name %in% daily_prob) {
    ret <- plot_flex_distribution_duration(project, name, normalise = FALSE)
  }
  
  return(ret)
}

#------------------------------------------------
# plot duration distribution
#' @import ggplot2 dplyr
#' @importFrom grDevices grey
#' @noRd
plot_flex_distribution_duration <- function(project, name, normalise = FALSE) {
  
  # required to remove notes
  x <- NULL
  
  # TODO - check epi model parameters
  #check_epi_model_params(project)
  
  # make plotting dataframe
  y <- project$epi_model_parameters[[name]]
  
  plot_df <- mapply(function(i) {
    ret <- data.frame(x = seq_along(y[[i]]), y = y[[i]], rep = i)
    if (normalise) {
      ret$y <- ret$y / sum(ret$y)
    }
    return(ret)
  }, seq_along(y), SIMPLIFY = FALSE) %>%
    bind_rows()
  
  # produce plot
  ret <- ggplot(data = plot_df) + theme_bw() +
    theme(panel.grid.major = element_blank()) +
    geom_bar(aes(x = x, y = y), color = NA, fill = simplegen_cols()[4],
                      stat = "identity", width = 0.9) +
    scale_y_continuous(limits = c(0, 1.05*max(plot_df$y)), expand = c(0, 0)) +
    xlab("Day") + ylab("Probability") +
    ggtitle(name)
  
  # facet if multiple distributions
  if (length(y) > 1) {
    ret <- ret + facet_wrap(. ~rep) +
      theme(strip.background = element_rect(fill = grey(1.0)))
  }
  
  return(ret)
}

#------------------------------------------------
# plot transition probabilities
#' @import ggplot2
#' @noRd
plot_flex_distribution_transition <- function(project, name) {
  
  # required to remove notes
  x <- NULL
  
  # TODO - check epi model parameters
  #check_epi_model_params(project)
  
  # make plotting dataframe
  y <- project$epi_model_parameters[[name]]
  plot_df <- data.frame(x = seq_along(y), y = y)
  
  # produce plot
  ret <- ggplot(data = plot_df) + theme_bw() +
    geom_bar(aes(x = x, y = y), color = NA, fill = simplegen_cols()[4],
                      stat = "identity", width = 0.9) +
    scale_y_continuous(limits = c(0, 1.05*max(plot_df$y)), expand = c(0, 0)) +
    xlab("exposure") + ylab("probability") +
    ggtitle(name)
  
  return(ret)
}

#------------------------------------------------
#' @title Plot distribution of treatment seeking in the population
#'
#' @description Within the inbuilt SIMPLEGEN transmission model, each host has
#'   its own treatment seeking parameter defining the probability of actively
#'   seeking treatment once disease reaches a point that this becomes necessary.
#'   These parameter values are drawn independently for each host from a Beta
#'   distribution with user-defined mean and standard deviation. This function
#'   plots this Beta distribution for the current project parameter values.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#'
#' @importFrom stats dbeta
#' @importFrom grDevices grey
#' @export

plot_treatment_seeking <- function(project) {
  
  # avoid "no visible binding" errors
  x <- y <- NULL
  
  # check inputs
  assert_class(project, "simplegen_project")
  
  # TODO - check epi model parameters
  # check_epi_model_params(project)
  
  # extract treatment seeking parameters
  mu <- project$epi_model_parameters$treatment_seeking_mean
  sigma <- project$epi_model_parameters$treatment_seeking_sd
  alpha <- mu^2*(1 - mu) / sigma^2 - mu
  beta <- mu*(1 - mu)^2 / sigma^2 - (1 - mu)
  
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
  ggplot(df_plot) + theme_bw() +
    geom_area(aes(x = x, y = y), colour = "black", fill = simplegen_cols()[4]) +
    geom_vline(xintercept = mu, linetype = "dashed") +
    xlab("Pr(seek treatment)") + ylab("Probability density") +
    ggtitle("Treatment seeking distribution")
  
}
