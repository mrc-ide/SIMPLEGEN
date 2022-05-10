
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
#'   model can be defined flexibly as either a vector or a list of vectors. This
#'   function produces simple visualisations of these distributions.
#'
#' @param project a SIMPLEGEN project, as produced by the
#'   \code{simplegen_project()} function.
#' @param name the name of the distribution to be plotted, which must be one of:
#'   \itemize{
#'     \item \code{"duration_acute"}
#'     \item \code{"duration_chronic"}
#'     \item \code{"time_treatment_acute"}
#'     \item \code{"time_treatment_chronic"}
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
#' Plotting SIMPLEGEN outputs against epidemiological data
#'
#' @param model_out list of SIMPLEGEN model outputs from running sim_epi function
#' @param dataname  Data to compare SIMPLEGEN output to. "micro_vs_PCR" for plotting
#'  the relationship between prevalence by microscopy and prevalence by PCR.
#'  "prevalence_EIR" is the relationship between Annual EIR and prevalence,
#'  "prevalence_incidence" plots the relationship between Prevalence in 2-10 year olds
#'   against incidence in 0-5 year olds, and the modelled relationship between these 
#'   from Griffin et a. 2014
#'   
#'   TODO - Need to add Griffin et al. 2014  model fit data and example list of SIMPLEGEN
#'   model outputs 
#'
#'
#' @export 
#' @import tidyverse, ggpubr, ggplot2, 
#'
#' @examples
#'
#' plot_epi_data(model_out, "micro_vs_PCR")
#' plot_epi_data(model_out, "prevalence_EIR")
#' plot_epi_data(model_out, "prevalence_incidence")

plot_epi_data <- function(model_out, dataname) {
  
  #check input
  assert_class(model_out[[1]], "simplegen_project")
  assert_class(model_out, "list")
  assert_in(dataname, c("micro_vs_PCR", "prevalence_EIR", "prevalence_incidence"))
  
  
  # Microscopy detection plot
  if (dataname == "micro_vs_PCR") {
    out <- lapply(
      FUN = function(x) {
        prev <- model_out[[x]]$epi_output$surveys %>%
          filter(measure == "prevalence" &
                   diagnostic == "microscopy") %>%
          select(value) %>%
          rename(prev_micro = value)
        
        incid <- model_out[[x]]$epi_output$surveys %>%
          filter(measure == "prevalence" & diagnostic == "PCR") %>%
          select(value) %>%
          rename(prev_PCR = value)
        
        res <- cbind(prev, incid)
        return(res)
      },
      X = 1:length(model_out)
    )
    
    
    data("PCR_micro_full_whittaker2021")
    SIMPLEGEN_dat <- do.call(rbind.data.frame, out)
    
    
    # error if nothing to plot
    if (nrow(SIMPLEGEN_dat) == 0) {
      stop("SIMPLEGEN output doesn't exist or incompatible with function")
    }
    
    # error if nothing to plot
    if (nrow(PCR_micro_full_whittaker2021) == 0) {
      stop("data missing or formatted incorrectly")
    }
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = PCR_micro_full_whittaker2021,
        aes(x = PCR_Prev, y = Micro_Prev, col = "Whittaker et al. 2021"),
        size = 2
      ) +
      ggplot2:: geom_point(
        data = SIMPLEGEN_dat,
        aes(
          x = prev_PCR / 100,
          y = prev_micro / 100,
          col = "SIMPLEGEN"
        ),
        size = 2
      ) +
      ggplot2::geom_abline(intercept = 0, lwd = 1) +
      ggplot2::theme_bw() +
      ggpubr::labs_pubr(base_size = 14) +
      ggplot2::xlab("Prevalence by PCR") +
      ggplot2::ylab("Prevalence by microscopy") +
      ggplot2::labs(colour = "Source")
    p
    set_palette(p, "npg")
    
    return(p)
  }
  
  
  # Microscopy detection plot
  if (dataname == "prevalence_EIR") {
    # check ages for reported prevalence
    
    data("EIR_prev_hay2005")
    #data("EIR_prev_beier1999")
    
    out <- lapply(function(x) {
      # get prevalence
      prev <- model_out[[x]]$epi_output$surveys %>%
        filter(measure == 'prevalence' &
                 diagnostic == "microscopy") %>%
        select(value) %>%
        rename(prev = value)
      
      
      # get EIR
      EIR <- model_out[[x]]$epi_output$daily %>%
        filter(
          measure == "EIR" &
            time == unique(model_out[[x]]$epi_output$surveys$reporting_time)
        ) %>%
        select(value) %>%
        rename(EIR = value)
      
      res <- data.frame(EIR = EIR,
                        prev = mean(prev$prev, na.rm = TRUE))
      return(res)
      
    }, X = 1:length(model_out))
    
    SIMPLEGEN_dat <- do.call(rbind.data.frame, out)
    
    # error if nothing to plot
    if (nrow(SIMPLEGEN_dat) == 0) {
      stop("SIMPLEGEN output doesn't exist or incompatible with function")
    }
    
    # error if nothing to plot
    if (nrow(EIR_prev_hay2005) == 0) {
      stop("data missing or formatted incorrectly")
    }
    
    
    p <- ggplot() +
      geom_point(
        data = EIR_prev_hay2005,
        aes(x = annual_EIR, y = prevalence, col = "Hay et al. 2005"),
        size = 2
      ) +
      geom_point(data =  SIMPLEGEN_dat,
                 aes(x = EIR, y = prev / 100, col = "SIMPLEGEN"),
                 size = 2) +
      theme_bw() +
      labs_pubr(base_size = 14) +
      xlim(c(0, 750)) +
      xlab("Annual EIR") +
      ylab("Prevalence by Microscopy") +
      labs(colour = "Source")
    p
    set_palette(p, "npg")
    return(p)
  }
  
  
  # Prev-incidence relationship by age
  if (dataname == "prevalence_incidence") {
    # check ages for reported prevalence
    
    data("prev_inc_griffin2014")
    
    # Griffin Model data must be added to package
    # data("griffinmod")
    # griffin_mod$site_name <- "griffin_model"
    # griffin_mod <- griffin_mod[, c(3, 1, 2)]
    # colnames(griffin_mod) <- c("site_name", "Prev_2_10", "incid_0_5")
    # griffin_mod$Prev_2_10 <- as.numeric(griffin_mod$Prev_2_10) / 100
    
    # adding data used to fit griffin et al
    prev_inc_griffin2014$value <-
      prev_inc_griffin2014$numer / prev_inc_griffin2014$denom
    
    # separate out prevalence and incidence
    prev <- prev_inc_griffin2014 %>% filter(type == "prevalence")
    incid <- prev_inc_griffin2014 %>% filter(type == "incidence")
    
    # for prevalence, active or passive detection
    passive_det <- prev %>% filter(case_detection == "PCD")
    active_det <- prev %>% filter(case_detection != "PCD")
    
    # filter to ages 2-10
    active_prev <- prev_inc_griffin2014 %>%
      filter(site_index %in% unique(active_det$site_index)) %>%
      filter(age1 < 11) %>%
      filter(age0 > 1) %>%
      filter(type == 'prevalence')
    
    
    # filter to age 0 - 5
    incid <- prev_inc_griffin2014 %>%
      filter(age1 < 6) %>%
      filter(type == 'incidence')
    
    dat_reshaped <- as.data.frame(matrix(ncol = 3, nrow = 3))
    colnames(dat_reshaped) <-
      c("site_name", "Prev_0_10", "incid_0_5")
    
    for (i in 1:length(unique(active_prev$site_name))) {
      site_n <- unique(active_prev$site_name)[i]
      i_tmp <- incid %>% filter(site_name == site_n)
      ap_tmp <- active_prev %>% filter(site_name  == site_n)
      
      dat_reshaped[i, "site_name"] <- site_n
      dat_reshaped[i, "Prev_0_10"] <-
        sum(ap_tmp$numer) / sum(ap_tmp$denom)
      dat_reshaped[i, "incid_0_5"] <-
        sum(i_tmp$numer) / sum(i_tmp$denom)
      
    }
    
    
    
    out <- lapply(function(x) {
      prev210 <- model_out[[x]]$epi_output$surveys %>%
        filter(measure == 'prevalence' &
                 age_max == 10 & diagnostic == "PCR") %>%
        select(value) %>%
        rename(prev210 = value)
      
      
      incid05 <- model_out[[x]]$epi_output$surveys %>%
        filter(measure == 'incidence' & age_max == 5) %>%
        select(value) %>%
        rename(incid05 = value)
      
      
      # incid05 <- model_out[[x]]$epi_output$daily %>%
      #   filter(
      #     measure == 'incidence' &
      #       time <  model_out[[x]]$epi_output$surveys$reporting_time &
      #       time >  (
      #         model_out[[x]]$epi_output$surveys$reporting_time[1] - model_out[[x]]$epi_output$surveys$reporting_interval[1]
      #       )
      #   ) %>%
      #   #select(value)%>%
      #   summarize(incid05 = sum(value, na.rm = TRUE)/sum(denom%>%
      #   as.data.frame()
      # colnames(incid05)<- c("incid05")
      
      
      res <- cbind(prev210, incid05)
      return(res)
      
    }, X = 1:length(model_out))
    
    # combine model output and griffin et al 2014 data
    SIMPLEGEN_dat <- do.call(rbind.data.frame, out)
    
    
    SIMPLEGEN_dat$site_name <-
      rep('SIMPLEGEN', length.out = nrow(SIMPLEGEN_dat))
    dat_reshaped <-
      dat_reshaped %>% mutate(site_name = gsub('_.*', '', site_name))
    
    # error if nothing to plot
    if (nrow(SIMPLEGEN_dat) == 0) {
      stop("SIMPLEGEN output doesn't exist or incompatible with function")
    }
    
    # error if nothing to plot
    if (nrow(dat_reshaped) == 0) {
      stop("data missing or formatted incorrectly")
    }
    p <-
      ggplot2::ggplot(dat_reshaped,
                      ggplot2::aes(x =  Prev_0_10, y = incid_0_5, col = site_name)) +
      ggplot2::theme_bw() +
      # ggplot2::theme(legend.position = "none") +
      ggplot2::geom_point(size = 2) +
      labs_pubr(base_size = 14) +
      ggplot2::xlab("Prevalence in 2 - 10 year olds") +
      ggplot2::ylab("Incidence in 0 - 5 year olds") +
      ggplot2::ggtitle("Prevalence in 2-10 year olds plotted against incidence in 0-5 year olds")
    # p <-
    #   p + geom_line(data = griffin_mod,
    #                 aes(x = Prev_2_10, y = incid_0_5),
    #                 size = 1)
    p <- p + geom_point(data = SIMPLEGEN_dat,
                        aes(x = prev210 / 100, y = incid05),
                        size = 2)
    
    set_palette(p, "npg")
    return(p)
    
  }
  
}




