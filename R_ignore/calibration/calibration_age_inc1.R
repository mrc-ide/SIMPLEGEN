# calibration_age_inc1.R
#
# Author: Bob Verity
# Date: 2023-03-26
#
# Purpose:
# Explores patterns of age vs. clinical incidence over multiple sites. Attempts
# to find a mosquito density (and implied EIR) for each site that reproduces the
# observed pattern with age. This is very difficult under the current model
# structure, partly due to some missing element, but also potentially due to
# site-specific differences that need to be taken into account.
#
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# create sampling data.frame representing incidence in 1-year age bands. Also
# store raw population counts in these same age bands to provide denominator for
# incidence calculation. Finally, store EIR
daily_df <- rbind(data.frame(name = "inc",
                              measure = "incidence_active",
                              state = "A",
                              diagnostic = "Microscopy",
                              deme = -1,
                              age_min = 0:49),
                  data.frame(name = "pop",
                             measure = "count",
                             state = "H",
                             diagnostic = NA,
                             deme = -1,
                             age_min = 0:49)) %>%
  mutate(age_max = age_min + 1,
         name = sprintf("%s_%s_%s", name, age_min, age_max)) %>%
  rbind(data.frame(name = "EIR", measure = "EIR", state = NA, diagnostic = NA, deme = -1, age_min = 0, age_max = 100))

# create project and load model and sampling parameters
s <- simplegen_project() %>%
  define_epi_model_parameters() %>%
  define_epi_sampling_parameters(daily = daily_df)

# manually define mosquito population density per site
df_m <- data.frame(m = c(1.1, 1.1, 1.2, 1.3,
                         1.3, 1.7, 1.7, 2, 
                         2, 2.5, 2.4, 3.1,
                         3, 3, 12, 4,
                         12),
                   EIR = NA)

# run model over all sites
ret_list <- list()
for (i in 1:nrow(df_m)) {
  
  # run model
  H <- 1e4
  s <- s %>%
    define_epi_model_parameters(H = H,
                                seed_infections = floor(H * 0.1),
                                M = round(H * df_m$m[i])) %>%
    sim_epi(max_time = 50*365)
  
  # get EIR
  df_m$EIR[i] <- s$epi_output$daily %>%
    dplyr::filter(time > 30*365) %>%
    pull(EIR) %>%
    mean() %>%
    "*"(365)
  
  # take mean over time, get into long format
  df_long <- s$epi_output$daily %>%
    dplyr::filter(time > 30*365) %>%
    pivot_longer(-c(time, EIR), names_to = "group", values_to = "value") %>%
    group_by(group) %>%
    summarise(value = mean(value)) %>%
    separate(group, into = c("type", "age_min", "age_max"), sep = "_")
  
  # divide by denominator population
  df_summary <- df_long %>%
    pivot_wider(names_from = type, values_from = value) %>%
    mutate(age_min = as.numeric(age_min),
           age_max = as.numeric(age_max),
           inc = inc / pop) %>%
    dplyr::select(age_min, age_max, inc)
  
  # store results for this site (yearly incidence)
  ret_list[[i]] <- data.frame(site = i,
                              age_min = df_summary$age_min,
                              age_max = df_summary$age_max,
                              inc = 365* df_summary$inc)
}
df_sim <- ret_list %>%
  bind_rows()

# subset data to ACD only and get into plotting format
df_data <- prev_inc_griffin2014 %>%
  dplyr::filter(type == "incidence") %>%
  dplyr::filter(case_detection != "PCD") %>%
  mutate(inc = numer / denom,
         lower = CI_inc(numer, denom)$lower,
         upper = CI_inc(numer, denom)$upper,
         age_mid = (age0 + age1) / 2) %>%
  dplyr::filter(age_mid < 20)

# get order of increasing incidence in data
df_order <- df_data %>%
  group_by(plot_title) %>%
  summarise(max_inc = max(inc)) %>%
  arrange(max_inc)

# add EIR to plot title
df_order <- df_order %>%
  mutate(EIR = df_m$EIR,
         plot_title_EIR = sprintf("%s\nEIR=%s", plot_title, signif(EIR, 2)),
         plot_title_EIR = factor(plot_title_EIR, levels = plot_title_EIR))

# assume sim results are in same order
df_sim <- df_sim %>%
  mutate(plot_title_EIR = df_order$plot_title_EIR[site])

# plot data
df_data %>%
  left_join(df_order) %>%
  ggplot() + theme_bw() +
  geom_hline(yintercept = 0, alpha = 0) +
  geom_errorbar(aes(x = age_mid, ymin = lower, ymax = upper, color = case_detection)) +
  geom_point(aes(x = age_mid, y = inc, color = case_detection)) +
  xlim(c(0, 20)) + xlab("Age, years") + ylab("Incidence per year") +
  facet_wrap(~plot_title_EIR, ncol = 4, scales = "free_y") +
  geom_line(aes(x = age_min, y = inc), size = 0.5, data = df_sim)

# save model output to inst/extdata directory
if (FALSE) {
  saveRDS(list(df_sim = df_sim,
               df_order = df_order),
          file = "inst/extdata/calibration_age_inc1.rds")
}
