# calibration_micro_PCR1.R
#
# Author: Bob Verity
# Date: 2023-03-26
#
# Purpose:
# Explores data on parasite prevalence by PCR vs. microscopy, broken down by
# age. This gives information on how much less sensitive microscopy is against
# the gold standard of PCR, and also if there are any trends with prevalence
# and/or age. Compares model predictions by running a sweep of increasing
# mosquito density, each time running for 20 years and taking the average
# results over the final 10 years.
#
# ------------------------------------------------------------------

library(tidyverse)

# define parameter combinations to explore
H <- 1e3
df_model <- data.frame(H = H,
                       m = 10^seq(0, 2, l = 21),
                       Micro_0_5 = NA,
                       Micro_5_15 = NA,
                       Micro_15 = NA,
                       PCR_0_5 = NA,
                       PCR_5_15 = NA,
                       PCR_15 = NA) %>%
  mutate(M = round(m*H)) %>%
  dplyr::select(H, m, M, Micro_0_5, Micro_5_15, Micro_15, PCR_0_5, PCR_5_15, PCR_15)

# create a basic model and load model and sampling parameters
s <- simplegen_project() %>%
  define_epi_model_parameters() %>%
  define_epi_sampling_parameters(daily = rbind(data.frame(name = "Micro_0_5_A", measure = "prevalence", state = "A", diagnostic = "Microscopy", deme = -1, age_min = 0, age_max = 5),
                                               data.frame(name = "Micro_0_5_C", measure = "prevalence", state = "C", diagnostic = "Microscopy", deme = -1, age_min = 0, age_max = 5),
                                               data.frame(name = "Micro_5_15_A", measure = "prevalence", state = "A", diagnostic = "Microscopy", deme = -1, age_min = 5, age_max = 15),
                                               data.frame(name = "Micro_5_15_C", measure = "prevalence", state = "C", diagnostic = "Microscopy", deme = -1, age_min = 5, age_max = 15),
                                               data.frame(name = "Micro_15_A", measure = "prevalence", state = "A", diagnostic = "Microscopy", deme = -1, age_min = 15, age_max = 100),
                                               data.frame(name = "Micro_15_C", measure = "prevalence", state = "C", diagnostic = "Microscopy", deme = -1, age_min = 15, age_max = 100),
                                               data.frame(name = "PCR_0_5_A", measure = "prevalence", state = "A", diagnostic = "PCR", deme = -1, age_min = 0, age_max = 5),
                                               data.frame(name = "PCR_0_5_C", measure = "prevalence", state = "C", diagnostic = "PCR", deme = -1, age_min = 0, age_max = 5),
                                               data.frame(name = "PCR_5_15_A", measure = "prevalence", state = "A", diagnostic = "PCR", deme = -1, age_min = 5, age_max = 15),
                                               data.frame(name = "PCR_5_15_C", measure = "prevalence", state = "C", diagnostic = "PCR", deme = -1, age_min = 5, age_max = 15),
                                               data.frame(name = "PCR_15_A", measure = "prevalence", state = "A", diagnostic = "PCR", deme = -1, age_min = 15, age_max = 100),
                                               data.frame(name = "PCR_15_C", measure = "prevalence", state = "C", diagnostic = "PCR", deme = -1, age_min = 15, age_max = 100)
                                               ))

# run model over parameter combinations
for (i in 1:nrow(df_model)) {
  
  # tweak parameters and run model
  s <- s %>%
    define_epi_model_parameters(H = H,
                                seed_infections = floor(H * 0.1),
                                M = df_model$M[i]) %>%
    sim_epi(max_time = 20*365,
            silent = TRUE)
  
  # get mean over time range
  res <- s$epi_output$daily %>%
    dplyr::filter(time > 10*365) %>%
    colMeans()
  
  # store results
  df_model$Micro_0_5[i] <- res["Micro_0_5_A"] + res["Micro_0_5_C"]
  df_model$Micro_5_15[i] <- res["Micro_5_15_A"] + res["Micro_5_15_C"]
  df_model$Micro_15[i] <- res["Micro_15_A"] + res["Micro_15_C"]
  df_model$PCR_0_5[i] <- res["PCR_0_5_A"] + res["PCR_0_5_C"]
  df_model$PCR_5_15[i] <- res["PCR_5_15_A"] + res["PCR_5_15_C"]
  df_model$PCR_15[i] <- res["PCR_15_A"] + res["PCR_15_C"]
}

# get into long format over age groups
df_long <- rbind.data.frame(data.frame(Age_Group = "0-5", PCR = df_model$PCR_0_5, Micro = df_model$Micro_0_5),
                            data.frame(Age_Group = "5-15years", PCR = df_model$PCR_5_15, Micro = df_model$Micro_5_15),
                            data.frame(Age_Group = "15+", PCR = df_model$PCR_15, Micro = df_model$Micro_15)) %>%
  mutate(Age_Group = factor(Age_Group, levels = c("0-5", "5-15years", "15+")))

# fit linear model to each age group
slope_df <- mapply(function(x) {
  mod <- lm(x$Micro_Prev ~ x$PCR_Prev + 0)
  data.frame(Age_Group = x$Age_Group[1], slope = mod$coefficients[1])
}, split(PCR_micro_age_whittaker2021, f = PCR_micro_age_whittaker2021$Age_Group), SIMPLIFY = FALSE) %>%
  bind_rows()

# plot data and model
PCR_micro_age_whittaker2021  %>%
  mutate(Age_Group = factor(Age_Group, levels = c("0-5", "5-15years", "15+"))) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = PCR_Prev, y = Micro_Prev, color = Age_Group)) +
  geom_abline(aes(intercept = 0, slope = slope, color = Age_Group), linetype = "dashed", data = slope_df) +
  geom_line(aes(x = PCR, y = Micro, color = Age_Group), data = df_long) +
  scale_color_manual(values = simplegen_cols())

# save model output to inst/extdata directory
if (FALSE) {
  saveRDS(df_long, file = "inst/extdata/calibration_micro_PCR1.rds")
}
