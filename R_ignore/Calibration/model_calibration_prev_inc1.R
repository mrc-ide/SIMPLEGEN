# model_calibration_prev_inc1.R
#
# Author: Bob Verity
# Date: 2023-03-25
#
# Purpose:
# Compares model output against EIR vs. prevalence data. Runs a sweep of
# increasing mosquito density, each time running for 20 years and taking the
# average results over the final 10 years. Prevalence is calculated as the sum
# of acute and chronic microscopy prevalence in children (<15) to best match the
# description of the data.
#
# ------------------------------------------------------------------

# define parameter combinations to explore
H <- 1e3
df_model <- data.frame(H = H,
                       m = 10^seq(0, 2, l = 51),
                       annual_EIR = NA,
                       prev = NA) %>%
  mutate(M = round(m*H)) %>%
  dplyr::select(H, m, M, annual_EIR, prev)

# create project, load parameters and sampling design
s <- simplegen_project() %>%
  define_epi_model_parameters(H = H,
                              seed_infections = floor(H * 0.1)) %>%
  define_epi_sampling_parameters(daily = rbind(data.frame(name = "EIR", measure = "EIR", state = NA, diagnostic = NA, deme = -1, age_min = 0, age_max = 100),
                                               data.frame(name = "prev_A", measure = "prevalence", state = "A", diagnostic = "Microscopy", deme = -1, age_min = 0, age_max = 15),
                                               data.frame(name = "prev_C", measure = "prevalence", state = "C", diagnostic = "Microscopy", deme = -1, age_min = 0, age_max = 15)))

# run model over parameter combinations
for (i in 1:nrow(df_model)) {
  
  # finalise parameters and run model
  s <- s %>%
    define_epi_model_parameters(M = df_model$M[i]) %>%
    sim_epi(max_time = 20*365,
            silent = TRUE)
  
  # get mean over time range
  res <- s$epi_output$daily %>%
    dplyr::filter(time > 10*365) %>%
    mutate(prev = prev_A + prev_C) %>%
    dplyr::select(EIR, prev) %>%
    colMeans()
  
  # store results
  df_model$annual_EIR[i] <- 365 * res["EIR"]
  df_model$prev[i] <- res["prev"]
}

# filter model output to positive EIR
df_model <- df_model %>%
  dplyr::filter(annual_EIR > 0)

# plot against Beier data
EIR_prev_beier1999 %>%
  dplyr::filter(annual_EIR > 0) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = annual_EIR, y = prevalence), col = simplegen_cols()[3]) +
  scale_x_log10(limits = c(1e-2, 1e3)) +
  ylim(c(0, 1)) + xlab("Annual EIR") + ylab("Prevalence") +
  geom_line(aes(x = annual_EIR, y = prev), size = 0.8, col = simplegen_cols()[1], data = df_model)

# plot against Hay data
EIR_prev_hay2005 %>%
  dplyr::filter(annual_EIR > 0) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = annual_EIR, y = prevalence), col = simplegen_cols()[3]) +
  scale_x_log10(limits = c(1e-2, 1e3)) +
  ylim(c(0, 1)) + xlab("Annual EIR") + ylab("Prevalence") +
  geom_line(aes(x = annual_EIR, y = prev), size = 0.8, col = simplegen_cols()[1], data = df_model)

# save model output to inst/extdata directory
if (FALSE) {
  saveRDS(df_model, file = "inst/extdata/model_calibration_prev_inc1.rds")
}
