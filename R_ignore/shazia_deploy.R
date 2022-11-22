#####################################################
##### Title: deploy.R script to run simulation code
##### Author: Shazia Ruybal-Pesántez
##### Date: 2022-11-22
#####################################################

library(SIMPLEGEN)
library(tidyverse)

set.seed(1)

daily_df <- data.frame(name = "prev_acute", deme = 1, state = "A", measure = "prevalence", diagnostic = "true", age_min = 0, age_max = 100) %>% 
                  add_row(name = "prev_chronic", deme = 1, state = "C", measure = "prevalence", diagnostic = "true", age_min = 0, age_max = 100) %>% 
                  add_row(name = "inc_acute", deme = 1, state = "A", measure = "incidence", diagnostic = "true", age_min = 0, age_max = 100) %>% 
                  add_row(name = "inc_chronic", deme = 1, state = "C", measure = "incidence", diagnostic = "true", age_min = 0, age_max = 100) 
  
daily_df

# sweeps_df <- rbind(data.frame(name = "prev_acute", time = 365, deme = 1, measure = "prevalence", state = "A", diagnostic = "true", age_min = seq(0, 100, 5), age_max = seq(0, 100, 5) + 4),
#                          data.frame(name = "inc_acute", time = 365, deme = 1, measure = "incidence",  state = "A", diagnostic = "true", age_min = seq(0, 100, 5), age_max = seq(0, 100, 5) + 4),
#                          make.row.names = FALSE)
# 
# head(sweeps_df)

survey_df <- data.frame(t_start = 2000, t_end = 2000, reporting_interval = 7, deme = 1, measure = "prevalence", sampling = NA, diagnostic = "true", age_min = 0, age_max = 100, sample_size = 100)

survey_df


simproj <- simplegen_project() %>%
  define_epi_model_parameters(H = 1e3,
                              M = 5e4) %>% 
  define_epi_sampling_parameters(daily = daily_df,
                                 # sweeps = sweeps_df, 
                                 surveys = survey_df) %>% 
  sim_epi(max_time = 3650,
          save_transmission_record = TRUE, 
          transmission_record_location = "ignore/transmission_record.csv", 
          overwrite_transmission_record = TRUE, 
          pb_markdown = TRUE)

simproj$epi_sampling_parameters

simproj$epi_output
