
#------------------------------------------------
test_that("epi model fails to run if input is incomplete", {
  
  # create empty project
  p <- simplegen_project()
  
  # should fail with no model or sampling parameters loaded
  expect_error(sim_epi(p))
  
  # should fail with only model parameters loaded
  p <- define_epi_model_parameters(p)
  expect_error(sim_epi(p))
  
  # should fail with only sampling parameters loaded
  daily_df <- data.frame(deme = 1, state = "S", measure = "prevalence", diagnostic = NA, age_min = 0, age_max = 100)
  p <- simplegen_project() %>%
    define_epi_sampling_parameters(daily = daily_df)
  expect_error(sim_epi(p)) 
  
})

#------------------------------------------------
test_that("default epi model runs and maintains expected project structure", {
  
  # create basic project
  p <- simplegen_project()
  
  # check for expected project slot names throughout
  expected_names <- c("epi_model_parameters", "epi_sampling_parameters", "epi_output", 
                      "genetic_parameters", "relatedness", "true_genotypes", "observed_genotypes")
  expect_equal(names(p), expected_names)
  
  # load default epi parameters
  p <- define_epi_model_parameters(p)
  expect_equal(names(p), expected_names)
  
  # load a simple daily output dataframe
  daily_df <- data.frame(deme = 1, state = "S", measure = "prevalence", diagnostic = NA, age_min = 0, age_max = 100)
  p <- define_epi_sampling_parameters(project = p, daily = daily_df)
  expect_equal(names(p), expected_names)
  
  # should now run without errors or warnings
  expect_error(sim_epi(p), regexp = NA)
  expect_warning(sim_epi(p), regexp = NA)
  
  # run model
  p <- sim_epi(p)
  expect_equal(names(p), expected_names)
  
  # check that expected elements are present in output
  expect_equal(names(p$epi_output), c("daily", "sweeps", "surveys_indlevel", "surveys_summary"))
  
})






