
#------------------------------------------------
test_that("default epi model fails with bad input", {
  
  # create basic project
  p <- simplegen_project() %>%
    define_epi_params()
  
  # run model with bad input
  expect_error(sim_epi(p, max_time = 10, output_age_times = c(0,10)))
  expect_error(sim_epi(p, max_time = 10, output_age_times = c(1,11)))
  
})

#------------------------------------------------
test_that("default epi model runs and produces correct output", {
  
  # define expected project slot names
  expected_names <- c("epi_parameters", "sampling_strategy", "epi_output", "sample_details", 
                      "relatedness", "true_genotypes", "observed_genotypes")
  
  # create basic project
  p <- simplegen_project()
  expect_equal(names(p), expected_names)
  
  # load default epi parameters
  p <- define_epi_params(p)
  expect_equal(names(p), expected_names)
  
  # run model with minimal output
  p <- sim_epi(p, output_daily_counts = FALSE, output_age_distributions = FALSE)
  
  # check that expected elements are present in output lists
  expect_equal(names(p), expected_names)
  expect_equal(names(p$epi_output), c("daily_values", "age_distributions"))
  
  # check that return objects are null if they were not requested to be output
  expect_null(p$epi_output$daily_values)
  expect_null(p$epi_output$age_distributions)
})
