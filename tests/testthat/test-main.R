
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
  expect_equal(names(p$epi_output), c("daily", "sweeps", "surveys"))
  
})

#------------------------------------------------
test_that("pruning generates expected output", {
  
  # make project
  p <- simplegen_project()
  
  # manually define surver details (shortcut simulation)
  p$epi_output$surveys <- data.frame(time = 365,
                                     deme = 1,
                                     host_ID = 1:2,
                                     positive = TRUE,
                                     inoc_IDs = c(15,16))
  
  # locations of transmission record (read from package), and pruned record
  # (write to temporary file)
  transmission_record_location <- system.file("extdata/", "trans_record_test1.txt", package = 'SIMPLEGEN', mustWork = TRUE)
  pruned_record_location <- sprintf("%s/pruned_record_test1.txt", tempdir())
  
  # prune transmission record
  p <- prune_transmission_record(project = p,
                                 transmission_record_location = transmission_record_location,
                                 pruned_record_location = pruned_record_location,
                                 overwrite_pruned_record = TRUE)
  
  # read pruned record, and compare against expected output
  pruned_record <- readLines(pruned_record_location)
  pruned_record_ecpected_location <- system.file("extdata/", "pruned_record_test1.txt", package = 'SIMPLEGEN', mustWork = TRUE)
  pruned_record_ecpected <- readLines(pruned_record_ecpected_location)
  expect_identical(pruned_record, pruned_record_ecpected)
  
})






