
#------------------------------------------------
test_that("default epi model fails with bad input", {
  
  # create empty project
  p <- simplegen_project()
  
  # should fail with no epi model parameters loaded
  expect_error(sim_epi(p))
  
  # load epi model parameters
  p <- define_epi_model_parameters(p)
  
  # should fail with no epi sampling parameters loaded
  expect_error(sim_epi(p))
  
  # define empty sampling parameters
  p <- define_epi_sampling_parameters(p)
  
  # should now run without errors or warnings
  expect_error(sim_epi(p), regexp = NA)
  expect_warning(sim_epi(p), regexp = NA)
  
})

#------------------------------------------------
test_that("default epi model runs and produces correct output", {
  
  # define expected project slot names
  expected_names <- c("epi_model_parameters", "epi_sampling_parameters", "epi_output", 
                      "genetic_parameters", "relatedness", "true_genotypes", "observed_genotypes")
  
  # create basic project
  p <- simplegen_project()
  expect_equal(names(p), expected_names)
  
  # load default epi parameters
  p <- define_epi_model_parameters(p)
  expect_equal(names(p), expected_names)
  
  # define NULL sampling parameters
  p <- define_epi_sampling_parameters(p)
  
  # run model
  p <- sim_epi(p)
  
  # check that expected elements are present in output lists
  expect_equal(names(p), expected_names)
  expect_equal(names(p$epi_output), c("daily", "sweeps", "surveys"))
  
  # check that return objects are null if they were not requested to be output
  expect_null(p$epi_output$daily)
  expect_null(p$epi_output$sweeps)
  expect_null(p$epi_output$surveys)
  
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






