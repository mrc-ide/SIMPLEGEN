
#------------------------------------------------
test_that("sampling strategy fails with bad input", {
  
  # create project and load epi parameters
  p <- simplegen_project() %>%
    define_epi_model_parameters()
  
  # define bad strategies
  df_sample <- data.frame(time = c(100, 10),
                          deme = 1,
                          case_detection = "passive",
                          diagnosis = "microscopy",
                          n = 5)
  expect_error(define_sampling_strategy(p, df_sample))
  
})

#------------------------------------------------
test_that("default epi model fails with bad input", {
  
  # create empty project
  p <- simplegen_project()
  
  # functions that should fail with no epi parameters loaded
  expect_error(sim_epi(p))
  
  # load epi parameters
  p <- define_epi_model_parameters(p)
  
  # run model with bad input
  expect_error(sim_epi(p, max_time = 100, output_age_times = c(0,100)))
  expect_error(sim_epi(p, max_time = 100, output_age_times = c(1,101)))
  
  # define sampling strategy
  df_sample <- data.frame(time = 100,
                          deme = 1,
                          case_detection = "passive",
                          diagnosis = "microscopy",
                          n = 5)
  p <- define_sampling_strategy(p, df_sample)
  
  # run model with input that clashes with sampling strategy
  expect_error(sim_epi(p, max_time = 10))
  p$sampling_strategy$deme <- 2
  expect_error(sim_epi(p, max_time = 100))
  
})

#------------------------------------------------
test_that("default epi model runs and produces correct output", {
  
  # define expected project slot names
  expected_names <- c("epi_parameters", "sampling_strategy", "epi_output", "sample_output", 
                      "genetic_parameters", "relatedness", "true_genotypes", "observed_genotypes")
  
  # create basic project
  p <- simplegen_project()
  expect_equal(names(p), expected_names)
  
  # load default epi parameters
  p <- define_epi_model_parameters(p)
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

#------------------------------------------------
test_that("pruning generates expected output", {
  
  # make project
  p <- simplegen_project()
  
  # manually define sample details (shortcut simulation)
  p$sample_output <- data.frame(time = 365,
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






