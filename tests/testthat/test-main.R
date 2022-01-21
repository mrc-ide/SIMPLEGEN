
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
                      "sample_details", "pruned_record", "genetic_parameters", "haplotype_tree", "block_tree",
                      "true_genotypes", "observed_genotypes")
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
test_that("block coalescent algorithm produces expected results", {
  
  # create basic project
  p <- simplegen_project() %>%
    define_genetic_params()
  p$genetic_parameters$contig_lengths <- NULL
  
  # test1: error if no contig lengths defined within genetic parameters
  expect_error(get_haplotype_coalescence(p, haplo_IDs = c(100, 99)))
  
  # define contig lengths
  p <- define_genetic_params(p, contig_lengths = 100)
  
  # test2: error if no block tree loaded
  expect_error(get_haplotype_coalescence(p, haplo_IDs = c(100, 99)))
  
  # test3: error if target IDs not present in tree
  p$block_tree <- simplegen_file("block_coalescence_1.csv")
  expect_error(get_haplotype_coalescence(p, haplo_IDs = c(1, 2)))
  
  # test4: coalesce directly to de-novo generation. Expect no coalescent output
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test5: contigs in tree go beyond those specified in genetic parameters
  p$block_tree <- simplegen_file("block_coalescence_2.csv")
  expect_error(get_haplotype_coalescence(p, haplo_IDs = c(100, 99)))
  
  # test5: coalesce to different contigs. Expect no coalescent output
  p$genetic_parameters$contig_lengths <- rep(100, 2)
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test6: stimulate get_overlap() config1. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_3.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test7: stimulate get_overlap() config2. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_4.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test8: stimulate get_overlap() config3. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_5.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test9: stimulate get_overlap() config4. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_6.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test10: stimulate get_overlap() config5. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_7.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test11: stimulate get_overlap() config6. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_8.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test12: completely distinct lineages, no opportunity for coalescence
  p$block_tree <- simplegen_file("block_coalescence_9.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test13: check_coalescence, config1. Expect no coalescent output
  p$block_tree <- simplegen_file("block_coalescence_10.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test14: check_coalescence, config2.
  p$block_tree <- simplegen_file("block_coalescence_11.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  z_expected <- data.frame(lineage1 = 99,
                           lineage2 = 100,
                           contig = 1,
                           left = 51,
                           right = 100,
                           time = 99)
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
  # test15: check_coalescence, config3.
  p$block_tree <- simplegen_file("block_coalescence_12.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  z_expected <- data.frame(lineage1 = 99,
                           lineage2 = 100,
                           contig = 1,
                           left = 51,
                           right = 60,
                           time = 99)
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
  # test16: check_coalescence, config4.
  p$block_tree <- simplegen_file("block_coalescence_13.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  z_expected <- data.frame(lineage1 = 99,
                           lineage2 = 100,
                           contig = 1,
                           left = 51,
                           right = 60,
                           time = 99)
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
  # test17: check_coalescence, config5.
  p$block_tree <- simplegen_file("block_coalescence_14.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  z_expected <- data.frame(lineage1 = 99,
                           lineage2 = 100,
                           contig = 1,
                           left = 1,
                           right = 50,
                           time = 99)
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
  # test18: check_coalescence, config6. Expect no coalescent output 
  p$block_tree <- simplegen_file("block_coalescence_15.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99))
  expect_equal(dim(z), c(0, 0))
  
  # test19: series of single-locus coalescences
  p$block_tree <- simplegen_file("block_coalescence_16.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99), generations = 99)
  z_expected <- data.frame(lineage1 = 99,
                           lineage2 = 100,
                           contig = 1,
                           left = 1:5,
                           right = 1:5,
                           time = 6:2)
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
  # test20: repeated back-crosses over 5 generations. Haplos are related in
  # blocks as follows:
  # 16 loci in generation 1 at time 98
  #  8 loci in generation 2 at time 95
  #  4 loci in generation 3 at time 90
  #  2 loci in generation 4 at time 86
  # 1 locus in generation 5 at time 82
  p$block_tree <- simplegen_file("block_coalescence_17.csv")
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99), generations = 99)
  z_expected <- data.frame(lineage1 = 99,
                           lineage2 = 100,
                           contig = 1,
                           left = c(1, 17, 25, 29, 31),
                           right = c(16, 24, 28, 30, 31),
                           time = c(98, 94, 90, 86, 82))
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
  # test21: as before, but check that generation limiter working
  z <- get_haplotype_coalescence(p, haplo_IDs = c(100, 99), generations = 2)
  expect_true(dplyr::all_equal(z, z_expected[1:2,], convert = TRUE))
  
  # test22: simulated example that was used to diagnose bugs at one point.
  # Solution was worked through manually in excel
  p$block_tree <- simplegen_file("block_coalescence_18.csv")
  p$genetic_parameters$contig_lengths <- 3291871
  z <- get_haplotype_coalescence(p, haplo_IDs = c(228, 230), generations = 99)
  z_expected <- simplegen_file("block_coalescence_18_solution.csv")
  expect_true(dplyr::all_equal(z, z_expected, convert = TRUE))
  
})




