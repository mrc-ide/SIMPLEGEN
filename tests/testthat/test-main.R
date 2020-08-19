
test_that("default epi model runs without error", {
  
  # run most basic epi simulation
  p <- simplegen_project() %>%
    define_epi_params() %>%
    sim_epi()
  
  # check output in correct format
  expect_true(assert_custom_class(p, "simplegen_project"))
})
