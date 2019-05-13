context("test-main")

test_that("dummy1 handles default/null values", {
  expect_equal(dummy1(), dummy1(1:5))
  expect_error(dummy1(NULL))
})
