test_that("lancaster works with basic p-value vector", {
  result <- lancaster(c(0.01, 0.05, 0.1), c(1, 1, 1))
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("lancaster handles single p-value", {
  result <- lancaster(0.05, 1)
  expect_equal(result, 0.05)
})

test_that("lancaster returns NA for all-NA input", {
  result <- lancaster(c(NA), c(1))
  expect_true(is.na(result))
})

test_that("lancaster errors when lengths differ", {
  expect_error(lancaster(c(0.01, 0.05), c(1)), "Length of weights not equal")
})

test_that("lancaster handles NA p-values by removing them", {
  result <- lancaster(c(0.01, NA, 0.1), c(1, 1, 1))
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
})

test_that("lancaster warns about extreme p-values", {
  expect_warning(
    lancaster(c(1e-321, 0.05), c(1, 1)),
    "Extreme p-values"
  )
})

test_that("lancaster gives consistent results", {
  result1 <- lancaster(c(0.01, 0.05, 0.1), c(1, 2, 3))
  result2 <- lancaster(c(0.01, 0.05, 0.1), c(1, 2, 3))
  expect_equal(result1, result2)
})

test_that("lancaster result makes sense for significance", {
  small_pvals <- c(0.001, 0.002, 0.003)
  result_small <- lancaster(small_pvals, c(1, 1, 1))

  large_pvals <- c(0.5, 0.6, 0.7)
  result_large <- lancaster(large_pvals, c(1, 1, 1))

  expect_true(result_small < result_large)
})

test_that("lancaster handles zero-weight entries", {
  result <- lancaster(c(0.01, 0.05, 0.1), c(1, 0, 1))
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
})
