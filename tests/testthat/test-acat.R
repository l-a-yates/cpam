test_that("acat works with basic p-value vector", {
  # Basic functionality
  p_values <- c(0.02, 0.0004, 0.2, 0.1, 0.8)
  result <- acat(p_values)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("acat works with equal weights (NULL)", {
  p_values <- c(0.01, 0.05, 0.1)
  result_no_weights <- acat(p_values, weights = NULL)
  result_equal_weights <- acat(p_values, weights = c(1, 1, 1))

  expect_equal(result_no_weights, result_equal_weights)
})

test_that("acat works with custom weights", {
  p_values <- c(0.01, 0.05, 0.1)
  weights <- c(0.5, 0.3, 0.2)
  result <- acat(p_values, weights = weights)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("acat handles very small p-values", {
  p_values <- c(1e-20, 1e-16, 0.5)
  result <- acat(p_values)

  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
})

test_that("acat works with matrix input", {
  p_matrix <- matrix(runif(100), ncol = 10)
  result <- acat(p_matrix)

  expect_type(result, "double")
  expect_length(result, 10)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("acat errors with NA values", {
  p_values <- c(0.01, NA, 0.1)

  expect_error(
    acat(p_values),
    "Cannot have NAs in the p-values!"
  )
})

test_that("acat errors with p-values outside [0,1]", {
  expect_error(
    acat(c(0.01, 1.5, 0.1)),
    "P-values must be between 0 and 1!"
  )

  expect_error(
    acat(c(-0.01, 0.5, 0.1)),
    "P-values must be between 0 and 1!"
  )
})

test_that("acat warns with exact 0 p-values", {
  p_values <- c(0, 0.01, 0.1)

  expect_warning(
    acat(p_values),
    "There are p-values that are exactly 0!"
  )
})

test_that("acat warns with exact 1 p-values", {
  p_values <- c(1, 0.01, 0.1)

  expect_warning(
    acat(p_values),
    "There are p-values that are exactly 1!"
  )
})

test_that("acat errors with both 0 and 1 in same column", {
  p_values <- c(0, 1, 0.5)

  expect_error(
    acat(p_values),
    "Cannot have both 0 and 1 p-values in the same column!"
  )
})

test_that("acat errors when weight dimensions don't match", {
  p_values <- c(0.01, 0.05, 0.1)
  weights <- c(0.5, 0.5)  # Wrong length

  expect_error(
    acat(p_values, weights = weights),
    "The dimensions of weights and Pvals must be the same!"
  )
})

test_that("acat errors with negative weights", {
  p_values <- c(0.01, 0.05, 0.1)
  weights <- c(0.5, -0.3, 0.8)

  expect_error(
    acat(p_values, weights = weights),
    "All the weights must be nonnegative!"
  )
})

test_that("acat errors when all weights are zero", {
  p_values <- c(0.01, 0.05, 0.1)
  weights <- c(0, 0, 0)

  expect_error(
    acat(p_values, weights = weights),
    "At least one weight should be positive in each column!"
  )
})

test_that("acat can skip validity checks", {
  # Should run faster without checks (though we can't test speed easily)
  p_values <- c(0.01, 0.05, 0.1)
  result <- acat(p_values, is.check = FALSE)

  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("acat gives consistent results", {
  set.seed(123)
  p_values <- runif(10)

  result1 <- acat(p_values)
  result2 <- acat(p_values)

  expect_equal(result1, result2)
})

test_that("acat result makes sense for significance", {
  # Very small p-values should give small combined p-value
  small_pvals <- c(0.001, 0.002, 0.003)
  result_small <- acat(small_pvals)

  # Large p-values should give large combined p-value
  large_pvals <- c(0.5, 0.6, 0.7)
  result_large <- acat(large_pvals)

  expect_true(result_small < result_large)
})
