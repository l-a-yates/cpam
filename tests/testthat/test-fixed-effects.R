# tests/testthat/test-validate-fixed-effects

test_that("validate_fixed_effects handles NULL input", {
  exp_design <- data.frame(
    time = 1:10,
    treatment = rep(c("A", "B"), each = 5),
    block = rep(1:5, 2)
  )

  expect_null(validate_fixed_effects(NULL, exp_design))
})

test_that("validate_fixed_effects validates formula class", {
  exp_design <- data.frame(
    time = 1:10,
    treatment = rep(c("A", "B"), each = 5),
    block = rep(1:5, 2)
  )

  expect_error(
    validate_fixed_effects("treatment + block", exp_design),
    "fixed_effects must be a formula"
  )
})

test_that("validate_fixed_effects handles empty formula", {
  exp_design <- data.frame(
    time = 1:10,
    treatment = rep(c("A", "B"), each = 5),
    block = rep(1:5, 2)
  )

  expect_null(validate_fixed_effects(~ 1, exp_design))
})

test_that("validate_fixed_effects checks for missing terms", {
  exp_design <- data.frame(
    time = 1:10,
    treatment = rep(c("A", "B"), each = 5),
    block = rep(1:5, 2)
  )

  expect_error(
    validate_fixed_effects(~ treatment + missing_var, exp_design),
    "The following fixed effects variables are missing from exp_design: missing_var"
  )

  expect_error(
    validate_fixed_effects(~ missing_var1 + missing_var2, exp_design),
    "The following fixed effects variables are missing from exp_design: missing_var1, missing_var2"
  )
})

test_that("validate_fixed_effects returns correct string for valid input", {
  exp_design <- data.frame(
    time = 1:10,
    treatment = rep(c("A", "B"), each = 5),
    block = rep(1:5, 2)
  )

  expect_equal(
    validate_fixed_effects(~ treatment, exp_design),
    "treatment"
  )

  expect_equal(
    validate_fixed_effects(~ block, exp_design),
    "block"
  )
})

test_that("validate_fixed_effects warns about time collinearity", {
  # Create design where treatment perfectly correlates with time
  exp_design <- data.frame(
    time = 1:10,
    treatment = rep(c("A", "B"), 5),
    block = rep(1:5, 2)
  )

  # Force collinearity by making treatment fully determined by time
  exp_design$collinear_var <- exp_design$time

  expect_error(
    validate_fixed_effects(~ collinear_var, exp_design),
    "The model matrix is not full rank. This is likely due to collinear fixed effects."
  )
})


test_that("check_perfect_separation warns about singleton levels", {
  exp_design <- data.frame(
    time = 1:10,
    treatment = c("A", "A", "A", "A", "A", "B", "B", "B", "B", "C"),
    block = rep(1:5, 2)
  )

  expect_warning(
    check_perfect_separation(exp_design, c("treatment")),
    "Perfect separation detected in the following terms: treatment"
  )

  # No warning when all levels have multiple observations
  exp_design$balanced <- rep(c("X", "Y"), 5)
  expect_no_warning(check_perfect_separation(exp_design, c("balanced")))
})

test_that("check_perfect_separation handles multiple categorical variables", {
  exp_design <- data.frame(
    time = 1:10,
    cat1 = c("A", "A", "A", "A", "A", "B", "B", "B", "B", "C"),
    cat2 = c(rep("X", 9), "Y"),
    num = 1:10
  )

  # Should warn about both problematic variables
  expect_warning(
    check_perfect_separation(exp_design, c("cat1", "cat2", "num")),
    "Perfect separation detected in the following terms: cat1, cat2"
  )

})

