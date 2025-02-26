test_that("results() returns expected basic output", {
  # Basic call without filters
  res <- results(cpo_example)

  # Check return type
  expect_s3_class(res, "tbl_df")

  # Check that p-values are sorted
  expect_true(all(diff(res$p) >= 0))
})

test_that("results() handles p-value threshold filtering", {
  # Test different p-value thresholds
  res_default <- results(cpo_example)
  res_strict <- results(cpo_example, p_threshold = 0.01)

  # Stricter threshold should return fewer rows
  expect_true(nrow(res_strict) <= nrow(res_default))
})

test_that("results() handles p-value type selection", {
  # Test both p-value types
  res_gam <- results(cpo_example, p_type = "p_gam")
  res_mvn <- results(cpo_example, p_type = "p_mvn")

  # Ensure both can be called without error
  expect_s3_class(res_gam, "tbl_df")
  expect_s3_class(res_mvn, "tbl_df")
})

test_that("results() handles minimum log fold change filtering", {
  # Test minimum log fold change filter
  res_default <- results(cpo_example)
  res_lfc <- results(cpo_example, min_lfc = 1)

  # Stricter LFC threshold should return fewer rows
  expect_true(nrow(res_lfc) <= nrow(res_default))
})

test_that("results() handles minimum count filtering", {
  # Assuming the example object has pred slot filled
  res_default <- results(cpo_example)
  res_count <- results(cpo_example, min_count = 10)

  # Stricter count threshold should return fewer rows
  expect_true(nrow(res_count) <= nrow(res_default))
})

test_that("results() handles additional column addition", {
  # Test adding log fold changes
  res_lfc <- results(cpo_example, add_lfc = TRUE)
  expect_true(any(grepl("^lfc\\.", names(res_lfc))))

  # Test adding counts
  res_counts <- results(cpo_example, add_counts = TRUE)
  expect_true(any(grepl("^counts\\.", names(res_counts))))
})


test_that("results() handles shape type selection", {
  # Test different shape types
  res_shape1 <- results(cpo_example, shape_type = "shape1")
  res_shape2 <- results(cpo_example, shape_type = "shape2")

  # Both should return a tibble
  expect_s3_class(res_shape1, "tbl_df")
  expect_s3_class(res_shape2, "tbl_df")
})

test_that("results() handles null target removal", {
  # Test removing null targets
  res_with_null <- results(cpo_example, remove_null_targets = FALSE)
  res_without_null <- results(cpo_example, remove_null_targets = TRUE)

  # Removing null targets should reduce row count
  expect_true(nrow(res_without_null) <= nrow(res_with_null))
})

