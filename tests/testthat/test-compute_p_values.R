# tests/testthat/test-compute_p_values.R


# Helper function to create a mock cpam object
create_mock_cpam <- function(model_type = "case-only", n_times = 5, n_reps = 2, n_targets = 10) {
  # Create sample data
  times <- seq(0,(n_times-1))

  conditions <- c("case")
  if(model_type == "case-control") conditions <- c("case", "treatment")

  n_samples <- n_times * n_reps
  exp_design <- dplyr::tibble(
    sample = paste0("s", 1:n_samples),
    time = rep(times,each = n_reps),
    condition = rep_len(conditions, n_samples),
    case = as.numeric(condition=="treatment")
  )


  # Create mock count data
  target_ids <- paste0("target", 1:n_targets)
  count_matrix <- matrix(
    rpois(n_samples * n_targets, lambda = 100),
    nrow = n_targets,
    dimnames = list(target_ids, paste0("s", 1:n_samples))
  )

  # Create mock t2g mapping
  t2g <- dplyr::tibble(
    target_id = target_ids,
    gene_id = paste0("gene", ceiling(seq_along(target_ids)/2))
  )

  # Create long format data
  data_long <- dplyr::tibble(
    target_id = rep(target_ids, each = n_samples),
    sample = rep(exp_design$sample, n_targets),
    time = rep(exp_design$time, n_targets),
    condition = rep(exp_design$condition, n_targets),
    case = rep(exp_design$case, n_targets),
    counts = as.vector(count_matrix) + time * 10,
    norm_factor = 1,
    disp = 0.1
  )

  # Create mock cpam object
  structure(
    list(
      exp_design = exp_design,
      count_matrix_raw = count_matrix,
      t2g = t2g,
      data_long = data_long,
      target_to_keep = target_ids,
      model_type = model_type,
      regularize = TRUE,
      fixed_effects = NULL,
      intercept_cc = if(model_type == "case-control") "condition" else NULL,
      aggregate_to_gene = TRUE,
      num_cores = 1
    ),
    class = "cpam"
  )
}

# Test basic functionality
test_that("compute_p_values works with basic case-only model", {
  # Create mock data
  cpo <- create_mock_cpam(model_type = "case-only")

  # Run function
  result <- compute_p_values(cpo)

  # Tests
  expect_s3_class(result, "cpam")
  expect_true(!is.null(result$p_table))
  expect_true(all(c("target_id", "gene_id", "counts_mean",
                    "p_val_target", "q_val_target") %in%
                    names(result$p_table)))
  expect_true(all(result$p_table$p_val_target >= 0 &
                    result$p_table$p_val_target <= 1))
  expect_true(all(result$p_table$q_val_target >= 0 &
                    result$p_table$q_val_target <= 1))
})

# Test case-control model
test_that("compute_p_values works with case-control model", {
  # Create mock data
  cpo <- create_mock_cpam(model_type = "case-control")

  # Run function
  result <- compute_p_values(cpo)

  # Tests
  expect_s3_class(result, "cpam")
  expect_true(!is.null(result$p_table))
  expect_true(!is.null(result$p_table$counts_mean))
})

# Test subset functionality
test_that("compute_p_values handles subset parameter correctly", {
  cpo <- create_mock_cpam()
  subset_targets <- c("target1", "target2")

  # Test valid subset
  result <- compute_p_values(cpo, subset = subset_targets)
  expect_equal(sort(unique(result$p_table$target_id)), sort(subset_targets))

  # Test invalid subset
  expect_error(
    compute_p_values(cpo, subset = c("invalid_target1", "invalid_target2")),
    "Subset contains invalid targets"
  )

  # Test invalid subset type
  expect_error(
    compute_p_values(cpo, subset = 1:2),
    "subset must be character vector of target_id values"
  )
})

# Test p-value adjustment methods
test_that("compute_p_values handles different p-value adjustment methods", {
  cpo <- create_mock_cpam(n_targets = 30)

  # Test different adjustment methods
  result_bh <- compute_p_values(cpo, p_adj_method = "BH")
  result_bonf <- compute_p_values(cpo, p_adj_method = "bonferroni")

  expect_false(identical(
    result_bh$p_table$q_val_target,
    result_bonf$p_table$q_val_target
  ))
})

# Test gene-level aggregation
test_that("compute_p_values deals gene-level aggregated counts correctly", {
  cpo <- create_mock_cpam()
  cpo$aggregate_to_gene <- TRUE

  result <- compute_p_values(cpo)

  expect_true(all(c("p_val_gene", "q_val_gene") %in% names(result$p_table)))
  expect_true(length(unique(result$p_table$gene_id)) <=
                length(unique(result$p_table$target_id)))
})

# Test error handling
test_that("compute_p_values handles errors gracefully", {
  cpo <- create_mock_cpam()

  # Test with invalid data
  cpo$data_long$counts <- NA
  result <- compute_p_values(cpo)
  expect_true(all(is.na(result$p_table$p_val_target)))

})

# Test parallel processing
test_that("compute_p_values works with different numbers of cores", {
  skip_on_cran()
  cpo <- create_mock_cpam()

  # Test with multiple cores
  cpo$num_cores <- 2
  result_parallel <- compute_p_values(cpo)

  # Test with single core
  cpo$num_cores <- 1
  result_serial <- compute_p_values(cpo)

  # Results should be the same
  expect_equal(
    result_parallel$p_table$p_val_target,
    result_serial$p_table$p_val_target,
    tolerance = 1e-6
  )
})

# Test fixed effects
test_that("compute_p_values handles fixed effects correctly", {
  cpo <- create_mock_cpam()
  cpo$fixed_effects <- "batch"
  cpo$data_long$batch <- factor(rep(c("A", "B"), length.out = nrow(cpo$data_long)))

  result <- compute_p_values(cpo)
  expect_true(!is.null(result$p_table))

  cpo <- create_mock_cpam(model_type = "case-control")
  cpo$fixed_effects <- "batch"
  cpo$data_long$batch <- factor(rep(c("A", "B"), length.out = nrow(cpo$data_long)))

  result <- compute_p_values(cpo)
  expect_true(!is.null(result$p_table))
})

