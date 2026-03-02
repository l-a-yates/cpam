# test-pipeline.R
# Integration tests for the full cpam pipeline.
# These test the PUBLIC API contract — they should NOT need updating
# during internal refactors unless the user-facing behavior changes.


# --- Helper: run the full pipeline on example data ---
run_pipeline <- function(subset_genes = paste0("g", sprintf("%03d", 1:30)),
                         q_threshold = 0.05) {
  load(system.file("extdata", "count_matrix_example.rda", package = "cpam"))
  load(system.file("extdata", "exp_design_example.rda", package = "cpam"))

  cpo <- prepare_cpam(
    exp_design = exp_design_example,
    count_matrix = count_matrix_example[subset_genes, ],
    gene_level = TRUE,
    num_cores = 1
  )

  cpo <- compute_p_values(cpo)

  sig <- cpo$p_table$target_id[cpo$p_table$q_val_target < q_threshold]
  if (length(sig) > 0) {
    set.seed(42)
    cpo <- estimate_changepoint(cpo, subset = sig)
    cpo <- select_shape(cpo)
  }

  cpo
}


# ====================================================================
# 1. Structural validation of the pre-computed example object
# ====================================================================

test_that("cpo_example has expected structure", {
  load(system.file("extdata", "cpo_example.rda", package = "cpam"))

  expect_s3_class(cpo_example, "cpam")

  # Core data slots exist
  expect_true(!is.null(cpo_example$exp_design))
  expect_true(!is.null(cpo_example$count_matrix_raw))
  expect_true(!is.null(cpo_example$data_long))
  expect_true(!is.null(cpo_example$target_to_keep))

  # Pipeline result slots exist
  expect_true(!is.null(cpo_example$p_table))
  expect_true(!is.null(cpo_example$changepoints))
  expect_true(!is.null(cpo_example$shapes))
})

test_that("cpo_example results() returns valid output", {
  load(system.file("extdata", "cpo_example.rda", package = "cpam"))
  res <- results(cpo_example)

  expect_s3_class(res, "tbl_df")
  expect_true(nrow(res) > 0)
  expect_true(all(c("target_id", "p", "cp", "shape") %in% names(res)))

  # p-values sorted ascending
  expect_equal(res$p, sort(res$p))

  # Values in valid ranges
  expect_true(all(res$p >= 0 & res$p <= 1))
  expect_true(all(res$cp %in% cpo_example$times))
})


# ====================================================================
# 2. Full pipeline from raw data — structural invariants
# ====================================================================

test_that("prepare_cpam returns valid structure", {
  load(system.file("extdata", "count_matrix_example.rda", package = "cpam"))
  load(system.file("extdata", "exp_design_example.rda", package = "cpam"))

  subset_genes <- paste0("g", sprintf("%03d", 1:30))
  cpo <- prepare_cpam(
    exp_design = exp_design_example,
    count_matrix = count_matrix_example[subset_genes, ],
    gene_level = TRUE,
    num_cores = 1
  )

  expect_s3_class(cpo, "cpam")
  expect_equal(cpo$model_type, "case-only")
  expect_true(cpo$regularize)
  expect_true(length(cpo$target_to_keep) > 0)
  expect_true(length(cpo$target_to_keep) <= 30)

  # data_long has required columns
  expect_true(all(c("target_id", "sample", "counts", "time", "norm_factor")
                  %in% names(cpo$data_long)))

  # Count matrix dimensions consistent with design
  expect_equal(ncol(cpo$count_matrix_raw), nrow(cpo$exp_design))
})

test_that("compute_p_values returns valid p-values", {
  cpo <- run_pipeline()

  expect_true(!is.null(cpo$p_table))
  expect_true(all(c("target_id", "p_val_target", "q_val_target")
                  %in% names(cpo$p_table)))

  # All kept targets have p-values
  expect_equal(nrow(cpo$p_table), length(cpo$target_to_keep))

  # p-values in valid range
  pvals <- cpo$p_table$p_val_target
  expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))

  # q-values >= p-values (BH adjustment inflates)
  qvals <- cpo$p_table$q_val_target
  non_na <- !is.na(pvals) & !is.na(qvals)
  expect_true(all(qvals[non_na] >= pvals[non_na] - 1e-15))
})

test_that("estimate_changepoint returns valid changepoints", {
  cpo <- run_pipeline()

  if (is.null(cpo$changepoints)) skip("No significant targets for changepoint estimation")

  expect_s3_class(cpo$changepoints, "tbl_df")
  expect_true(all(c("target_id", "cp_min", "cp_1se") %in% names(cpo$changepoints)))

  # Changepoints within observed time range
  time_range <- range(cpo$exp_design$time)
  expect_true(all(cpo$changepoints$cp_min >= time_range[1] &
                    cpo$changepoints$cp_min <= time_range[2]))
  expect_true(all(cpo$changepoints$cp_1se >= time_range[1] &
                    cpo$changepoints$cp_1se <= time_range[2]))
})

test_that("select_shape returns valid shapes", {
  cpo <- run_pipeline()

  if (is.null(cpo$shapes)) skip("No significant targets for shape selection")

  expect_s3_class(cpo$shapes, "tbl_df")
  expect_true(all(c("target_id", "cp", "shape1", "shape2") %in% names(cpo$shapes)))

  valid_shapes <- c("micv", "mdcx", "cv", "cx", "lin", "ilin", "dlin", "tp", "null")
  expect_true(all(cpo$shapes$shape1 %in% valid_shapes))
  expect_true(all(cpo$shapes$shape2 %in% valid_shapes))
})

test_that("results returns valid output from fresh pipeline", {
  cpo <- run_pipeline()

  if (is.null(cpo$shapes)) skip("No significant targets for results")

  res <- results(cpo)
  expect_s3_class(res, "tbl_df")
  expect_true(nrow(res) > 0)
  expect_true(all(c("target_id", "p", "cp", "shape") %in% names(res)))

  # p-values sorted ascending
  expect_equal(res$p, sort(res$p))
})


# ====================================================================
# 3. Cross-step consistency
# ====================================================================

test_that("targets flow consistently through pipeline steps", {
  cpo <- run_pipeline()

  if (is.null(cpo$shapes)) skip("No significant targets")

  # Every target in shapes must be in changepoints
  expect_true(all(cpo$shapes$target_id %in% cpo$changepoints$target_id))

  # Every target in changepoints must be in p_table
  expect_true(all(cpo$changepoints$target_id %in% cpo$p_table$target_id))

  # Changepoints used in shapes must match
  for (tid in cpo$shapes$target_id) {
    shape_cp <- cpo$shapes$cp[cpo$shapes$target_id == tid]
    cp_1se <- cpo$changepoints$cp_1se[cpo$changepoints$target_id == tid]
    expect_equal(shape_cp, cp_1se)
  }
})


# ====================================================================
# 4. Numerical reference values (snapshot)
# ====================================================================

test_that("pipeline p-values match reference", {
  cpo <- run_pipeline()

  expect_snapshot_value(
    round(cpo$p_table$p_val_target, digits = 8),
    style = "serialize"
  )
})

test_that("pipeline changepoints match reference", {
  cpo <- run_pipeline()

  if (is.null(cpo$changepoints)) skip("No significant targets")

  expect_snapshot_value(
    cpo$changepoints[, c("target_id", "cp_min", "cp_1se")],
    style = "serialize"
  )
})

test_that("pipeline shapes match reference", {
  cpo <- run_pipeline()

  if (is.null(cpo$shapes)) skip("No significant targets")

  expect_snapshot_value(
    cpo$shapes[, c("target_id", "cp", "shape1", "shape2")],
    style = "serialize"
  )
})

test_that("pipeline final results match reference", {
  cpo <- run_pipeline()

  if (is.null(cpo$shapes)) skip("No significant targets")

  res <- results(cpo)

  expect_snapshot_value(
    data.frame(
      target_id = res$target_id,
      p = round(res$p, 8),
      cp = res$cp,
      shape = res$shape
    ),
    style = "serialize"
  )
})
