

#' Generate minimal test data for cpam testing
#'
#' @param n_samples Number of samples
#' @param n_genes Number of genes
#' @param n_transcripts Number of transcripts per gene
#' @return A list containing test data structures
#'
generate_test_data <- function(n_samples = 9, n_genes = 10, n_transcripts = 2) {
  # Create experimental design with minimal required columns
  exp_design <- dplyr::tibble(
    sample = paste0("sample_", 1:n_samples),
    time = rep(c(0, 2, 5), each = n_samples/3),
    path = paste0("path/to/sample_", 1:n_samples, "/abundance.h5")
  )

  # Create transcript to gene mapping
  t2g <- dplyr::tibble(
    target_id = paste0("transcript_", 1:(n_genes * n_transcripts)),
    gene_id = rep(paste0("gene_", 1:n_genes), each = n_transcripts)
  )

  # Create a small synthetic count matrix
  count_matrix <- matrix(
    rpois(n_samples * n_genes * n_transcripts, lambda = 10),
    nrow = n_genes * n_transcripts,
    ncol = n_samples,
    dimnames = list(
      paste0("transcript_", 1:(n_genes * n_transcripts)),
      paste0("sample_", 1:n_samples)
    )
  )

  list(
    exp_design = exp_design,
    t2g = t2g,
    count_matrix = count_matrix
  )
}

test_that("prepare_cpam validation functions work correctly", {
  test_data <- generate_test_data()

  # Test validate_inputs function
  expect_error(
    validate_inputs(
      exp_design = test_data$exp_design[, -which(names(test_data$exp_design) == "time")],
      count_matrix = test_data$count_matrix,
      t2g = test_data$t2g,
      import_type = NULL,
      model_type = "case-only",
      condition_var = "condition",
      case_value = "treatment",
      fixed_effects = NULL,
      aggregate_to_gene = TRUE,
      gene_level = FALSE,
      import = FALSE
    ),
    "'exp_design' must include a 'time' column"
  )

  # Test validate_cores function
  expect_warning(
    result <- validate_cores(1000)
  )
  expect_true(result < 1000)
})

test_that("prepare_cpam helper functions process data correctly", {
  test_data <- generate_test_data()

  # Test convert_to_long_format
  long_data <- convert_to_long_format(test_data$count_matrix, test_data$exp_design)
  expect_equal(nrow(long_data), nrow(test_data$count_matrix) * ncol(test_data$count_matrix))
  expect_true(all(c("target_id", "sample", "counts_raw", "counts", "time") %in% names(long_data)))

  # Test handle_gene_aggregation
  agg_data <- handle_gene_aggregation(long_data, test_data$t2g, aggregate_to_gene = TRUE)
  expect_true("gene_id" %in% names(agg_data))
  expect_equal(
    length(unique(agg_data$gene_id)),
    length(unique(test_data$t2g$gene_id))
  )

  # Test compute_normalization_factors
  norm_factors <- compute_normalization_factors(test_data$count_matrix, TRUE)
  expect_equal(length(norm_factors), ncol(test_data$count_matrix))
  expect_equal(names(norm_factors), colnames(test_data$count_matrix))
})



#' prepare_cpam with data and mimimal computation
#'
test_that("prepare_cpam with case-only data and minimal computation", {
  test_data <- generate_test_data()


      # Run the function with minimal computation
      result <- prepare_cpam(
        exp_design = test_data$exp_design,
        count_matrix = test_data$count_matrix,
        t2g = test_data$t2g,
        bootstrap = FALSE,
        regularize = FALSE,
        gene_level = FALSE
      ) %>% suppressMessages()

      # Validate the result structure
      expect_s3_class(result, "cpam")
      expect_equal(result$model_type, "case-only")
      expect_equal(result$gene_level, FALSE)
      expect_equal(result$bootstrap, FALSE)
      expect_equal(result$regularize, FALSE)

})




#' prepare_cpam with case-control data and mimimal computation
#'
test_that("prepare_cpam with case-control data and minimal computation", {
  test_data <- generate_test_data()
  test_data$exp_design$condition <- rep(c("treatment", "control"), length.out = nrow(test_data$exp_design))


  # Run the function with minimal computation
  result <- prepare_cpam(
    exp_design = test_data$exp_design,
    count_matrix = test_data$count_matrix,
    t2g = test_data$t2g,
    model_type = "case-control",
    condition_var = "condition",
    case_value = "treatment",
    bootstrap = FALSE,
    regularize = FALSE,
    gene_level = FALSE
  ) %>% suppressMessages()

  # Verify case-control specific outputs
  expect_equal(result$model_type, "case-control")
  expect_equal(result$condition_var, "condition")
  expect_equal(result$case_value, "treatment")
  expect_true("case" %in% names(result$exp_design))
  expect_equal(sum(result$exp_design$case), sum(test_data$exp_design$condition == "treatment"))

})

