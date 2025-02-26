res_example <- results(cpo_example)

test_that("plot_cluster() handles various input scenarios", {

  plot_result <- plot_cluster(cpo_example, res_example, changepoints = 2, shapes = "ilin")
  expect_s3_class(plot_result, "ggplot")

  expect_warning(
    result <- plot_cluster(cpo_example, res_example, changepoints = 999, shapes = "nonexistent"),
    "No targets found for the selected shapes and timepoints"
  )
  expect_null(result)

  multi_result <- plot_cluster(cpo_example, res_example,
                               changepoints = c(2),
                               shapes = c("ilin", "dlin"))
  expect_s3_class(multi_result, "ggplot")

  alpha_result <- plot_cluster(cpo_example, res_example,
                               changepoints = 2,
                               shapes = "ilin",
                               alpha = 0.5)
  expect_s3_class(alpha_result, "ggplot")

  plot_layers <- alpha_result$layers
  expect_length(plot_layers, 1)

  line_layer <- plot_layers[[1]]
  expect_equal(line_layer$aes_params$alpha, 0.5)
})

test_that("plot_cluster() handles edge cases", {
  # Test invalid alpha values
  expect_error(
    plot_cluster(cpo_example, res_example,
                 changepoints = 2,
                 shapes = "ilin",
                 alpha = 1.5),
    "alpha must be between 0 and 1"
  )

  # Test empty results tibble
  empty_res <- dplyr::tibble(
    target_id = character(0),
    cp = numeric(0),
    shape = character(0)
  )
  expect_warning(
    result <- plot_cluster(cpo_example, empty_res,
                           changepoints = 2,
                           shapes = "ilin"),
    "No targets found for the selected shapes and timepoints"
  )
  expect_null(result)
})


test_that("plot_cluster() provides informative feedback", {
  # Test message when targets are found
  expect_message(
    plot_cluster(cpo_example, res_example,
                 changepoints = 2,
                 shapes = "ilin"),
    "Plotting [0-9]+ targets"
  )
})

test_that("plot_cpam() works with gene_id", {
  expect_s3_class(
    plot_cpam(cpo_example, gene_id = "g003"),
    "ggplot"
  )
})

test_that("plot_cpam() works with target_id", {
  expect_s3_class(
    plot_cpam(cpo_example, target_id = "g003"),
    "ggplot"
  )
})

test_that("plot_cpam() handles invalid inputs", {
  expect_error(
    plot_cpam(cpo_example),
    "gene_id and target_id cannot both be null"
  )

  expect_error(
    plot_cpam(cpo_example, gene_id = "nonexistent_gene"),
    "Invalid gene_id"
  )

  expect_error(
    plot_cpam(cpo_example, target_id = "nonexistent_target"),
    "Invalid target_id"
  )
})


test_that("plot_cpam() handles various parameters", {
  expect_s3_class(
    plot_cpam(cpo_example,
              gene_id = "g003",
              cp_type = "cp_min"),
    "ggplot"
  )
  expect_s3_class(
    plot_cpam(cpo_example,
              gene_id = "g003",
              shape_type = "shape2"),
    "ggplot"
  )
  expect_s3_class(
    plot_cpam(cpo_example,
              gene_id = "g003",
              bs = "lin"),
    "ggplot"
  )
  expect_s3_class(
    plot_cpam(cpo_example,
              gene_id = "g003",
              cp_fix = 2),
    "ggplot"
  )
})

test_that("return_fits_only returns model fits", {
  # Should return a single gam object or named list of gam objects
  fits <- plot_cpam(cpo_example,
                    gene_id = "g003",
                    bs = "cv",
                    return_fits_only = TRUE)

  expect_true(inherits(fits, "scam"))
})

test_that("plot_cpam() visualization options work", {
  p1 <- plot_cpam(cpo_example,
                  gene_id = "g003",
                  show_fit = FALSE,
                  show_data = FALSE,
                  show_fit_ci = FALSE,
                  show_data_ci = FALSE)
  expect_s3_class(p1, "ggplot")

  # Test faceting
  p2 <- plot_cpam(cpo_example,
                  gene_id = "g003",
                  facet = TRUE)
  expect_s3_class(p2, "ggplot")
})

test_that("Additional parameter edge cases", {
  expect_error(
    plot_cpam(cpo_example,
              gene_id = "g003",
              cp_fix = "not_numeric"),
    "The fixed changepoint must be numeric"
  )
  expect_s3_class(
    plot_cpam(cpo_example,
              gene_id = "g003",
              base_size = 14),
    "ggplot"
  )
})
