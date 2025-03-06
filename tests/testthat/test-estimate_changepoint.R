test_that("estimate_changepoint", {
  load(system.file("extdata", "cpo_example.rda", package = "cpam"))

  targets <- cpo_example$p_table %>% dplyr::filter(.data$q_val_target <0.05) %>%
    dplyr::pull("target_id")
  cpo_example <- estimate_changepoint(cpo_example, subset = targets[1:10])

  expect_s3_class(cpo_example,"cpam")
  expect_type(cpo_example$changepoints,"list")
})

