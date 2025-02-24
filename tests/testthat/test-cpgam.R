create_test_data_cpgam <- function(model_type = "case-only", n_times = 6, n_reps = 4) {
  # Create sample data
  times <- seq(0,(n_times-1))

  conditions <- c("case")
  if(model_type == "case-control") conditions <- c("case", "treatment")

  n_samples <- n_times * n_reps
  dplyr::tibble(
    sample = paste0("s", 1:n_samples),
    time = rep(times,each = n_reps),
    condition = rep_len(conditions, n_samples),
    case = as.numeric(condition=="treatment"),
    counts = rpois(n_samples, lambda = 100),
    disp = 0.1,
    norm_factor = 1
  )

}

test_that("cpgam handles case-only model", {
  # Generate test data
  data <- create_test_data_cpgam(model = "case-only")

  # Fit model
  fit <- cpgam(
    data = data,
    cp = 1,
    regularize = TRUE,
    model_type = "case-only",
    family = "nb"
  )

  # Tests
  expect_s3_class(fit, "gam")
  expect_equal(fit$cp, 1)
  expect_equal(fit$model_type, "case-only")
  expect_equal(fit$bs, "tp")  # default basis
  expect_true(!is.null(fit$data))
})


test_that("cpgam handles case-control model", {
  data <- create_test_data_cpgam(model = "case-control")

  fit <- cpgam(
    data = data,
    cp = 0,
    regularize = TRUE,
    model_type = "case-control",
    family = "nb"
  )

  expect_s3_class(fit, "gam")
  expect_true("time" %in% names(fit$data))
  expect_true("case" %in% names(fit$data))
  expect_equal(fit$model_type, "case-control")
})


test_that("cpgam handles different basis types", {
  data <- create_test_data_cpgam()

  # Test GAM bases
  gam_bases <- c("tp", "null", "lin", "ilin", "dlin")
  for(bs in gam_bases) {
    fit <- cpgam(
      data = data,
      cp = 1,
      regularize = TRUE,
      bs = bs
    )
    expect_s3_class(fit, "gam")
    expect_equal(fit$bs, bs)
  }

  # Test SCAM bases
  scam_bases <- c("micv", "mdcx", "cv", "cx", "micx", "mdcv")
  for(bs in scam_bases) {
    fit <- cpgam(
      data = data,
      cp = 1,
      regularize = TRUE,
      bs = bs
    )
    expect_s3_class(fit, "scam")
    expect_equal(fit$bs, bs)
  }
})

test_that("cpgam handles gaussian family", {
  data <- create_test_data_cpgam()

  fit <- cpgam(
    data = data,
    cp = 1,
    regularize = TRUE,
    family = "gaussian"
  )

  expect_s3_class(fit, "gam")
  expect_equal(fit$family$family, "gaussian")
})

test_that("cpgam handles fixed effects", {
  data <- create_test_data_cpgam()
  data$batch <- factor(rep(c("A", "B"), length.out = nrow(data)))

  fit <- cpgam(
    data = data,
    cp = 1,
    regularize = TRUE,
    fixed_effects = "batch"
  )

  expect_s3_class(fit, "gam")
  expect_true(grepl("batch", deparse(fit$f)))
})


test_that("cpgam handles sp", {
  data <- create_test_data_cpgam()

  fit <- cpgam(
    data = data,
    cp = 1,
    regularize = TRUE,
    sp = 0.1
  )

  expect_s3_class(fit, "gam")
  expect_equal(as.numeric(fit$full.sp), 0.1)
})


test_that("cpgam handles errors gracefully", {
  data <- create_test_data_cpgam()

  # Test with invalid basis
  expect_error(
    cpgam(
      data = data,
      cp = 5,
      regularize = TRUE,
      bs = "invalid"
    )
  )

  # Test with missing required columns
  bad_data <- data
  bad_data$time <- NULL
  expect_error(
    cpgam(
      data = bad_data,
      cp = 5,
      regularize = TRUE
    )
  )
})

test_that("cpgam works without regularize", {
  data <- create_test_data_cpgam()

  fit <- cpgam(
    data = data,
    cp = 0,
    gam_optimizer = "efs",
    regularize = FALSE,
  );fit

  expect_s3_class(fit, "gam")

  expect_warning(cpgam(
    data = data,
    cp = 1,
    bs = "cv",
    regularize = FALSE
  ))
})




