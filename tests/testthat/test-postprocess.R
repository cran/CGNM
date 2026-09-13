test_that("acceptedApproximateMinimizers returns a data.frame of parameter sets", {
  res <- fit_flipflop()

  accepted <- acceptedApproximateMinimizers(res)

  expect_s3_class(accepted, "data.frame")
  expect_equal(ncol(accepted), 3)
  expect_true(nrow(accepted) >= 1)
  expect_true(nrow(accepted) <= 50)
  expect_equal(colnames(accepted), c("Ka", "V1", "CL"))
})

test_that("acceptedIndices are a valid subset of the rows of X", {
  res <- fit_flipflop()

  idx <- acceptedIndices(res)

  expect_true(all(idx >= 1 & idx <= nrow(res$X)))
  expect_equal(length(idx), length(unique(idx)))
})

test_that("table_parameterSummary summarizes one row per parameter", {
  res <- fit_flipflop()

  tab <- table_parameterSummary(res)

  expect_equal(nrow(tab), 3)
  expect_equal(rownames(tab), c("Ka", "V1", "CL"))
  expect_true(all(c("Median") %in% colnames(tab)))
  expect_false(anyNA(tab))
})

test_that("table_parameterSummary(pretty = TRUE) shows the best fit value, with RSE% when bootstrap is available", {
  res <- fit_flipflop()

  pretty_no_boot <- table_parameterSummary(res, pretty = TRUE)
  expect_equal(rownames(pretty_no_boot), c("Ka", "V1", "CL"))
  expect_equal(colnames(pretty_no_boot), "Estimate")
  expect_false(any(grepl("%", pretty_no_boot$Estimate)))

  set.seed(7)
  boot <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    res, nonlinearFunction = flipflop_model, num_bootstrapSample = 15))

  pretty_boot <- table_parameterSummary(boot, pretty = TRUE)
  expect_equal(colnames(pretty_boot), "Estimate (RSE%)")
  expect_true(all(grepl("^[0-9.]+ \\([0-9.]+%\\)$", pretty_boot$Estimate)))

  best <- bestApproximateMinimizers(res, numParameterSet = 1)
  expect_equal(
    as.numeric(sub(" .*", "", pretty_boot$Estimate)),
    as.numeric(signif(unlist(best[rownames(pretty_boot)]), digits = 3))
  )

  # pretty = FALSE output is unaffected by this feature
  expect_false(anyNA(table_parameterSummary(res)))
})
