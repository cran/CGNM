test_that("Cluster_Gauss_Newton_Bootstrap_method adds bootstrap fields to the result", {
  res <- fit_flipflop()

  set.seed(42)
  boot <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    res,
    nonlinearFunction = flipflop_model,
    num_bootstrapSample = 15
  ))

  expect_true(all(c("bootstrapX", "bootstrapY") %in% names(boot)))
  expect_equal(dim(boot$bootstrapX), c(15, 3))
  expect_false(anyNA(boot$bootstrapX))
})

test_that("table_parameterSummary switches to the bootstrap columns once bootstrapX is present", {
  res <- fit_flipflop()

  set.seed(42)
  boot <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    res,
    nonlinearFunction = flipflop_model,
    num_bootstrapSample = 15
  ))

  tab <- table_parameterSummary(boot)

  expect_equal(nrow(tab), 3)
  expect_equal(rownames(tab), c("Ka", "V1", "CL"))
  expect_true("CGNM Bootstrap: Minimum" %in% colnames(tab))
  expect_false(anyNA(tab))
})
