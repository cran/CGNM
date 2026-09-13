test_that("outputS4 defaults to FALSE and does not change the classic list output", {
  res <- fit_flipflop()

  expect_type(res, "list")
  expect_false(methods::is(res, "CGNM_result"))
})

test_that("outputS4 = TRUE returns a CGNM_result S4 object with identical field values", {
  res_list <- fit_flipflop()
  res_s4 <- as_CGNM_result_S4(res_list)

  expect_true(methods::is(res_s4, "CGNM_result"))
  # names(res_s4) lists every slot the S4 class defines (unused ones are NULL);
  # names(res_list) only lists the fields this particular function populated.
  expect_true(all(names(res_list) %in% names(res_s4)))
  expect_null(res_s4$bootstrapX)
  expect_identical(res_s4$X, res_list$X)
  expect_identical(res_s4$Y, res_list$Y)
  expect_identical(res_s4$residual_history, res_list$residual_history)
  expect_identical(res_s4$runSetting$ParameterNames, res_list$runSetting$ParameterNames)
})

test_that("as_CGNM_result_S4() is idempotent", {
  res_s4 <- as_CGNM_result_S4(fit_flipflop())
  res_s4_again <- as_CGNM_result_S4(res_s4)

  expect_identical(res_s4, res_s4_again)
})

test_that("plot_* functions accept the S4 result just like the classic list", {
  skip_if_not_installed("ggplot2")
  library(ggplot2)

  res_s4 <- as_CGNM_result_S4(fit_flipflop())

  p1 <- plot_Rank_SSR(res_s4)
  p2 <- plot_goodnessOfFit(res_s4, plotType = 1,
    independentVariableVector = c(0.1, 0.2, 0.4, 0.6, 1, 2, 3, 6, 12),
    plotRank = seq(1, 10))
  p3 <- plot_paraDistribution_byHistogram(res_s4)
  p4 <- plot_paraDistribution_byViolinPlots(res_s4)
  p5 <- plot_SSR_parameterValue(res_s4)
  p6 <- plot_parameterValue_scatterPlots(res_s4)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
})

test_that("postprocessing functions give identical results for list vs. S4 input", {
  res_list <- fit_flipflop()
  res_s4 <- as_CGNM_result_S4(res_list)

  expect_identical(acceptedApproximateMinimizers(res_s4), acceptedApproximateMinimizers(res_list))
  expect_identical(acceptedIndices(res_s4), acceptedIndices(res_list))
  expect_identical(acceptedIndices_binary(res_s4), acceptedIndices_binary(res_list))
  expect_identical(acceptedMaxSSR(res_s4), acceptedMaxSSR(res_list))
  expect_identical(topIndices(res_s4, 5), topIndices(res_list, 5))
  expect_identical(table_parameterSummary(res_s4), table_parameterSummary(res_list))
  expect_identical(
    bestApproximateMinimizers(res_s4, numParameterSet = 1),
    bestApproximateMinimizers(res_list, numParameterSet = 1)
  )
})

test_that("Cluster_Gauss_Newton_EBE_method accepts an S4 CGNM_result directly and returns S4", {
  res_s4 <- as_CGNM_result_S4(fit_flipflop())

  set.seed(11)
  ebe_s4 <- suppressWarnings(Cluster_Gauss_Newton_EBE_method(
    res_s4,
    nonlinearFunction = flipflop_model,
    individualIndices_vec = seq_len(length(flipflop_observation))
  ))

  expect_true(methods::is(ebe_s4, "CGNM_result"))
  expect_false(is.null(ebe_s4$EBE_X))
})

test_that("EBE on a classic list input still returns a classic list by default", {
  res_list <- fit_flipflop()

  set.seed(11)
  ebe_list <- suppressWarnings(Cluster_Gauss_Newton_EBE_method(
    res_list,
    nonlinearFunction = flipflop_model,
    individualIndices_vec = seq_len(length(flipflop_observation))
  ))

  expect_type(ebe_list, "list")
  expect_false(methods::is(ebe_list, "CGNM_result"))
})

test_that("Cluster_Gauss_Newton_method(outputS4 = TRUE) returns the S4 class directly", {
  set.seed(1)
  res_s4 <- suppressWarnings(Cluster_Gauss_Newton_method(
    nonlinearFunction = flipflop_model,
    targetVector = flipflop_observation,
    initial_lowerRange = rep(0.01, 3),
    initial_upperRange = rep(100, 3),
    num_minimizersToFind = 10,
    num_iteration = 3,
    saveLog = FALSE,
    outputS4 = TRUE
  ))

  expect_true(methods::is(res_s4, "CGNM_result"))
  expect_equal(dim(res_s4$X), c(10, 3))
})

test_that("an S4 result passed into the bootstrap method comes back as S4 without setting outputS4", {
  res_s4 <- as_CGNM_result_S4(fit_flipflop())

  set.seed(7)
  boot_s4 <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    res_s4,
    nonlinearFunction = flipflop_model,
    num_bootstrapSample = 10
  ))

  expect_true(methods::is(boot_s4, "CGNM_result"))
  expect_equal(dim(boot_s4$bootstrapX), c(10, 3))

  tab <- table_parameterSummary(boot_s4)
  expect_true("CGNM Bootstrap: Minimum" %in% colnames(tab))
  expect_false(anyNA(tab))
})

test_that("bootstrap on a classic list input still returns a classic list by default", {
  res_list <- fit_flipflop()

  set.seed(7)
  boot_list <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    res_list,
    nonlinearFunction = flipflop_model,
    num_bootstrapSample = 10
  ))

  expect_type(boot_list, "list")
  expect_false(methods::is(boot_list, "CGNM_result"))
})
