# plot_profileLikelihood, compare_profileLikelihood,
# table_profileLikelihoodConfidenceInterval, plot_SSRsurface,
# suggestInitialLowerRange, and suggestInitialUpperRange all read the
# iteration-by-iteration log CGNM writes to disk when saveLog = TRUE, so
# unlike the other postprocessing functions they need a real saveLog = TRUE
# fit (not just an in-memory result) to exercise.

test_that("table_profileLikelihoodConfidenceInterval gives identical results for a folder path, a classic list, and an S4 result", {
  skip_if_not_installed("ggplot2")

  result <- with_temp_wd({
    set.seed(1)
    res_list <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = TRUE,
      ParameterNames = c("Ka", "V1", "CL")
    ))
    res_s4 <- as_CGNM_result_S4(res_list)

    list(
      path = table_profileLikelihoodConfidenceInterval("CGNM_log", silent = TRUE),
      list = table_profileLikelihoodConfidenceInterval(res_list, silent = TRUE),
      s4   = table_profileLikelihoodConfidenceInterval(res_s4, silent = TRUE)
    )
  })

  expect_identical(result$path, result$list)
  expect_identical(result$path, result$s4)
})

test_that("suggestInitialLowerRange/UpperRange give identical results for a folder path and an S4 result", {
  result <- with_temp_wd({
    set.seed(1)
    res_list <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = TRUE,
      ParameterNames = c("Ka", "V1", "CL")
    ))
    res_s4 <- as_CGNM_result_S4(res_list)

    list(
      lower_path = suggestInitialLowerRange("CGNM_log"),
      lower_s4   = suggestInitialLowerRange(res_s4),
      upper_path = suggestInitialUpperRange("CGNM_log"),
      upper_s4   = suggestInitialUpperRange(res_s4)
    )
  })

  expect_identical(result$lower_path, result$lower_s4)
  expect_identical(result$upper_path, result$upper_s4)
})

test_that("plot_SSRsurface and plot_profileLikelihood accept an S4 result", {
  skip_if_not_installed("ggplot2")
  library(ggplot2)

  plots <- with_temp_wd({
    set.seed(1)
    res_s4 <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = TRUE,
      ParameterNames = c("Ka", "V1", "CL"),
      outputS4 = TRUE
    ))

    list(
      surface = plot_SSRsurface(res_s4),
      profile = plot_profileLikelihood(res_s4)
    )
  })

  expect_s3_class(plots$surface, "ggplot")
  expect_s3_class(plots$profile, "ggplot")
})

test_that("compare_profileLikelihood accepts a list mixing classic-list and S4 results", {
  skip_if_not_installed("ggplot2")
  library(ggplot2)

  plot <- with_temp_wd({
    set.seed(1)
    res_list <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = TRUE,
      runName = "listrun",
      ParameterNames = c("Ka", "V1", "CL")
    ))
    res_s4 <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = TRUE,
      runName = "s4run",
      ParameterNames = c("Ka", "V1", "CL"),
      outputS4 = TRUE
    ))

    compare_profileLikelihood(list(res_list, res_s4))
  })

  expect_s3_class(plot, "ggplot")
})
