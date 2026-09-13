# outputS4 = TRUE retains the full per-iteration X/Y history (normally
# discarded once the final result is assembled) so that profile-likelihood
# functions can be used directly on the in-memory object, with no dependency
# on saveLog/disk/working directory, and so that they behave identically to
# the classic disk-based (saveLog = TRUE) workflow when a bootstrap result is
# also available.

test_that("outputS4 = TRUE retains X_history/Y_history; default does not", {
  res_s4 <- suppressWarnings(Cluster_Gauss_Newton_method(
    nonlinearFunction = flipflop_model,
    targetVector = flipflop_observation,
    initial_lowerRange = rep(0.01, 3),
    initial_upperRange = rep(100, 3),
    num_minimizersToFind = 20,
    num_iteration = 8,
    saveLog = FALSE,
    ParameterNames = c("Ka", "V1", "CL"),
    outputS4 = TRUE
  ))

  expect_equal(length(res_s4$X_history), 9) # num_iteration + 1
  expect_equal(length(res_s4$Y_history), 9)
  expect_identical(res_s4$X_history[[1]], res_s4$initialX)
  expect_identical(res_s4$X_history[[9]], res_s4$X)

  res_list <- fit_flipflop()
  expect_false("X_history" %in% names(res_list))
})

test_that("bootstrap on an S4 result captures its own iteration history as bootstrapIterationHistory", {
  res_s4 <- suppressWarnings(Cluster_Gauss_Newton_method(
    nonlinearFunction = flipflop_model,
    targetVector = flipflop_observation,
    initial_lowerRange = rep(0.01, 3),
    initial_upperRange = rep(100, 3),
    num_minimizersToFind = 20,
    num_iteration = 8,
    saveLog = FALSE,
    ParameterNames = c("Ka", "V1", "CL"),
    outputS4 = TRUE
  ))

  boot_s4 <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    res_s4, nonlinearFunction = flipflop_model, num_bootstrapSample = 10))

  expect_true(methods::is(boot_s4, "CGNM_result"))
  expect_false(is.null(boot_s4$bootstrapIterationHistory))
  expect_false(is.null(boot_s4$bootstrapIterationHistory$X_history))

  boot_list <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
    fit_flipflop(), nonlinearFunction = flipflop_model, num_bootstrapSample = 10))
  expect_false("bootstrapIterationHistory" %in% names(boot_list))
})

test_that("profile-likelihood functions work on a saveLog = FALSE, outputS4 = TRUE result with no disk logs at all", {
  skip_if_not_installed("ggplot2")
  library(ggplot2)

  plots <- with_temp_wd({
    set.seed(21)
    res_s4 <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = FALSE,
      ParameterNames = c("Ka", "V1", "CL"),
      outputS4 = TRUE
    ))
    boot_s4 <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
      res_s4, nonlinearFunction = flipflop_model, num_bootstrapSample = 10))

    no_disk_logs <- length(list.files(".", pattern = "CGNM_log")) == 0

    list(
      no_disk_logs = no_disk_logs,
      surface = plot_SSRsurface(boot_s4),
      profile = plot_profileLikelihood(boot_s4),
      table = table_profileLikelihoodConfidenceInterval(boot_s4, silent = TRUE),
      lower = suggestInitialLowerRange(boot_s4),
      upper = suggestInitialUpperRange(boot_s4)
    )
  })

  expect_true(plots$no_disk_logs)
  expect_s3_class(plots$surface, "ggplot")
  expect_s3_class(plots$profile, "ggplot")
  expect_false(anyNA(plots$table))
})

test_that("in-memory (no disk) profile likelihood matches the classic disk-based (saveLog = TRUE) workflow", {
  results <- with_temp_wd({
    set.seed(31)
    res_disk <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = TRUE,
      ParameterNames = c("Ka", "V1", "CL")
    ))
    set.seed(32)
    boot_disk <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
      res_disk, nonlinearFunction = flipflop_model, num_bootstrapSample = 10))

    tab_disk <- table_profileLikelihoodConfidenceInterval(
      c("CGNM_log", "CGNM_log_bootstrap"), silent = TRUE)
    # passing the bootstrap-augmented object (rather than the bare "CGNM_log"
    # string) so both the main and bootstrap folders are included, matching
    # what passing the in-memory bootstrap object below does automatically.
    lower_disk <- suggestInitialLowerRange(boot_disk)
    upper_disk <- suggestInitialUpperRange(boot_disk)

    list(tab_disk = tab_disk, lower_disk = lower_disk, upper_disk = upper_disk)
  })

  in_memory <- with_temp_wd({
    set.seed(31)
    res_mem <- suppressWarnings(Cluster_Gauss_Newton_method(
      nonlinearFunction = flipflop_model,
      targetVector = flipflop_observation,
      initial_lowerRange = rep(0.01, 3),
      initial_upperRange = rep(100, 3),
      num_minimizersToFind = 20,
      num_iteration = 8,
      saveLog = FALSE,
      ParameterNames = c("Ka", "V1", "CL"),
      outputS4 = TRUE
    ))
    set.seed(32)
    boot_mem <- suppressWarnings(Cluster_Gauss_Newton_Bootstrap_method(
      res_mem, nonlinearFunction = flipflop_model, num_bootstrapSample = 10))

    tab_mem <- table_profileLikelihoodConfidenceInterval(boot_mem, silent = TRUE)
    lower_mem <- suggestInitialLowerRange(boot_mem)
    upper_mem <- suggestInitialUpperRange(boot_mem)

    list(tab_mem = tab_mem, lower_mem = lower_mem, upper_mem = upper_mem)
  })

  expect_identical(results$tab_disk, in_memory$tab_mem)
  expect_identical(results$lower_disk, in_memory$lower_mem)
  expect_identical(results$upper_disk, in_memory$upper_mem)
})
