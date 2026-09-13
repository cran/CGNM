test_that("Cluster_Gauss_Newton_method returns the documented list fields with consistent shapes", {
  res <- fit_flipflop()

  n_par <- 3
  n_obs <- length(flipflop_observation)
  num_minimizersToFind <- 50
  num_iteration <- 25

  expect_type(res, "list")
  expect_true(all(c("X", "Y", "residual_history", "initialX", "runSetting") %in% names(res)))

  expect_equal(dim(res$X), c(num_minimizersToFind, n_par))
  expect_equal(dim(res$Y), c(num_minimizersToFind, n_obs))
  expect_equal(dim(res$initialX), c(num_minimizersToFind, n_par))
  expect_equal(nrow(res$residual_history), num_minimizersToFind)
  expect_equal(ncol(res$residual_history), num_iteration + 1)
})

test_that("the fit converges towards the known solution for the flip-flop kinetics example", {
  res <- fit_flipflop()

  ssr <- rowSums(sweep(res$Y, 2, flipflop_observation)^2)
  expect_lt(min(ssr), 0.01)

  best <- bestApproximateMinimizers(res, numParameterSet = 1)
  expect_equal(dim(best), c(1, 3))
  expect_equal(as.numeric(best[1, ]), c(0.9265, 19.072, 9.877), tolerance = 0.01)
})

test_that("saveLog = FALSE does not write anything to disk", {
  old_wd <- getwd()
  tmp_dir <- tempfile("cgnm-test-")
  dir.create(tmp_dir)
  on.exit({
    setwd(old_wd)
    unlink(tmp_dir, recursive = TRUE)
  }, add = TRUE)
  setwd(tmp_dir)

  set.seed(1)
  suppressWarnings(Cluster_Gauss_Newton_method(
    nonlinearFunction = flipflop_model,
    targetVector = flipflop_observation,
    initial_lowerRange = rep(0.01, 3),
    initial_upperRange = rep(100, 3),
    num_minimizersToFind = 10,
    num_iteration = 3,
    saveLog = FALSE
  ))

  expect_length(list.files(tmp_dir), 0)
})
