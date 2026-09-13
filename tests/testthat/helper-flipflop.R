# Shared fixture used across test files: the flip-flop kinetics model from
# the package vignette and ?Cluster_Gauss_Newton_method (an example known to
# have two distinct best-fit solutions). Fits are cached per parameter
# combination so repeated calls across test files don't re-run CGNM.

flipflop_model <- function(x) {
  observation_time <- c(0.1, 0.2, 0.4, 0.6, 1, 2, 3, 6, 12)
  Dose <- 1000
  F <- 1
  ka <- x[1]
  V1 <- x[2]
  CL_2 <- x[3]
  t <- observation_time

  Cp <- ka * F * Dose / (V1 * (ka - CL_2 / V1)) * (exp(-CL_2 / V1 * t) - exp(-ka * t))
  log10(Cp)
}

flipflop_observation <- log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))

.flipflop_fit_cache <- new.env(parent = emptyenv())

fit_flipflop <- function(num_minimizersToFind = 50, num_iteration = 25, seed = 123) {
  key <- paste(num_minimizersToFind, num_iteration, seed, sep = "-")

  if (is.null(.flipflop_fit_cache[[key]])) {
    set.seed(seed)
    .flipflop_fit_cache[[key]] <- suppressWarnings(
      Cluster_Gauss_Newton_method(
        nonlinearFunction = flipflop_model,
        targetVector = flipflop_observation,
        initial_lowerRange = rep(0.01, 3),
        initial_upperRange = rep(100, 3),
        num_minimizersToFind = num_minimizersToFind,
        num_iteration = num_iteration,
        saveLog = FALSE,
        ParameterNames = c("Ka", "V1", "CL")
      )
    )
  }

  .flipflop_fit_cache[[key]]
}
