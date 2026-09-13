# CGNM: Cluster Gauss-Newton Method

An R package to find **multiple** approximate minimizers of a nonlinear least
squares problem

```
argmin_x || f(x) - y* ||
```

without assuming the minimizer is unique. In the context of model fitting,
`f` is the model, `x` is the parameter to estimate, and `y*` is the observed
data. Because CGNM searches from a *range* of initial iterates rather than a
single starting point, and returns many approximate minimizers at once, it
can also reveal when a model's parameters are not practically identifiable
(the returned minimizers won't converge to a single point).

See the papers below for the algorithm and comparisons with conventional
multi-start optimization:

- Aoki, Y., Hayami, K., Toshimoto, K., & Sugiyama, Y. (2020). Cluster
  Gauss–Newton method. *Optimization and Engineering*, 1-31.
  [doi:10.1007/s11081-020-09571-2](https://doi.org/10.1007/s11081-020-09571-2)
- Aoki, Y. & Sugiyama, Y. (2024). Cluster Gauss-Newton method for a quick
  approximation of profile likelihood: With application to
  physiologically-based pharmacokinetic models. *CPT Pharmacometrics Syst
  Pharmacol*. 13(1):54-67.
  [doi:10.1002/psp4.13055](https://doi.org/10.1002/psp4.13055)

If you use CGNM in your research, please cite the relevant paper(s) above.

## When to use CGNM

- You want to fit a relatively complex parameterized model to data but are
  unsure of a good initial guess — CGNM takes an initial *range* instead of a
  point.
- You want to know whether the fit is practically identifiable (one best-fit
  parameter set, or many) — CGNM returns multiple approximate minimizers so
  non-identifiability shows up directly.
- Your model is a black box (not necessarily continuous or differentiable in
  its parameters) — CGNM makes no assumptions about the model's internal form.

## When not to use CGNM

- You already know roughly where the best-fit parameter is and just need one
  answer — plain `nls()` will be faster.
- The model is cheap to evaluate and computational cost isn't a concern —
  CGNM's main advantage is reducing the number of model evaluations needed
  compared to repeated multi-start `nls()`.

## Installation

```r
# from a local checkout of this repository
install.packages("path/to/CGNM", repos = NULL, type = "source")

# or, from within the checked-out directory
# R CMD INSTALL .
```

## Quick example

```r
library(CGNM)

# flip-flop kinetics: a model known to have two distinct best-fit solutions
model_function <- function(x) {
  observation_time <- c(0.1, 0.2, 0.4, 0.6, 1, 2, 3, 6, 12)
  Dose <- 1000; F <- 1
  ka <- x[1]; V1 <- x[2]; CL_2 <- x[3]
  t <- observation_time
  Cp <- ka * F * Dose / (V1 * (ka - CL_2 / V1)) * (exp(-CL_2 / V1 * t) - exp(-ka * t))
  log10(Cp)
}

observation <- log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))

CGNM_result <- Cluster_Gauss_Newton_method(
  nonlinearFunction  = model_function,
  targetVector       = observation,
  initial_lowerRange = rep(0.01, 3),
  initial_upperRange = rep(100, 3),
  saveLog            = FALSE
)

acceptedApproximateMinimizers(CGNM_result)
```

Then inspect the fit visually:

```r
library(ggplot2)
plot_Rank_SSR(CGNM_result)
plot_goodnessOfFit(CGNM_result, plotType = 1,
                    independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12),
                    plotRank = seq(1, 50))
```

![goodness of fit](goodness_of_fit.png)
![parameter distribution](param_histogram.png)
![rank vs SSR](rank_ssr.png)

For the full workflow — fitting, residual-resampling bootstrap, accepted
minimizers, summary tables, and profile likelihood — see
`vignette("CGNM-vignette", package = "CGNM")` or
[vignettes/CGNM-vignette.Rmd](vignettes/CGNM-vignette.Rmd).

## Documentation

- `?CGNM` — package overview and workflow map.
- `?Cluster_Gauss_Newton_method` — the main fitting function.

By default, `Cluster_Gauss_Newton_method()` and the bootstrap/EBE variants
return a plain list. Pass `outputS4 = TRUE` to get a `CGNM_result` S4 object
instead (see `?"CGNM_result-class"`) — every field is still reachable via
`$`/`[[` exactly as on the list, so this is a drop-in, fully backward-
compatible opt-in.


## License

MIT. See [LICENSE](LICENSE).
