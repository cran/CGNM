#' @title CGNM_result-class
#' @name CGNM_result-class
#' @aliases CGNM_result-class
#' @description
#' An optional, opt-in S4 alternative to the classic plain-\code{list} output
#' of \code{\link{Cluster_Gauss_Newton_method}},
#' \code{\link{Cluster_Gauss_Newton_Bootstrap_method}}, and
#' \code{\link{Cluster_Gauss_Newton_EBE_method}} (set \code{outputS4 = TRUE}
#' on any of those functions to receive one instead of a list).
#'
#' Every field that the classic list output ever contains (\code{X},
#' \code{Y}, \code{residual_history}, \code{initialX}, \code{initialY},
#' \code{lambda_history}, \code{finalParameterCombinations}, \code{runSetting},
#' the bootstrap fields, and the EBE fields) is a slot here, and \code{$},
#' \code{$<-}, \code{[[}, \code{[[<-}, and \code{names()} all work exactly as
#' they do on the classic list. This means every existing postprocessing or
#' plotting function in this package (\code{acceptedApproximateMinimizers()},
#' \code{table_parameterSummary()}, \code{plot_goodnessOfFit()}, etc.) accepts
#' a \code{CGNM_result} S4 object as a drop-in replacement for the list, with
#' no change in behavior. The only difference from the classic list is that
#' \code{inherits(x, "CGNM_result")}/\code{is(x, "CGNM_result")} lets calling
#' code (including automated tooling) confirm the object's identity, which a
#' plain list cannot do.
#'
#' Existing code that relies on the classic list output is unaffected:
#' \code{outputS4} defaults to \code{FALSE} everywhere.
#'
#' When \code{outputS4 = TRUE}, two further slots are populated so that
#' profile-likelihood functions (\code{\link{plot_profileLikelihood}},
#' \code{\link{compare_profileLikelihood}},
#' \code{\link{table_profileLikelihoodConfidenceInterval}},
#' \code{\link{plot_SSRsurface}}, \code{\link{suggestInitialLowerRange}},
#' \code{\link{suggestInitialUpperRange}}) can be used directly on the object
#' with no dependency on \code{saveLog}/disk/working directory:
#' \code{X_history}/\code{Y_history} hold every intermediate iteration of the
#' main fit (normally discarded once the final result is assembled), and
#' \code{bootstrapIterationHistory} holds the equivalent history from
#' \code{\link{Cluster_Gauss_Newton_Bootstrap_method}}, if one was run. Both
#' are \code{NULL} unless \code{outputS4 = TRUE} was used, so the memory cost
#' of retaining them is opt-in.
#' @param x a \code{CGNM_result} object.
#' @param name \emph{string} the slot/field name, as used with \code{$} and \code{$<-} (e.g. \code{"X"}).
#' @param i \emph{string} the slot/field name, as used with \code{[[} and \code{[[<-} (e.g. \code{"X"}).
#' @param j not used; present for consistency with the \code{[[} generic.
#' @param ... not used; present for consistency with the \code{[[} generic.
#' @param value the replacement value, as used with \code{$<-} and \code{[[<-}.
#' @exportClass CGNM_result
setClass(
  "CGNM_result",
  representation(
    X = "ANY",
    Y = "ANY",
    residual_history = "ANY",
    lambda_history = "ANY",
    initialX = "ANY",
    initialY = "ANY",
    finalParameterCombinations = "ANY",
    runSetting = "ANY",
    bootstrapX = "ANY",
    bootstrapY = "ANY",
    bootstrapTheta = "ANY",
    bootstrapParameterCombinations = "ANY",
    EBE_X = "ANY",
    EBE_Y = "ANY",
    EBE_ParameterCombinations = "ANY",
    EBE_ParameterCombinations_list = "ANY",
    X_history = "ANY",
    Y_history = "ANY",
    bootstrapIterationHistory = "ANY"
  ),
  prototype(
    X = NULL,
    Y = NULL,
    residual_history = NULL,
    lambda_history = NULL,
    initialX = NULL,
    initialY = NULL,
    finalParameterCombinations = NULL,
    runSetting = NULL,
    bootstrapX = NULL,
    bootstrapY = NULL,
    bootstrapTheta = NULL,
    bootstrapParameterCombinations = NULL,
    EBE_X = NULL,
    EBE_Y = NULL,
    EBE_ParameterCombinations = NULL,
    EBE_ParameterCombinations_list = NULL,
    X_history = NULL,
    Y_history = NULL,
    bootstrapIterationHistory = NULL
  )
)

#' @rdname CGNM_result-class
methods::setMethod("$", "CGNM_result", function(x, name) {
  methods::slot(x, name)
})

#' @rdname CGNM_result-class
methods::setMethod("$<-", "CGNM_result", function(x, name, value) {
  methods::slot(x, name) <- value
  x
})

#' @rdname CGNM_result-class
methods::setMethod("[[", "CGNM_result", function(x, i, ...) {
  methods::slot(x, i)
})

#' @rdname CGNM_result-class
methods::setMethod("[[<-", "CGNM_result", function(x, i, ..., value) {
  methods::slot(x, i) <- value
  x
})

#' @rdname CGNM_result-class
methods::setMethod("names", "CGNM_result", function(x) {
  methods::slotNames(x)
})

methods::setMethod("show", "CGNM_result", function(object) {
  dims <- function(m) if (is.null(m)) "not set" else paste(dim(m), collapse = " x ")
  cat("<CGNM_result> (S4; access fields with $ or [[ ]] exactly as with the classic list output)\n")
  cat("  X (parameters):   ", dims(object@X), "\n")
  cat("  Y (model output): ", dims(object@Y), "\n")
  cat("  residual_history: ", dims(object@residual_history), "\n")
  if (!is.null(object@bootstrapX)) cat("  bootstrapX:       ", dims(object@bootstrapX), "\n")
  if (!is.null(object@EBE_X)) cat("  EBE_X:            ", dims(object@EBE_X), "\n")
  invisible(object)
})

#' Convert a classic CGNM result list into a \code{CGNM_result} S4 object
#'
#' @description
#' \code{Cluster_Gauss_Newton_method()}, \code{Cluster_Gauss_Newton_Bootstrap_method()},
#' and \code{Cluster_Gauss_Newton_EBE_method()} all accept an \code{outputS4 = TRUE}
#' argument that does this conversion automatically. Use \code{as_CGNM_result_S4()}
#' directly when you have an existing classic list result (for example one
#' loaded from an older saved \code{.RDATA} file) that you want to convert
#' after the fact. Calling it on an object that is already a \code{CGNM_result}
#' S4 object returns that object unchanged.
#' @param CGNM_result_list (required input) \emph{a list, or a CGNM_result S4 object}
#' the classic list returned by \code{\link{Cluster_Gauss_Newton_method}} (or
#' the bootstrap/EBE variants).
#' @return a \code{CGNM_result} S4 object exposing the same fields via \code{$}/\code{[[}.
#' @examples
#' model_function <- function(x) x
#' CGNM_result <- Cluster_Gauss_Newton_method(
#'   nonlinearFunction = model_function, targetVector = c(1, 2, 3),
#'   initial_lowerRange = rep(0.1, 3), initial_upperRange = rep(10, 3),
#'   num_iteration = 2, num_minimizersToFind = 5, saveLog = FALSE)
#' CGNM_result_S4 <- as_CGNM_result_S4(CGNM_result)
#' CGNM_result_S4$X
#' @export
as_CGNM_result_S4 <- function(CGNM_result_list) {
  if (methods::is(CGNM_result_list, "CGNM_result")) {
    return(CGNM_result_list)
  }

  slotNames_available <- methods::slotNames("CGNM_result")
  present <- intersect(names(CGNM_result_list), slotNames_available)

  do.call(methods::new, c(list(Class = "CGNM_result"), CGNM_result_list[present]))
}
