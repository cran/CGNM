test_that("every exported function is actually defined in the namespace", {
  ns <- asNamespace("CGNM")
  exported <- getNamespaceExports(ns)

  missing <- exported[!vapply(exported, exists, logical(1), envir = ns, inherits = FALSE)]

  expect_length(missing, 0)
})
