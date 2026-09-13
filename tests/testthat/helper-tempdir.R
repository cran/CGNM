# Runs expr with the working directory switched to a fresh temp directory,
# restoring the original working directory and removing the temp directory
# afterwards regardless of success or failure. Used for tests of functions
# that read/write CGNM_log folders relative to the working directory.
with_temp_wd <- function(expr) {
  old_wd <- getwd()
  tmp <- tempfile("cgnm-test-")
  dir.create(tmp)
  on.exit({
    setwd(old_wd)
    unlink(tmp, recursive = TRUE)
  })
  setwd(tmp)
  force(expr)
}
