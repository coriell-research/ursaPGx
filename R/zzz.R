star_caller <- NULL

.onLoad <- function(libname, pkgname) {
  # Declare Cyrius' Python dependencies (see inst/python/requirements.txt).
  # reticulate resolves these with uv into a cached ephemeral environment the
  # first time Python is used, unless the user has already selected a Python
  # (e.g. via use_condaenv(), use_python(), or RETICULATE_PYTHON).
  reticulate::py_require(c(
    "numpy>=1.16",
    "scipy>=1.2",
    "pysam>=0.15.3",
    "statsmodels>=0.9"
  ))

  path <- system.file("python", package = "ursaPGx")
  star_caller <<- reticulate::import_from_path("star_caller", path = path, delay_load = TRUE)
}
