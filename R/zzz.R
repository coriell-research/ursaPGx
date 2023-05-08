star_caller <- NULL

.onLoad <- function(libname, pkgname) {
  path <- system.file("python", package = "ursaPGx")
  star_caller <<- reticulate::import_from_path("star_caller", path = path, delay_load = TRUE)
}
