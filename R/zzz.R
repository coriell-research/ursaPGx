cyrius <- NULL
.onLoad <- function(libname, pkgname) {
    reticulate::configure_environment(pkgname)
    cyrius <<- import("cyrius", delay_load = TRUE)
}