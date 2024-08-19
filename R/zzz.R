.onLoad <- function(libname, pkgname) {
  ## set options using the approach used by devtools
  opts <- options()
  opts.scfmutils <- list( # nolint
    scfmutils.driver.plot.scam = TRUE
  )
  toset <- !(names(opts.scfmutils) %in% names(opts))
  if (any(toset)) options(opts.scfmutils[toset])

  invisible()
}
