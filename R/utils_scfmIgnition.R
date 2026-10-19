#' `scfmIgnition`: draw this year's ignitions
#'
#' Each pixel ignites independently with probability `pIg`.
#'
#' @param pIg ignition probability: either a `SpatRaster` of per-pixel probabilities
#'   (pixels with `NA` never ignite), or a single number applied to every non-`NA` pixel
#'   of `fireRegimeRas`.
#'
#' @param fireRegimeRas `SpatRaster` of fire regime polygon IDs (`PolyID`).
#'
#' @return integer vector of the indices of the pixels that ignite, in increasing order.
#'
#' @export
#' @importFrom stats runif
#' @importFrom terra ncell
calcIgnitions <- function(pIg, fireRegimeRas) {
  if (inherits(pIg, "SpatRaster")) {
    stopifnot(ncell(fireRegimeRas) == ncell(pIg))
    igs <- which(runif(ncell(fireRegimeRas)) < as.vector(pIg)) ## NA are allowed and omitted
  } else {
    fireRegimePixels <- which(!is.na(as.vector(fireRegimeRas)))
    igs <- fireRegimePixels[which(runif(length(fireRegimePixels)) <= pIg)]
  }

  return(igs)
}
