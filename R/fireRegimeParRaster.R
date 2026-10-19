utils::globalVariables(c(
  "flam", "value"
))

#' Map a fire regime polygon attribute onto flammable pixels
#'
#' Joins one column of `fireRegimePolys` (e.g., `pIgnition`, `p0`, `pSpread`) onto the pixels of
#' `fireRegimeRas` by `PolyID`.
#' Used by the `scfmIgnition`, `scfmEscape`, and `scfmSpread` modules.
#'
#' @template fireRegimePolys
#'
#' @param column character, name of the column of `fireRegimePolys` to map.
#'
#' @param fireRegimeRas `SpatRaster` of fire regime polygon IDs (`PolyID`).
#'
#' @template flammableMap
#'
#' @return `SpatRaster` with the geometry of `fireRegimeRas`, holding the polygon's value of
#'   `column` in each pixel, and `NA` where `flammableMap` is `0` (any value other than `1`).
#'   Pixels where `flammableMap` is `NA` keep their polygon's value.
#'
#' @export
#' @importFrom data.table := data.table
#' @importFrom terra rast setValues
fireRegimeParRaster <- function(fireRegimePolys, column, fireRegimeRas, flammableMap) {
  polyAttrs <- as.data.frame(fireRegimePolys)
  if (!column %in% names(polyAttrs)) {
    stop("fireRegimePolys has no column '", column, "'")
  }

  polyValues <- data.table(PolyID = polyAttrs[["PolyID"]], value = polyAttrs[[column]])
  pixValues <- data.table(PolyID = as.vector(fireRegimeRas), flam = as.vector(flammableMap))
  pixValues <- polyValues[pixValues, on = c("PolyID")]
  pixValues[flam != 1, value := NA]

  return(setValues(rast(fireRegimeRas), pixValues$value))
}
