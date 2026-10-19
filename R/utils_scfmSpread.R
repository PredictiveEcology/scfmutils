utils::globalVariables(c(
  ".N", "grp", "initialPixels", "pixels"
))

#' `scfmSpread`: count burned pixels inside and outside the reporting area
#'
#' Codes each burned pixel by whether its fire was ignited inside the reporting area (SAR,
#' `studyAreaReporting`) and whether the pixel itself is inside it, then counts the pixels of
#' each fire in each group.
#'
#' @details
#' The groups, as read by `comparePredictions_summaryDT()`:
#'
#' - `grp` 1: pixels from fires ignited in SAR & spread in SAR;
#' - `grp` 2: pixels from fires ignited in SAR & spread outside SAR;
#' - `grp` 3: pixels from fires ignited outside SAR & spread in SAR;
#' - `grp` 4: pixels from fires ignited outside SAR & spread outside SAR.
#'
#' @param burnDT `data.table` of burned pixels, as returned by [SpaDES.tools::spread2()] with
#'   `asRaster = FALSE`, with columns `initialPixels` (the pixel where each fire was ignited)
#'   and `pixels`.
#'
#' @param pixInSAR integer vector of the indices of the pixels inside the reporting area.
#'
#' @return `data.table` with columns `initialPixels`, `grp`, and `N` (number of pixels).
#'
#' @export
#' @importFrom data.table := copy
countBurnedPixelsInSAR <- function(burnDT, pixInSAR) {
  tempDT <- copy(burnDT)

  tempDT[(initialPixels %in% pixInSAR) & (pixels %in% pixInSAR), grp := 1L]
  tempDT[(initialPixels %in% pixInSAR) & !(pixels %in% pixInSAR), grp := 2L]
  tempDT[!(initialPixels %in% pixInSAR) & (pixels %in% pixInSAR), grp := 3L]
  tempDT[!(initialPixels %in% pixInSAR) & !(pixels %in% pixInSAR), grp := 4L]
  tempDT <- tempDT[, .N, by = c("initialPixels", "grp")]

  return(tempDT)
}
