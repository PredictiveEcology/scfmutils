utils::globalVariables(c(
  "cellSize"
))

#' @keyword internal
#'
#' @importFrom purrr transpose
#' @importFrom raster extract focal getValues
#' @importFrom stats na.omit
.makeLandscapeAttr <- function(flammableMap, weight, fireRegimePolys) {
  neighMap <- raster::focal(x = flammableMap, w = weight, na.rm = TRUE) # default function is sum(...,na.rm)

  # extract table for each polygon
  valsByPoly <- raster::extract(neighMap, fireRegimePolys, cellnumbers = TRUE)
  valsByPoly <- lapply(valsByPoly, na.omit)
  names(valsByPoly) <- fireRegimePolys$PolyID
  uniqueZoneNames <- unique(fireRegimePolys$PolyID) # get unique zones.
  valsByZone <- lapply(uniqueZoneNames, function(ecoName) {
    aa <- valsByPoly[names(valsByPoly) == ecoName]
    if (is.list(aa)) {
      aa <- do.call(rbind, aa)
    }
    return(aa)
  })
  names(valsByZone) <- uniqueZoneNames

  # Derive frequency tables of number of flammable cells, per polygon type, currently ECOREGION
  nNbrs <- lapply(valsByZone, function(x) {
    nNbrs <- tabulate(x[, 2] + 1, 9) # depends on sfcmLandCoverInit
    names(nNbrs) <- 0:8
    return(nNbrs)
  })

  nFlammable <- lapply(valsByZone, function(x) {
    sum(getValues(flammableMap)[x[, 1]], na.rm = TRUE) # sums flammable pixels in FRI polygons
  })

  landscapeAttr <- purrr::transpose(list(
    cellSize = rep(list(cellSize), length(nFlammable)),
    nFlammable = nFlammable,
    nNbrs = nNbrs,
    cellsByZone = lapply(valsByZone, function(x) x[, 1])
  ))

  landscapeAttr <- lapply(landscapeAttr, function(x) {
    append(x, list(burnyArea = x$cellSize * x$nFlammable))
  })
  names(landscapeAttr) <- names(valsByZone)

  return(landscapeAttr)
}

#' `scfmLandCoverInit`: `genFireMapAttr`
#'
#' @param flammableMap `RasterLayer`. TODO.
#' @param fireRegimePolys TODO
#' @param neighbours TODO
#'
#' @return TODO
#'
#' @export
#' @importFrom raster res
genFireMapAttr <- function(flammableMap, fireRegimePolys, neighbours) {
  # calculate the cell size, total area, and number of flammable cells, etc.
  # All areas in ha
  cellSize <- prod(res(flammableMap)) / 1e4 # in ha

  if (neighbours == 8) {
    w <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3, ncol = 3)
  } else if (neighbours == 4) {
    w <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
  } else {
    stop("illegal neighbours specification")
  }

  landscapeAttr <- .makeLandscapeAttr(flammableMap, w, fireRegimePolys)

  return(invisible(landscapeAttr))
}
