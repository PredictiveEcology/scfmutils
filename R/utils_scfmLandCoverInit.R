utils::globalVariables(c(
  ".", ".N", "cell", "cellSize", "flam", "focal_sum"
))

#' `scfmLandCoverInit`: `.makeLandscapeAttr`
#'
#' Calculate the cell size, total area, and number of flammable cells, etc.
#' All areas in `ha`.
#'
#' @keywords internal
#'
#' @importFrom data.table := as.data.table data.table dcast
#' @importFrom dplyr left_join
#' @importFrom purrr transpose
#' @importFrom terra extract focal values res
#' @importFrom stats na.omit
.makeLandscapeAttr <- function(flammableMap, weight, fireRegimePolys, neighbours) {

  cellSize <- prod(res(flammableMap)) / 1e4 # in ha
  neighMap <- focal(x = flammableMap, w = weight, na.rm = TRUE) # default function is sum(..., na.rm)

  # extract table for each polygon
  valsByPoly <- extract(neighMap, fireRegimePolys, cells = TRUE, ID = TRUE) ## TODO: use terra
  valsByPoly <- as.data.table(valsByPoly)
  valsByPoly[, flam := values(flammableMap, mat = FALSE)[cell]]
  valsByPoly <- valsByPoly[flam == 1]

  #get the FRP ID
  tempDT <- data.table(PolyID = fireRegimePolys$PolyID, ID = 1:nrow(fireRegimePolys))
  valsByPoly <- valsByPoly[tempDT, on = c("ID")]

  valsByZone <- lapply(fireRegimePolys$PolyID, FUN = function(x, df = valsByPoly){
    df[PolyID == x]
  })

  #there are occasional NAs - rasterize/extract differences?
  valsByZone <- lapply(valsByZone, na.omit)
  names(valsByZone) <- fireRegimePolys$PolyID

  #Derive frequency tables of number of flammable cells, per polygon type, currently ECOREGION
  nNbrs <- lapply(valsByZone, function(x) {
    nNbrs <- x[, .N, .(focal_sum, PolyID)] # depends on sfcmLandCoverInit

    possibleNbrs <- data.table(nbr = 0:neighbours)
    nNbrs <- nNbrs[possibleNbrs, on = c("focal_sum" = "nbr")]
    nNbrs[, focal_sum := paste0("nNbr_", focal_sum)]
    nNbrs[is.na(N), N := 0]
    nNbrs <- dcast(nNbrs, formula = PolyID ~ focal_sum, value.var = "N")

    return(nNbrs)
  })

  nNbrs <- rbindlist(nNbrs)
  #find total flammable pixels in cell
  flamByPoly <- valsByPoly[, .(flam = sum(flam, na.rm = TRUE)), PolyID]

  #assign nFlammable, cellSize, burnyArea, and nBrs to fireRegimePolys
  fireRegimePolys$nFlammable <- flamByPoly$flam
  fireRegimePolys$cellSize <- cellSize
  fireRegimePolys$burnyArea <- cellSize * fireRegimePolys$nFlammable

  fireRegimePolys <- left_join(fireRegimePolys, nNbrs, by = "PolyID")

  return(fireRegimePolys)
}

#' `scfmLandCoverInit`: `genFireMapAttr`
#'
#' @param flammableMap `SpatRaster`. TODO.
#' @param fireRegimePolys TODO
#' @param neighbours TODO
#'
#' @return TODO
#'
#' @export
#' @importFrom terra res
genFireMapAttr <- function(flammableMap, fireRegimePolys, neighbours) {
  if (neighbours == 8) {
    w <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3, ncol = 3)
  } else if (neighbours == 4) {
    w <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
  } else {
    stop("illegal neighbours specification")
  }

  landscapeAttr <- .makeLandscapeAttr(flammableMap, w, fireRegimePolys, neighbours = neighbours)

  return(invisible(landscapeAttr))
}
