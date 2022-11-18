utils::globalVariables(c(
  "geometry", "USETHIS"
))

#' `prepInputsFireRegimePolys`
#'
#' Create fire regime polygons for `scfmRegime`.
#'
#' @param url character. URL from which to download and prepare fire regime polygons.
#'            Defaults are provided for Canadian ecodistrict, ecoregion, ecoprovince, and ecozone.
#'
#' @param destinationPath character. Path to directory where data will be downloaded.
#'
#' @param studyArea `sf` object corresponding to the study area of interest.
#'
#' @param rasterToMatch TODO
#'
#' @param type character. The polygon type to use:
#'             Must be one of "ECODISTRICT", "ECOREGION" (default), "ECOPROVINCE", or "ECOZONE".
#'             If `url` to BEC zones shapefile is provided, can also be one of "BECSUBZONE" or
#'             "BECZONE".
#'
#' @export
#' @importFrom dplyr %>% group_by summarise ungroup
#' @importFrom reproducible prepInputs
#' @importFrom sf st_as_sf st_collection_extract st_union
prepInputsFireRegimePolys <- function(url = NULL, destinationPath = tempdir(),
                                      studyArea = NULL, rasterToMatch = NULL, type = "ECOREGION") {
  type <- toupper(type)
  allowedTypes <- c("BECSUBZONE", "BECZONE", "ECODISTRICT", "ECOREGION", "ECOPROVINCE", "ECOZONE")

  stopifnot(type %in% allowedTypes)

  if (is.null(url)) {
    if (grepl("BEC", type)) {
      ## no public url available? user must pass their own, e.g. google drive link
      stop("url must be provided when using type 'BECSUBZONE' or 'BECZONE'.")
    } else {
      urlList <- list(
        ecodistrict = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
        ecoregion = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
        ecoprovince = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
        ecozone = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
      )
      url <- urlList[[tolower(type)]]
    }
  }

  if (!is.null(studyArea) && is(studyArea, "SpatialPolygons")) {
    studyArea <- sf::st_as_sf(studyArea)
  }

  tmp <- prepInputs(url = url,
                    destinationPath = destinationPath,
                    studyArea = studyArea,
                    rasterToMatch = rasterToMatch,
                    fun = "sf::st_read",
                    overwrite = TRUE)

  if (grepl("^ECO", type)) {
    cols2keep <- substr(type, 1, 10) ## colname abbrev to 10 chars
  } else {
    cols2keep <- c("ZONE", "SUBZONE")
  }

  tmp <- tmp[, cols2keep]

  tmp$USETHIS <- switch(type,
                        BECSUBZONE = paste0(tmp[["ZONE"]], "_", tmp[["SUBZONE"]]),
                        BECZONE = tmp[["ZONE"]],
                        ECODISTRICT = tmp[[cols2keep[1]]],
                        ECOREGION = tmp[[cols2keep[1]]],
                        ECOPROVINCE = tmp[[cols2keep[1]]],
                        ECOZONE = tmp[[cols2keep[1]]])
  tmp$USETHIS <- as.factor(tmp$USETHIS)

  tmp2 <- group_by(tmp, USETHIS) %>% summarise(geometry = sf::st_union(geometry)) %>% ungroup()
  polys <- sf::st_collection_extract(tmp2)
  polys[["PolyID"]] <- 1:nrow(polys)

  return(polys)
}
