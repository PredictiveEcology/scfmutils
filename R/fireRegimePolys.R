utils::globalVariables(c(
  "geometry", "USETHIS"
))

#' @export
#' @rdname prepInputsFireRegimePolys
fireRegimePolyTypes <- function() {
  c("BECNDT", "BECSUBZONE", "BECZONE", "ECODISTRICT", "ECOREGION", "ECOPROVINCE", "ECOZONE",
    "FRT", "FRU")
}

#' `prepInputsFireRegimePolys`
#'
#' Create fire regime polygons for `scfmRegime`.
#'
#' @param url character. URL from which to download and prepare fire regime polygons.
#'            Defaults are provided for Canadian ecodistrict, ecoregion, ecoprovince, and ecozone,
#'            as well as national Fire Regime Types and Fire Regime Units from Erni et al. (2020)
#'            \doi{10.1139/cjfr-2019-0191}.
#'
#' @param destinationPath character. Path to directory where data will be downloaded.
#'
#' @param studyArea `sf` object corresponding to the study area of interest.
#'
#' @param rasterToMatch TODO
#'
#' @param type character. The polygon type to use:
#'             Must be one of "ECODISTRICT", "ECOREGION" (default), "ECOPROVINCE", "ECOZONE",
#'             "FRT", or "FRU".
#'             If `url` to BEC shapefile is provided, can also be one of:
#'             "BECNDT", "BECSUBZONE", or "BECZONE".
#'
#' @export
#' @importFrom dplyr %>% group_by summarise ungroup
#' @importFrom raster crs
#' @importFrom reproducible prepInputs
#' @importFrom sf st_as_sf st_collection_extract st_union
#'
#' @examples
#' library(sf)
#' library(sp)
#'
#' ## random study area in central Alberta
#' studyArea <- SpatialPoints(data.frame(lon = -115, lat = 55), proj4string = CRS("EPSG:4326")) |>
#' st_as_sf() |>
#'   st_transform(paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
#'                    "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")) |>
#'   as_Spatial() |>
#'   SpaDES.tools::randomStudyArea(center = _, seed = 60, size = 1e10) |>
#'   st_as_sf()
#'
#' frpEcoregion <- prepInputsFireRegimePolys(studyArea = studyArea, type = "ECOREGION")
#' plot(frpEcoregion)
#'
#' frpFRT <- prepInputsFireRegimePolys(studyArea = studyArea, type = "FRT")
#' plot(frpFRT)
#'
#' frpFRU <- prepInputsFireRegimePolys(studyArea = studyArea, type = "FRU")
#' plot(frpFRU)
prepInputsFireRegimePolys <- function(url = NULL, destinationPath = tempdir(),
                                      studyArea = NULL, rasterToMatch = NULL, type = "ECOREGION") {
  type <- toupper(type)
  allowedTypes <- fireRegimePolyTypes()

  stopifnot(type %in% allowedTypes)

  if (is.null(url)) {
    if (grepl("BEC", type)) {
      ## no public url available? user must pass their own, e.g. google drive link
      stop("url must be provided when using type 'BECNDT', 'BECSUBZONE' or 'BECZONE'.")
    } else {
      urlList <- list(
        ecodistrict = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
        ecoregion = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
        ecoprovince = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
        ecozone = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
        frt = "https://zenodo.org/record/4458156/files/FRT.zip",
        fru = "https://zenodo.org/record/4458156/files/FRU.zip"
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
                    overwrite = TRUE) ## TODO: doesn't reproject -- fix upstream?

  ## workaround issues with prepInputs() not reprojecting:
  if (!is.null(rasterToMatch)) {
    tmp <- sf::st_transform(tmp, raster::crs(rasterToMatch))
  } else if (is.null(rasterToMatch) && !is.null(studyArea)) {
    tmp <- sf::st_transform(tmp, sf::st_crs(studyArea))
  }

  if (grepl("^ECO", type)) {
    cols2keep <- substr(type, 1, 10) ## colname abbrev to 10 chars
  } else if (grepl("^BEC.*ZONE", type)) {
    cols2keep <- c("ZONE", "SUBZONE")
  } else if (type == "BECNDT") {
    cols2keep <- names(tmp)[names(tmp) %in% c("NTRL_DSTRD", "NTRLDSTRBN")]
  } else if (type == "FRT") {
    cols2keep <- "Cluster"
  } else if (type == "FRU") {
    cols2keep <- "GRIDCODE"
  }

  tmp <- tmp[, cols2keep]

  tmp$USETHIS <- switch(type,
                        BECNDT = tmp[[cols2keep[1]]],
                        BECSUBZONE = paste0(tmp[["ZONE"]], "_", tmp[["SUBZONE"]]),
                        BECZONE = tmp[["ZONE"]],
                        tmp[[cols2keep[1]]])
  tmp$USETHIS <- as.factor(tmp$USETHIS)

  tmp2 <- group_by(tmp, USETHIS) %>% summarise(geometry = sf::st_union(geometry)) %>% ungroup()
  polys <- sf::st_collection_extract(tmp2)
  polys[["PolyID"]] <- as.integer(1:nrow(polys))
  polys[["USETHIS"]] <- NULL

  return(polys)
}

#' Check for various issues with `fireRegimePolys`
#'
#' @param fireRegimePolys TODO
#' @param studyArea TODO
#' @param rasterToMatch TODO
#' @param flammableMap TODO
#' @param sliverThresh TODO
#' @param cacheTag TODO
#'
#' @return a cleaned up `fireRegimePolys` object
#'
#' @export
#' @importFrom raster compareCRS
#' @importFrom reproducible Cache
#' @importFrom sf st_area st_is_longlat
checkForIssues <- function(fireRegimePolys, studyArea, rasterToMatch, flammableMap, sliverThresh, cacheTag) {
  compareCRS(rasterToMatch, flammableMap, fireRegimePolys) ## TODO: is there a better check?

  if (is.null(fireRegimePolys[["PolyID"]])) {
    stop("please supply fireRegimePolys with a PolyID")
  }

  if (st_is_longlat(fireRegimePolys)) {
    stop("scfm requires projected coordinate systems - lat/lon too prone to error.")
  }
  fireRegimePolys$trueArea <- round(st_area(fireRegimePolys), digits = 0)

  if (any(as.numeric(fireRegimePolys$trueArea) < sliverThresh)) {
    message("sliver polygon(s) detected. Merging to their nearest valid neighbour")
    fireRegimePolys <- Cache(deSliver, fireRegimePolys, threshold = sliverThresh, userTags = cacheTag)
  }

  return(fireRegimePolys)
}

#' Merge sliver polygons into non-sliver neighbours
#'
#' The threshold is applied to the area of the multipolygon object, not each individual polygon.
#' Non-sliver polygons keep their original attributes.
#' Intended to be used when it is important to retain the original extent of an
#' area while removing sliver polygons.
#'
#' @param x an `sf` POLYGONS or MULTIPLOYGONS object
#'
#' @param threshold the minimum area below which a polygon is considered a sliver
#'
#' @return an object of class `sf` with sliver polygons merged to their nearest valid neighbour.
#'
#' @export
#' @importFrom sf st_area st_buffer st_cast st_is_valid st_nearest_feature st_union
deSliver <- function(x, threshold) {
  x$tempArea <- as.numeric(st_area(x))

  ## determine slivers by area
  xSlivers <- x[x$tempArea < threshold, ]
  xNotSlivers <- x[x$tempArea >= threshold, ]
  if (nrow(xNotSlivers) < 1) {
    stop("Threshold exceeds the area of every polygon. Please select a smaller number")
  }

  ## split slivers from multipolygon, or nearest feature may be incorrect
  xSlivers <- suppressWarnings(st_cast(xSlivers, "POLYGON"))

  ## find nearest non-sliver
  nearestFeature <- st_nearest_feature(xSlivers, xNotSlivers)

  ## merge each sliver polygon into nearest neighbour
  mergeSlivers <- lapply(
    unique(nearestFeature),
    FUN = function(i,
                   ns = xNotSlivers,
                   s = xSlivers,
                   nf = nearestFeature) {
      featurePolys <- nearestFeature == i
      xMerge <- st_union(s[featurePolys, ])
      yMerge <- ns[i, ]
      out <- st_union(x = xMerge, y = yMerge) # convert slivers back to multipolygon
      yMerge$geometry <- out # update the geometry
      return(yMerge)
    }
  )
  otherPolys <- xNotSlivers[!(1:nrow(xNotSlivers) %in% nearestFeature),]
  if (length(mergeSlivers) > 1) {
    ## these polygons must be tracked and merged.
    ## they may be nrow(0) if every feature was modified in some way
    if (nrow(otherPolys) != 0) {
      mergeSlivers <- do.call(rbind, mergeSlivers)
      m <- rbind(otherPolys, mergeSlivers)
    } else {
      m <- rbind(mergeSlivers)
    }
  } else { #lapply over length 1 is special
    if (nrow(otherPolys) != 0) {
      mergeSlivers <- mergeSlivers[[1]]
      m <- rbind(mergeSlivers, otherPolys)
    } else {
      m <- mergeSlivers[[1]]
    }
  }
  ## the geometry will be sfc
  m <- st_cast(m, to = "MULTIPOLYGON")
  m$tempArea <- NULL # remove the temporary column

  ## remove self-intersecting geometries
  if (any(!st_is_valid(m))) {
    #m <- rgeos::gBuffer(spgeom = m, byid = TRUE, width = 0) ## TODO: use sf here
    m <- st_buffer(m, dist = 0) ## done by geometry
  }

  return(m)
}
