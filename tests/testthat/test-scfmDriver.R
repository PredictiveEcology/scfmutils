test_that("scfmDriver works with on polygons with no fire", {
  dPath <- withr::local_tempdir("inputs_")

  withr::local_seed(123)

  center <- terra::vect(cbind(-1349980, 6986895))
  terra::crs(center) <- paste(
    "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
    "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  )
  SA <- terra::buffer(center, 20000)
  SA$studyArea <- "Fort McMurrayish"
  SA <- sf::st_as_sf(SA)
  SA$PolyID <- 1
  SA$ignitionRate <- NA
  RTM <- terra::rast(SA, resolution = 250)
  RTM[] <- 1
  RTM <- terra::mask(RTM, SA)

  out <- calibrateFireRegimePolys(
    polygonType = 1,
    targetN,
    fireRegimePolys = SA,
    buffDist = 400,
    pJmp = 0.24,
    pMin = 0.2,
    pMax = 0.28,
    flammableMap = RTM,
    plotPath = NULL,
    outputPath = NULL,
    optimizer = "bfgs"
  )

  expect_true(out$pSpread == 0)
})

## PredictiveEcology/scfm#22: calibrateFireRegimePolys() finds p0 with optimise() over
## [hatP0(pEscape, 8), hatP0(pEscape, floor(mean number of flammable neighbours))].
## When every flammable cell of a polygon has 8 flammable neighbours, the two ends are equal.
##
## The toy landscape is 20 x 20 pixels of 250 m (6.25 ha), all flammable. genFireMapAttr() counts
## only neighbours inside the raster, so cells on the raster edge have 3 or 5.

toyFlammableMap <- function() {
  terra::rast(
    nrows = 20,
    ncols = 20,
    xmin = 0,
    xmax = 5000,
    ymin = 0,
    ymax = 5000,
    crs = "EPSG:3005",
    vals = 1L
  )
}

## one fire regime polygon (PolyID 1) covering the square [lo, hi] x [lo, hi] (m),
## with neighbour counts from genFireMapAttr() and fixed fire regime attributes
toyFireRegimePoly <- function(lo, hi, pEscape = 0.2, neighbours = 8) {
  bbox <- sf::st_bbox(c(xmin = lo, ymin = lo, xmax = hi, ymax = hi), crs = sf::st_crs(3005))
  frp <- sf::st_sf(PolyID = 1L, geometry = sf::st_as_sfc(bbox))
  frp <- genFireMapAttr(toyFlammableMap(), frp, neighbours = neighbours)
  frp$ignitionRate <- 1e-4 ## fires / ha / yr
  frp$pEscape <- pEscape
  frp$xBar <- 25 ## ha
  frp$emfs_ha <- 100 ## ha
  frp
}

toyCalibrate <- function(frp) {
  withr::local_options(list(reproducible.useCache = FALSE, scfmutils.driver.plot.scam = FALSE))
  outPath <- withr::local_tempdir()
  suppressMessages(calibrateFireRegimePolys(
    polygonType = 1L,
    targetN = 200,
    fireRegimePolys = frp,
    buffDist = 1000,
    pJmp = 0.23,
    pMin = 0.185,
    pMax = 0.253,
    flammableMap = terra::wrap(toyFlammableMap()),
    plotPath = outPath,
    outputPath = outPath
  ))
}

## escape probability implied by p0, given the polygon's neighbour counts; optimise() finds p0
## to 1e-4, which puts this within 1e-3 of pEscape
toyEscapeProb <- function(p0, frp) {
  nNbrCols <- grep("^nNbr_", names(frp), value = TRUE)
  nNbrs <- unlist(sf::st_drop_geometry(frp)[nNbrCols])
  k <- as.integer(sub("nNbr_", "", nNbrCols))
  sum(nNbrs / sum(nNbrs) * (1 - (1 - p0)^k))
}

test_that("p0 is calibrated when every flammable cell has 8 flammable neighbours (scfm#22)", {
  withr::local_seed(123)
  frp <- toyFireRegimePoly(1250, 3750) ## central 10 x 10 block: no cell on the raster edge
  expect_equal(frp$nNbr_8, frp$nFlammable)

  out <- toyCalibrate(frp)
  expect_lt(abs(toyEscapeProb(out$p0, frp) - frp$pEscape), 1e-3)
})

test_that("p0 is calibrated for a fully flammable polygon that reaches the raster edge", {
  withr::local_seed(123)
  frp <- toyFireRegimePoly(0, 5000) ## the whole raster
  expect_lt(frp$nNbr_8, frp$nFlammable)

  out <- toyCalibrate(frp)
  expect_lt(abs(toyEscapeProb(out$p0, frp) - frp$pEscape), 1e-3)
})

test_that("p0 is calibrated from 4-neighbour counts", {
  withr::local_seed(123)
  frp <- toyFireRegimePoly(0, 5000, neighbours = 4) ## nNbr_0 to nNbr_4 only
  expect_false("nNbr_8" %in% names(frp))

  out <- toyCalibrate(frp)
  expect_lt(abs(toyEscapeProb(out$p0, frp) - frp$pEscape), 1e-3)
})
