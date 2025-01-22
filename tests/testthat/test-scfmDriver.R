test_that("scfmDriver works with on polygons with no fire", {
  dPath <- withr::local_tempdir("inputs_")

  withr::local_seed(123)

  center <- terra::vect(cbind(-1349980, 6986895))
  terra::crs(center) <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
                              "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  SA <- terra::buffer(center, 20000)
  SA$studyArea <- "Fort McMurrayish"
  SA <- sf::st_as_sf(SA)
  SA$PolyID <- 1
  SA$ignitionRate <- NA
  RTM <- terra::rast(SA, resolution = 250)
  RTM[] <- 1
  RTM <- terra::mask(RTM, SA)


  out <- calibrateFireRegimePolys(
    polygonType = 1, targetN, fireRegimePolys = SA,
    buffDist = 400, pJmp = 0.24, pMin = 0.2, pMax = 0.28, flammableMap = RTM,
    plotPath = NULL, outputPath = NULL, optimizer = "bfgs"
  )

  expect_true(out$pSpread == 0)
})
