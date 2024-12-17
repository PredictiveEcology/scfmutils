test_that("download of fireRegimePolys works", {
  skip_if_not(interactive())
  skip_if_not_installed(c("curl", "googledrive", "httr", "RCurl", "XML"))

  mainDir <- tempdir()
  dPath <- file.path(mainDir, "inputs")
  set.seed(123)
  center <- terra::vect(cbind(-1349980, 6986895))
  terra::crs(center) <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
                       "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  SA <- terra::buffer(center, 20000)
  SA$studyArea <- "Fort McMurrayish"
  RTM <- terra::rast(SA, res = 250)
  RTM[] <- 1
  RTM <- terra::mask(RTM, SA)

  frp1 <- prepInputsFireRegimePolys(destinationPath = dPath, studyArea = SA,
                                    rasterToMatch = RTM)
  frp2 <- prepInputsFireRegimePolys(destinationPath = dPath, studyArea = SA,
                                    rasterToMatch = RTM, type = "FRT")
  frp3 <- prepInputsFireRegimePolys(destinationPath = dPath, studyArea = SA,
                                    rasterToMatch = RTM, type = "FRU")
  #need BC study area for BEC
  # frp4 <- prepInputsFireRegimePolys(destinationPath = dPath, studyArea = SA, rasterToMatch = RTM,
  #                                   type = "BEC")
  #
  ## get all available species for 2001

  expect_true(inherits(frp1, "sf"))
  expect_true(inherits(frp2, "sf"))
  expect_true(inherits(frp3, "sf"))

  unlink(mainDir, recursive = TRUE)
})
