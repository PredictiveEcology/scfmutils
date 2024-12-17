test_that("download of fireRegimePolys works", {
  skip_if_not(interactive())
  skip_if_not_installed(c("curl", "googledrive", "httr", "RCurl", "XML"))

  mainDir <- tempdir()
  dPath <- file.path(mainDir, "inputs")
  set.seed(123)
  SA <- LandR::randomStudyArea(size = 10000000)
  RTM <- rast(SA, res = 250)

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
