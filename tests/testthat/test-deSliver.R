test_that("deSliver works", {
  skip_if_not(interactive())
  skip_if_not_installed(c("curl", "googledrive", "httr", "RCurl", "XML"))

  mainDir <- tempdir()
  dPath <- file.path(mainDir, "inputs")
  # set.seed(123)
  # SA <- randomStudyArea(size = 10000000)
  # RTM <- rast(SA, res = 250)

  studyArea <- prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
                          destinationPath = dPath)
  studyArea <- studyArea[studyArea$REGION_NAM == "Lac Seul Upland",]
  #buffer to get some slivers
  targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                         "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  studyArea <- sf::st_transform(studyArea, targetCRS)
  studyArea <- sf::st_buffer(studyArea, 10000)
  frp <- prepInputsFireRegimePolys(destinationPath = dPath, studyArea = studyArea)

  polyAreas <- sort(as.numeric(st_area(frp)), decreasing = TRUE)
  #TODO: correct this bug
  out <- scfmutils::deSliver(frp, threshold = 6.25e5)


  # expect_true(nrow(frp) != nrow(out))
  # expect_true(nrow(out) == 6)
  expect_no_error(scfmutils::deSliver(frp, threshold = 1),
                  message = "deSliver works when nothing is a sliver")
  expect_no_error(scfmutils::deSliver(frp, threshold = polyAreas[1] - 1),
                  message = "deSliver works when nothing is a sliver")
  out <- scfmutils::deSliver(frp, threshold = polyAreas[3]-1)
  expect_true(nrow(out) == 3) #because we used the 3rd largest area

  unlink(mainDir, recursive = TRUE)
})
