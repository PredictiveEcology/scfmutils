## comparePredictions_summaryDT(): one fire regime polygon of 2500 ha, one simulated year.
## scfmSpread codes burnSummary rows by where the fire was ignited and where its pixels burned, relative
## to the reporting area (SAR): grp 1 in -> in, grp 2 in -> out, grp 3 out -> in, grp 4 out -> out.

summaryInputs <- function() {
  frp <- data.frame(PolyID = 1L, cellSize = 6.25, burnyArea = 2500, ignitionRate = 4e-4, pEscape = 0.5,
                    xBar = 30, p0 = 0.1, pSpread = 0.23, empiricalBurnRate = 0.01)
  pts <- data.frame(PolyID = 1L, SIZE_HA = c(1, 20, 40))
  bs <- data.table::data.table(
    igLoc = c(1L, 2L, 3L), year = 1, PolyID = 1L, N = c(8, 5, 16),
    grp = c(1L, 2L, 3L),            # in -> in; in -> out; out -> in
    areaBurned = c(50, 30, 100)
  )
  list(fireRegimePoints = pts, burnSummary = bs, fireRegimePolys = frp, times = list(start = 1, end = 1))
}

test_that("the fire return interval counts the area burned inside the reporting area (groups 1 and 3)", {
  dt <- do.call(comparePredictions_summaryDT, summaryInputs())
  ## inside the area: 50 ha (grp 1) + 100 ha (grp 3); the 30 ha of grp 2 burned outside it
  expect_equal(as.numeric(dt$achievedFRI), 1 / ((50 + 100) / 2500))
})

test_that("the mean fire size counts fires ignited inside the reporting area (groups 1 and 2)", {
  dt <- do.call(comparePredictions_summaryDT, summaryInputs())
  expect_equal(dt$modMeanSize, mean(c(50, 30)))
})
