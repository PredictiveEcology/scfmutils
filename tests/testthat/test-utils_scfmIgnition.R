## 2 x 2 landscape; the last pixel is outside every fire regime polygon
fireRegimeRas <- function() {
  terra::rast(nrows = 2, ncols = 2, vals = c(1, 1, 2, NA))
}

test_that("calcIgnitions() with a single pIg ignites only pixels inside a fire regime polygon", {
  withr::local_seed(42)
  frr <- fireRegimeRas()

  expect_identical(calcIgnitions(pIg = 1, fireRegimeRas = frr), 1:3)
  expect_identical(calcIgnitions(pIg = 0, fireRegimeRas = frr), integer(0))
})

test_that("calcIgnitions() with a pIg raster ignites pixels by their own probability", {
  withr::local_seed(42)
  frr <- fireRegimeRas()
  pIg <- terra::rast(frr, vals = c(1, 0, NA, 1))

  ## pixel 3 has no pIg; pixel 4 has a pIg although it is outside every polygon
  expect_identical(calcIgnitions(pIg = pIg, fireRegimeRas = frr), c(1L, 4L))
})

test_that("calcIgnitions() rejects a pIg raster of a different size", {
  pIg <- terra::rast(nrows = 3, ncols = 3, vals = 1)

  expect_error(calcIgnitions(pIg = pIg, fireRegimeRas = fireRegimeRas()))
})
