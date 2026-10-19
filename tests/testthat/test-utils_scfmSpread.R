## pixels 1 to 10 are inside the reporting area (SAR).
## Fire A is ignited at pixel 1 (inside) and burns pixels 1, 2 (inside) and 11 (outside);
## fire B is ignited at pixel 12 (outside) and burns pixels 12 (outside) and 3 (inside).
burnDT <- function() {
  data.table::data.table(initialPixels = c(1L, 1L, 1L, 12L, 12L), pixels = c(1L, 2L, 11L, 12L, 3L))
}

test_that("countBurnedPixelsInSAR() counts each fire's pixels by where it was ignited and burned", {
  out <- countBurnedPixelsInSAR(burnDT(), pixInSAR = 1:10)
  data.table::setorderv(out, c("initialPixels", "grp"))

  expect_identical(out$initialPixels, c(1L, 1L, 12L, 12L))
  expect_identical(out$grp, 1:4)
  expect_identical(out$N, c(2L, 1L, 1L, 1L))
})

test_that("countBurnedPixelsInSAR() does not modify burnDT", {
  dt <- burnDT()
  countBurnedPixelsInSAR(dt, pixInSAR = 1:10)

  expect_identical(dt, burnDT())
})
