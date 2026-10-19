## 2 x 2 landscape with fire regime polygons 1 and 2; pixel 2 is not flammable,
## pixel 4 is outside every polygon.
parInputs <- function() {
  frr <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 1, 2, NA))
  names(frr) <- "PolyID"
  frp <- terra::as.polygons(frr)
  frp$pSpread <- c(0.2, 0.3)[match(frp$PolyID, c(1, 2))]
  flam <- terra::rast(frr, vals = c(1, 0, 1, NA))

  list(fireRegimePolys = frp, fireRegimeRas = frr, flammableMap = flam)
}

test_that("fireRegimeParRaster() maps a polygon column onto flammable pixels", {
  inputs <- parInputs()

  out <- fireRegimeParRaster(
    fireRegimePolys = inputs$fireRegimePolys,
    column = "pSpread",
    fireRegimeRas = inputs$fireRegimeRas,
    flammableMap = inputs$flammableMap
  )

  expect_true(terra::compareGeom(out, inputs$fireRegimeRas))
  expect_equal(as.vector(out), c(0.2, NA, 0.3, NA))
})

test_that("fireRegimeParRaster() gives the same result for sf and SpatVector polygons", {
  inputs <- parInputs()
  args <- list(
    column = "pSpread",
    fireRegimeRas = inputs$fireRegimeRas,
    flammableMap = inputs$flammableMap
  )

  outVect <- do.call(fireRegimeParRaster, c(list(fireRegimePolys = inputs$fireRegimePolys), args))
  outSf <- do.call(
    fireRegimeParRaster,
    c(list(fireRegimePolys = sf::st_as_sf(inputs$fireRegimePolys)), args)
  )

  expect_equal(as.vector(outSf), as.vector(outVect))
})

test_that("fireRegimeParRaster() errors on a missing column", {
  inputs <- parInputs()

  expect_error(
    fireRegimeParRaster(
      fireRegimePolys = inputs$fireRegimePolys,
      column = "pIgnition",
      fireRegimeRas = inputs$fireRegimeRas,
      flammableMap = inputs$flammableMap
    ),
    "no column 'pIgnition'"
  )
})
