## getFirePoints_NFDB_scfm() used `.../current_version/NFDB_point.zip`, which CFS renamed to
## `NFDB_point_shp.zip`; the old name returns HTTP 404.
test_that("the NFDB point default URL is the full shapefile archive", {
  url <- scfmutils:::nfdbPointUrl()
  expect_match(url, "/current_version/NFDB_point_shp\\.zip$")
  expect_no_match(url, "large_fires")
})

test_that("the NFDB point default URL is reachable", {
  skip_on_cran()
  skip_if_offline("cwfis.cfs.nrcan.gc.ca")
  expect_identical(attr(curlGetHeaders(scfmutils:::nfdbPointUrl(), timeout = 30), "status"), 200L)
})
