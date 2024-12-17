withr::local_options(list(
  warnPartialMatchArgs = TRUE,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE
))

library(testthat)
library(scfmutils)

test_check("scfmutils")
