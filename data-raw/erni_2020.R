## code to extract Erni et al.'s fire regime tables for FRUs (A1) and FRTs (A2)

library(pdftools)
library(stringr)

raw_text <- pdf_text("~/Downloads/erni_2020/erni_2020.pdf") |>
  str_split("\\n", simplify = TRUE)

# Fire Regime Unit (FRU) attributes -----------------------------------------------------------

## Table A1 (all of p12) has the fire regime parameters for FRUs
A1_start <- raw_text[12, ] |>
  str_which("Attributes and ﬁre metrics of the 60 Fire Regime Units \\(FRUs\\)")

A1_end <- raw_text[12, ] |>
  str_which("Note: Spring and summer columns provide the distribution")

A1_raw <- raw_text[12, (A1_start + 1):(A1_end - 1)]

A1_tbl <- A1_raw |>
  grep(pattern = "Downloaded from cdnsciencepub", x = _, invert = TRUE, value = TRUE) |>
  grep(pattern = "For personal use only", x = _, invert = TRUE, value = TRUE) |>
  str_subset(".+") |>
  str_replace_all("\\s{2,}", "|") |>
  str_replace_all("^\\|", "")

A1_head <- c("", str_split(A1_tbl[1], "\\|", simplify = TRUE)) |>
  paste(str_split(A1_tbl[2], "\\|", simplify = TRUE)) |>
  str_replace_all("\\(.*\\)", "") |>
  trimws()

fru_attr <- paste(A1_head, collapse = "|") |>
  c(A1_tbl[-(1:2)]) |>
  textConnection() |>
  read.csv(sep = "|", col.names = A1_head) |>
  as.data.frame()

usethis::use_data(fru_attr, overwrite = TRUE)

# Fire Regime Type (FRT) attributes -----------------------------------------------------------

## Table A2 (top of p13) has the fire regime parameters for FRTs
A2_start <- raw_text[13, ] |>
  str_which("Attributes and ﬁre metrics of the 15 Fire Regime Types \\(FRTs\\)")

A2_end <- raw_text[13, ] |>
  str_which("Note: Spring and summer columns provide the distribution")

A2_raw <- raw_text[13, (A2_start + 1):(A2_end - 1)]

A2_tbl <- A2_raw |>
  grep(pattern = "Downloaded from cdnsciencepub", x = _, invert = TRUE, value = TRUE) |>
  grep(pattern = "For personal use only", x = _, invert = TRUE, value = TRUE) |>
  str_subset(".+") |>
  str_replace_all("\\s{2,}", "|") |>
  str_replace_all("^\\|", "")

A2_head <- c("", str_split(A2_tbl[1], "\\|", simplify = TRUE)) |>
  paste(str_split(A2_tbl[2], "\\|", simplify = TRUE)) |>
  str_replace_all("\\(.*\\)", "") |>
  trimws()

frt_attr <- paste(A2_head, collapse = "|") |>
  c(A2_tbl[-(1:2)]) |>
  textConnection() |>
  read.csv(sep = "|", col.names = A2_head) |>
  as.data.frame()

usethis::use_data(frt_attr, overwrite = TRUE)
