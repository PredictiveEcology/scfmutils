utils::globalVariables(c(
  "PolyID"
))

#' Plot fire regime polygons
#'
#' @param fireRegimePolys `sf` polygon or multipolygon object defining the fire regime polygons
#'
#' @return a `ggplot` object
#'
#' @export
#' @importFrom ggplot2 aes geom_sf ggplot scale_fill_discrete theme_minimal
plot_fireRegimePolys <- function(fireRegimePolys) {
  if (!is.factor(fireRegimePolys$PolyID)) {
    fireRegimePolys$PolyID <- as.factor(fireRegimePolys$PolyID)
  }

  ggplot(fireRegimePolys) +
    geom_sf(aes(fill = PolyID)) +
    scale_fill_discrete() +
    theme_minimal()
}
