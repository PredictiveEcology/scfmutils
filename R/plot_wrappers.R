utils::globalVariables(c(
  "PolyID"
))

#' Plot fire regime polygons
#'
#' @template fireRegimePolys
#'
#' @param title character, the plot title
#'
#' @returns a `ggplot` object
#'
#' @export
#' @importFrom ggplot2 aes geom_sf ggplot scale_fill_discrete theme_minimal
plot_fireRegimePolys <- function(fireRegimePolys, title = NULL) {
  if (!is.factor(fireRegimePolys$PolyID)) {
    fireRegimePolys$PolyID <- as.factor(fireRegimePolys$PolyID)
  }
  gg_frp <- ggplot(fireRegimePolys) +
    geom_sf(aes(fill = PolyID)) +
    scale_fill_discrete() + ## TODO: use same palette as plot_fireRegimeRas ??
    theme_bw()
  
  if (!is.null(title)) {
    gg_frp <- gg_frp + ggtitle(title)
  }
  return(gg_frp)
}

#' Plot fire regime raster
#'
#' @param x `SpatRaster` object corresponding to a fire regime raster
#'
#' @param title character, the plot title
#'
#' @returns `ggplot` object
#'
#' @export
#' @importFrom ggplot2 ggplot ggtitle scale_fill_brewer theme_bw
#' @importFrom tidyterra geom_spatraster
plot_fireRegimeRas <- function(x, title) {
  gg_frr <- ggplot() +
    geom_spatraster(data = terra::as.factor(x)) +
    scale_fill_brewer(palette = "Paired", type = "qual", na.value = "transparent") +
    theme_bw()
  
  if (!is.null(title)) {
    gg_frr <- gg_frr + ggtitle(title)
  }
  return(gg_frr)
}

#' Plot age map
#'
#' @param x `SpatRaster` object corresponding to stand age or time since disturbance map
#'
#' @param title character, the plot title
#'
#' @param maxAge the maximum age to plot
#'
#' @returns `ggplot` object
#'
#' @export
#' @importFrom ggplot2 ggplot ggtitle scale_fill_distiller theme_bw
#' @importFrom tidyterra geom_spatraster
plot_ageMap <- function(x, title, maxAge) {
  x[x > maxAge] <- maxAge
  
  gg_am <- ggplot() +
    geom_spatraster(data = x) +
    scale_fill_distiller(palette = "Greens", direction = 1, na.value = "transparent") +
    theme_bw()
  
  if (!is.null(title)) {
    gg_am <- gg_am + ggtitle(title)
  }
  return(gg_am)
}

#' Plot burn maps
#'
#' @param x `SpatRaster` object corresponding to a current or cumulative burn map.
#'
#' @param title character, the plot title
#'
#' @returns `ggplot` object
#'
#' @export
#' @importFrom ggplot2 ggplot ggtitle theme_bw
#' @importFrom tidyterra geom_spatraster
#' @importFrom viridis scale_fill_viridis
plot_burnMap <- function(x, title) {
  gg_bm <- ggplot() +
    geom_spatraster(data = x) +
    scale_fill_viridis(na.value = "transparent") +
    theme_bw()
  
  if (!is.null(title)) {
    gg_bm <- gg_bm + ggtitle(title)
  }
  return(gg_bm) 
}

#' Plot flammable map
#'
#' @param x `SpatRaster` object corresponding to a flammability map.
#'
#' @param title character, the plot title
#'
#' @returns `ggplot` object
#'
#' @export
#' @importFrom ggplot2 ggplot ggtitle scale_fill_distiller theme_bw
#' @importFrom tidyterra geom_spatraster
plot_flammableMap <- function(x, title) {
  gg_fm <- ggplot() +
    geom_spatraster(data = x) +
    scale_fill_distiller(palette = "RdBu", na.value = "transparent") +
    theme_bw()
  
  if (!is.null(title)) {
    gg_fm <- gg_fm + ggtitle(title)
  }
  return(gg_fm)
}