utils::globalVariables(c(
  ".SD", "achievedFRI", "achievedIgnitions", "grp", "histMeanSize", "modMeanSize", "N", "PolyID",
  "targetFRI", "targetIgnitions"
))

#' Create `data.table` to compare scfm predictions with historical observations
#'
#' @param burnSummary `data.table`, produced by `scfmSpread` module
#' @param fireRegimePoints `SpatialPointsDataFrame`, produced by `scfmRegime` module
#' @param landscapeAttr list of landscape attributes for each polygon, produced by `scfmLandcoverInit` module
#' @param scfmDriverPars list of burn parameters for each polygon, produced by `scfmDriverPars` module
#' @param scfmRegimePars list of fire regime parameters, produced by `scfmRegime` module
#' @param times list of simulation start and end times (i.e., output from `times(sim)`)
#'
#' @return `comparePredictions_summaryDT` returns a `data.table` object;
#'         other functions return `ggplot` objects.
#'
#' @examples
#' \dontrun{
#' ## assumes user has run scfm to produce the simList `mySimOut`
#' dt <- comparePredictions_summaryDT(scfmDriverPars = mySimOut$scfmDriverPars,
#'                                    scfmRegimePars = mySimOut$scfmRegimePars,
#'                                    landscapeAttr = mySimOut$landscapeAttr,
#'                                    fireRegimePoints = mySimOut$fireRegimePoints,
#'                                    burnSummary = mySimOut$burnSummary)
#'
#' gg_mfs <- comparePredictions_meanFireSize(dt)
#' gg_fri <- comparePredictions_fireReturnInterval(dt)
#' gg_ign <- comparePredictions_annualIgnitions(dt)
#' gg_frp <- plot_fireRegimePolys(mySimOut$fireRegimePolys)
#'
#' grid.extra::grid.arrange(fps,  gg_mfs,  gg_fri,  gg_ign, nrow = 2, ncol = 2)
#' }
#'
#' @author Ian Eddy
#' @export
#' @importFrom data.table rbindlist copy
#' @importFrom SpaDES.core times
#' @rdname comparePredictions
comparePredictions_summaryDT <- function(scfmDriverPars = NULL,
                                         scfmRegimePars = NULL,
                                         landscapeAttr = NULL,
                                         fireRegimePoints = NULL,
                                         burnSummary = NULL,
                                         times = NULL) {
  if (any(is.null(scfmDriverPars), is.null(scfmRegimePars), is.null(landscapeAttr),
          is.null(fireRegimePoints), is.null(burnSummary), is.null(times))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  out <- lapply(names(scfmDriverPars), function(x) {
    regime <- scfmRegimePars[[x]]
    simLength <- times$end - times$start + 1
    driver <- scfmDriverPars[[x]]
    landscapeAttr <- landscapeAttr[[x]]
    fireRegimePoints <- fireRegimePoints[fireRegimePoints$PolyID == as.numeric(x), ]

    ## This is a long way of saying 'sum of fires / (flammable landscape * fire epoch)'.
    ## hist_mfs will be NaN if there were no fires larger than one pixel

    #median fire size is not used by scfm but is worth recording
    #regimes where mean is much greater than median will be hard to recreate

    medianFireSize <- median(fireRegimePoints$SIZE_HA) #should be no need for na.rm

    pSpread <- driver$pSpread
    pIg <- regime$ignitionRate

    if (!"grp" %in% names(burnSummary)) {
      stop("burnSummary data.table does not have a 'grp' column.\n",
           "Are you running a recent version of scfmSpread (>= 2.0.0)?")
    }

    ## burnSummary data.table from scfmSpread is coded as follows:
    ## grp 1: pixels from fires ignited in SAR & spread in SAR
    ## grp 2: pixels from fires ignited in SAR & spread outside SAR
    ## grp 3: pixels from fires ignited outside SAR & spread in SAR
    ## grp 4: pixels from fires ignited outside SAR & spread outside SAR
    burnSum <- burnSummary[PolyID == x, ]
    targetIgnitions <- pIg * landscapeAttr$burnyArea
    achievedIgnitions <- nrow(burnSum[grp %in% 1, ]) / simLength ## incl grp 2 would double count ignitions

    #escapes
    targetEscapes <- regime$pEscape * targetIgnitions
    achievedEscapes <- nrow(burnSum[grp %in% 1 & N > 1])/ simLength

    ## mean fire size: mean size of all fires ignited and escaped in SAR, regardless of where spread
    burnSum1 <- burnSum[grp %in% c(1, 2), lapply(.SD, sum), by = c("igLoc", "year"), .SDcols = "areaBurned"]
    burnSum1 <- burnSum1[areaBurned > landscapeAttr$cellSize, ]
    meanFireSize <- ifelse(nrow(burnSum1) == 0, 0, mean(burnSum1$areaBurned))

    ## Mean Annual Area Burned: total area of all burned pixels in SAR over n years of simulation
    burnSum2 <- burnSum[grp %in% c(1, 3), lapply(.SD, sum), by = c("igLoc", "year"), .SDcols = "areaBurned"]
    burnSum2 <- burnSum1[areaBurned > landscapeAttr$cellSize, ]
    MAAB <- sum(burnSum2$areaBurned) / simLength

    achievedFRI <- simLength / ( sum(0, burnSum2$areaBurned) / landscapeAttr$burnyArea)
    targetFRI  <- 1 / regime$empiricalBurnRate

    pred <- data.frame("PolyID" = x,
                       "histMeanSize" = regime$xBar, ## predicted (empirical) mean size of fires
                       "histMedianSize" = medianFireSize,
                       "modMeanSize" = meanFireSize,
                       "achievedFRI" = achievedFRI,
                       "targetFRI" = targetFRI,
                       "burnableArea_ha" = landscapeAttr$burnyArea,
                       "targetIgnitions" = targetIgnitions,
                       "achievedIgnitions" = achievedIgnitions,
                       "targetEscapes" = targetEscapes,
                       "achievedEscapes" = achievedEscapes,
                       "pEscape" = regime$pEscape, #escape prob (no. fires > cellSize / no. fires)
                       "pSpread" = pSpread, ## spread probability estimated from the SCAM model
                       "pIgnition" = pIg) ## ignition probability of a single pixel
    return(pred)
  })
  return(rbindlist(out))
}

#' @param dt scfm summary `data.table` produced by `comparePredictions_summaryDT()`
#'
#' @export
#' @importFrom ggplot2 aes geom_abline geom_point geom_text ggplot labs
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw
#' @rdname comparePredictions
comparePredictions_meanFireSize <- function(dt) {
  if (any(is.null(dt))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  ggplot(dt, aes(x = histMeanSize, y = modMeanSize)) +
    geom_point(aes(histMeanSize, modMeanSize)) +
    labs(x = "historical mean fire size", y = "modeled mean fire size") +
    theme_bw() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID), vjust = "inward", hjust = "inward") +
    geom_abline(slope = 1)
}

#' @export
#' @rdname comparePredictions
comparePredictions_fireReturnInterval <- function(dt, times) {
  if (any(is.null(dt), is.null(times))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  ## remove the infinite FRI caused by no simulated fires
  ## simulated FRI will likely be off where targetFRI is < 1/4 of the simulated time

  if (nrow(dt[targetFRI < c(times$end - times$start) * 4, ]) > 0) {
    ## TODO: confirm wording of this, with comment above
    warning("achievedFRI may be off where targetFRI is less than 4x the simulated time.")
  }

  ## TODO: remove the targetFRI filter below. plot those points differently to indicate poor estimates
  ggplot(dt[!is.infinite(achievedFRI) & targetFRI < c(times$end - times$start) * 4],
         aes(x = targetFRI, y = achievedFRI)) +
    geom_point() +
    labs(y = "simulation FRI (years)", x = "estimated FRI (years)") +
    theme_bw() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID, vjust = "inward", hjust = "inward"))
}

#' @export
#' @rdname comparePredictions
comparePredictions_annualIgnitions <- function(dt) {
  if (any(is.null(dt))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  dt <- copy(dt) #avoid adding per ha cols (or add them?)
  dt[, targetIgnitions_Mha := targetIgnitions/burnableArea_ha * 1e6]
  dt[, achiedIgnitions_Mha := targetIgnitions/burnableArea_ha * 1e6]


  ggplot(dt, aes(x = targetIgnitions_Mha, y = achievedIgnitions_Mha)) +
    geom_point() +
    labs(y = "simulation annual ignitions (per Mha)", x = "estimated annual ignitions (per Mha)") +
    theme_bw() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID, vjust = "inward", hjust = "inward"))
}

#' @export
#' @rdname comparePredictions
comparePredictions_annualEscapes <- function(dt) {
  if (any(is.null(dt))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  dt <- copy(dt) #avoid adding per ha cols (or add them?)
  dt[, targetEscapes_Mha := targetEscapes/burnableArea_ha * 1e6]
  dt[, achiedEscapes_Mha := targetEscapes/burnableArea_ha * 1e6]


  ggplot(dt, aes(x = targetEscapes_Mha, y = achievedEscapes_Mha)) +
    geom_point() +
    labs(y = "simulation annual escapes (per Mha)", x = "estimated annual escapes (per Mha)") +
    theme_bw() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID, vjust = "inward", hjust = "inward"))
}

#' @export
#' @rdname comparePredictions
comparePredictions_fireDistribution <- function(dt){
  if (any(is.null(dt))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  ggplot(dt, aes(x = histMedianSize, y = histMeanSize)) +
    geom_point() +
    labs(y = "estimated mean fire size (ha)", x = "estimated median fire size (ha)") +
    theme_bw() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID, vjust = "inward", hjust = "inward"))
}
