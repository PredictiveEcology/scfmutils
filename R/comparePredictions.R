utils::globalVariables(c(
  ".SD", "achievedEscapes", "achievedEscapes_Mha", "achievedFRI", "achievedIgnitions",
  "achievedIgnitions_Mha", "areaBurned", "burnableArea_ha", "pSpread", "burnyArea",
  "grp", "histMeanSize", "histMedianSize", "modMeanSize", "N", "PolyID",
  "targetEscapes", "targetEscapes_Mha", "targetFRI", "targetIgnitions", "targetIgnitions_Mha"
))

#' Create `data.table` to compare scfm predictions with historical observations
#'
#' @param burnSummary `data.table`, produced by `scfmSpread` module
#' @param fireRegimePoints `sf` object produced `scfmRegime` module
#' @param fireRegimePolys `sf` object modified by `scfm` modules
#' @param times list of simulation start and end times (i.e., output from `times(sim)`)
#'
#' @return `comparePredictions_summaryDT` returns a `data.table` object;
#'         other functions return `ggplot` objects.
#'
#' @examples
#' \dontrun{
#' ## assumes user has run scfm to produce the simList `mySimOut`
#' dt <- comparePredictions_summaryDT(fireRegimePoints = mySimOut$fireRegimePoints,
#'                                    burnSummary = mySimOut$burnSummary,
#'                                    fireRegimePolys = mySimOut$fireRegimePolys,
#'                                    times = times(mySimOut))
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
#' @importFrom stats median
#' @rdname comparePredictions
comparePredictions_summaryDT <- function(fireRegimePoints = NULL,
                                         burnSummary = NULL,
                                         fireRegimePolys = NULL,
                                         times = NULL) {

  if (any(is.null(fireRegimePolys$pSpread), is.null(fireRegimePolys$xBar),
          is.null(fireRegimePolys$burnyArea), is.null(fireRegimePoints),
          is.null(burnSummary), is.null(times))) {
    stop("fireRegimePolys is missing columns or insufficient args provided")
  }

  fireIDs <- unique(fireRegimePolys$PolyID)
  out <- lapply(fireIDs, function(x) {
    fireRegimePoly <- fireRegimePolys[fireRegimePolys$PolyID == x,]

    simLength <- times$end - times$start + 1
    fireRegimePoints <- fireRegimePoints[fireRegimePoints$PolyID == as.numeric(x), ]

    ## This is a long way of saying 'sum of fires / (flammable landscape * fire epoch)'.
    ## hist_mfs will be NaN if there were no fires larger than one pixel

    ## median fire size is not used by scfm but is worth recording
    ## regimes where mean is much greater than median will be hard to recreate
    escaped <- fireRegimePoints[fireRegimePoints$SIZE_HA > fireRegimePoly$cellSize,]
    medianFireSize <- median(escaped$SIZE_HA) #should be no need for na.rm

    pSpread <- fireRegimePoly$pSpread
    pIg <- fireRegimePoly$ignitionRate

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
    targetIgnitions <- pIg * fireRegimePoly$burnyArea
    achievedIgnitions <- nrow(burnSum[grp %in% 1, ]) / simLength ## incl grp 2 would double count ignitions

    #escapes
    targetEscapes <- fireRegimePoly$pEscape * targetIgnitions
    achievedEscapes <- nrow(burnSum[grp %in% 1 & N > 1]) / simLength

    ## mean fire size: mean size of all fires ignited and escaped in SAR, regardless of where spread
    burnSum1 <- burnSum[grp %in% c(1, 2), lapply(.SD, sum), by = c("igLoc", "year"), .SDcols = "areaBurned"]
    burnSum1 <- burnSum1[areaBurned > fireRegimePoly$cellSize, ]
    meanFireSize <- ifelse(nrow(burnSum1) == 0, 0, mean(burnSum1$areaBurned))

    ## Mean Annual Area Burned: total area of all burned pixels in SAR over n years of simulation
    burnSum2 <- burnSum[grp %in% c(1, 3), lapply(.SD, sum), by = c("igLoc", "year"), .SDcols = "areaBurned"]
    burnSum2 <- burnSum1[areaBurned > fireRegimePoly$cellSize, ]
    MAAB <- sum(burnSum2$areaBurned) / simLength

    achievedFRI <- simLength / ( sum(0, burnSum2$areaBurned) / fireRegimePoly$burnyArea)
    targetFRI  <- 1 / fireRegimePoly$empiricalBurnRate

    pred <- data.frame("PolyID" = x,
                       "histMeanSize" = fireRegimePoly$xBar, ## predicted (empirical) mean size of fires
                       "histMedianSize" = medianFireSize,
                       "modMeanSize" = meanFireSize,
                       "achievedFRI" = achievedFRI,
                       "targetFRI" = targetFRI,
                       "burnableArea_ha" = fireRegimePoly$burnyArea,
                       "targetIgnitions" = targetIgnitions,
                       "achievedIgnitions" = achievedIgnitions,
                       "targetEscapes" = targetEscapes,
                       "achievedEscapes" = achievedEscapes,
                       "pEscape" = fireRegimePoly$pEscape, ## escape prob (no. fires > cellSize / no. fires)
                       "p0" = fireRegimePoly$p0,  ## p0 and pEscape may indicate something incorrect
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
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw xlab ylab
#' @rdname comparePredictions
comparePredictions_meanFireSize <- function(dt) {
  if (any(is.null(dt))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  ggplot(dt, aes(x = histMeanSize, y = modMeanSize)) +
    geom_point(aes(histMeanSize, modMeanSize)) +
    labs(x = "historical mean fire size (ha)", y = "modeled mean fire size (ha)") +
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
  dt[, achievedIgnitions_Mha := achievedIgnitions/burnableArea_ha * 1e6]


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
  dt[, achievedEscapes_Mha := achievedEscapes/burnableArea_ha * 1e6]


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
comparePredictions_fireDistribution <- function(dt) {
  if (any(is.null(dt))) {
    stop("all arguments must be provided and cannot be NULL.")
  }

  ggplot(dt, aes(x = histMedianSize, y = histMeanSize)) +
    geom_point() +
    xlab("estimated median fire size (ha)") +
    ylab("estimated mean fire size (ha)") +
    theme_bw() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID, vjust = "inward", hjust = "inward"))
}
