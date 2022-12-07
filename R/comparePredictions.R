utils::globalVariables(c(
  "achievedFRI", "achievedIgnitions", "histMeanSize", "modMeanSize", "PolyID",
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
#' @importFrom data.table rbindlist
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

    pSpread <- driver$pSpread
    pIg <- regime$ignitionRate
    burnSum <- burnSummary[burnSummary$polyID == x, ]
    burnSum$N <- as.numeric(burnSum$N)
    targetIgnitions <- pIg * landscapeAttr$burnyArea
    achievedIgnitions <- nrow(burnSum) / simLength

    burnSum$areaBurned <- as.numeric(burnSum$areaBurned)
    burnSum <- burnSum[burnSum$N > 1]

    targetFRI  <- 1 / regime$empiricalBurnRate
    #the length of simulation over the proportion that actually burned
    achievedFRI <- simLength / (sum(burnSum$areaBurned) / landscapeAttr$burnyArea)

    #starting in year 1 and ending in year 10 means you have 10 years, so
    meanBurn <- ifelse(nrow(burnSum) == 0, yes = 0, no = mean(burnSum$areaBurned))
    pred <- data.frame("PolyID" = x,
                       "histMeanSize" = regime$xBar, ## predicted (empirical) mean size of fires
                       "modMeanSize" = meanBurn, ## modeled mean size of fires
                       "achievedFRI" = achievedFRI,
                       "targetFRI" = targetFRI,
                       "targetIgnitions" = targetIgnitions,
                       "achievedIgnitions" = achievedIgnitions,
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
    labs(y = "modeled mean fire size", x = "historical mean fire size") +
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

  ggplot(dt, aes(x = targetIgnitions, y = achievedIgnitions)) +
    geom_point() +
    labs(y = "simulation annual ignitions", x = "estimated annual ignitions") +
    theme_bw() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyID, vjust = "inward", hjust = "inward"))
}
