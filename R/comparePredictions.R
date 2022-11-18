utils::globalVariables(c(
  "achievedFRI", "achievedIgnitions", "histMeanSize", "modMeanSize", "PolyId",
  "targetFRI", "targetIgnitions"
))

#' Create `data.table` to compare scfm predictions with historical observations
#'
#' @param sim a `simList` object with completed scfm simulations.
#'   Must contain the following objects:
#'   - `burnSummary`
#'   - `fireRegimePoints`
#'   - `landscapeAttr`
#'   - `scfmDriverPars`
#'   - `scfmRegimePars`
#'
#' @return `comparePredictions_summaryDT` returns a `data.table` object;
#'         other functions return `ggplot` objects.
#'
#' @examples
#' \dontrun{
#' ## assumes user has run scfm to produce the simList `mySimOut`
#' df <- comparePredictions_summaryDT(mySimOut)
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
comparePredictions_summaryDT <- function(sim) {
  out <- lapply(names(sim$scfmDriverPars), function(x) {
    regime <- sim$scfmRegimePars[[x]]
    simLength <- times(sim)$end - times(sim)$start + 1
    driver <- sim$scfmDriverPars[[x]]
    landscapeAttr <- sim$landscapeAttr[[x]]
    fireRegimePoints <- sim$fireRegimePoints[sim$fireRegimePoints$PolyID == as.numeric(x), ]

    ## This is a long way of saying 'sum of fires / (flammable landscape * fire epoch)'.
    ## hist_mfs will be NaN if there were no fires larger than one pixel

    pSpread <- driver$pSpread
    pIg <- regime$ignitionRate
    burnSum <- sim$burnSummary[sim$burnSummary$polyID == x, ]
    burnSum$N <- as.numeric(burnSum$N)
    targetIgnitions <- regime$ignitionRate * landscapeAttr$burnyArea
    achievedIgnitions <- nrow(burnSum) / simLength

    burnSum$areaBurned <- as.numeric(burnSum$areaBurned)
    burnSum <- burnSum[burnSum$N > 1]

    targetFRI  <- 1 / regime$empiricalBurnRate
    #the length of simulation over the proportion that actually burned
    achievedFRI <- simLength / (sum(burnSum$areaBurned) / landscapeAttr$burnyArea)

    #starting in year 1 and ending in year 10 means you have 10 years, so
    meanBurn <- ifelse(nrow(burnSum) == 0, yes = 0, no = mean(burnSum$areaBurned))
    pred <- data.frame("PolyId" = x, #Polygon ID
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
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_minimal
#' @rdname comparePredictions
comparePredictions_meanFireSize <- function(dt) {
  ggplot(dt, aes(x = histMeanSize, y = modMeanSize)) +
    geom_point(aes(histMeanSize, modMeanSize)) +
    labs(y = "modeled mean fire size", x = "historical mean fire size") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyId), vjust = "inward", hjust = "inward") +
    geom_abline(slope = 1)
}

#' @export
#' @rdname comparePredictions
comparePredictions_fireReturnInterval <- function(dt) {
  ## remove the infinite FRI caused by no simulated fires
  ## simulated FRI will likely be off where targetFRI is < 1/4 of the simulated time

  if (nrow(dt[targetFRI < c(times$end - times$start) * 4, ]) > 0) {
    ## TODO: confirm wording of this, with comment above
    warning("achievedFRI may be off where targetFRI is less than 4x the simulated time.")
  }

  ## TODO: remove the targetFRI filter below?
  ggplot(dt[!is.infinite(achievedFRI) &
              targetFRI < c(times$end - times$start) * 4], aes(x = targetFRI, y = achievedFRI)) +
    geom_point() +
    labs(y = "simulation FRI (years)", x = "estimated FRI (years)") +
    theme_minimal() +
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyId, vjust = "inward", hjust = "inward"))
}

comparePredictions_annualIgnitions <- function(dt) {
  ggplot(dt, aes(x = targetIgnitions, y = achievedIgnitions)) +
    geom_point() +
    labs(y = "simulation annual ignitions", x = "estimated annual ignitions") +
    theme_minimal() +  #+
    geom_abline(slope = 1) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_text(aes(label = PolyId, vjust = "inward", hjust = "inward"))
}
