utils::globalVariables(c(
  "iter"
))

#' `scfmDriver`: `genSimLand`
#'
#' Buffers polygon, generates index raster
#'
#' @param coreLand TODO
#' @param buffDist distance to buffer `coreLand`
#' @param flammableMap `SpatRaster` with values 0 indicating non-flammable pixels, 1 flammable.
#'
#' @return list containing `fireRegimePoly`, `landscapeIndex`, `flammableMap` objects.
#'
#' @export
#' @importFrom terra rasterize
#' @importFrom reproducible Cache postProcess
#' @importFrom sf st_buffer st_cast st_difference st_is_valid
genSimLand <- function(coreLand, buffDist, flammableMap = NULL) {
  stopifnot(!is.null(buffDist))

  coreLand$fooField <- 1

  bfireRegimePoly <- st_buffer(coreLand, buffDist)
  bfireRegimePoly$fooField <- 0

  if (!all(st_is_valid(bfireRegimePoly))) {
    bfireRegimePoly <- st_buffer(bfireRegimePoly, width = 0)
  }

  message("creating polyLandscape...")
  #union doe not work, neither does default st_join
  polyLandscape <- rbind(coreLand, bfireRegimePoly)
  polyLandscape <- st_difference(polyLandscape)
  polyLandscape <- st_cast(polyLandscape, "MULTIPOLYGON")

  flammableMap <- Cache(postProcess, flammableMap, studyArea = polyLandscape,
                        useSAcrs = TRUE, filename2 = NULL)

  #Generate landscape Index raster
  landscapeIndex <- rasterize(polyLandscape, flammableMap, fun = 'min', "fooField")

  calibrationLandscape <- list(polyLandscape, landscapeIndex, flammableMap)
  names(calibrationLandscape) <- c("fireRegimePoly", "landscapeIndex", "flammableMap")

  return(calibrationLandscape)
}

#' `scfmDriver`: `makeDesign`
#'
#' @note This version of `makeDesign` is the simplest possible.
#'
#' @param indices TODO
#' @param targetN TODO
#' @param pEscape TODO
#' @param pmin TODO
#' @param pmax TODO
#' @param q TODO
#'
#' @return `data.frame` with columns `igLoc`, `p0`, and `p`.
#'
#' @export
#' @importFrom stats runif
makeDesign <- function(indices, targetN, pEscape = 0.1, pmin, pmax, q = 1) {
  ## TODO: Fix makeDesign to work if polygons have no fires
  sampleSize <- round(targetN / pEscape)
  cellSample <- sample(indices, sampleSize, replace = TRUE)
  pVec <- runif(sampleSize)^q
  pVec <- pVec * (pmax - pmin) + pmin

  ## derive p0 from escapeProb
  ## steal code from scfmRegime and friends.

  p0 <- 1 - (1 - pEscape)^0.125  ## assume 8 neighbours
  ## the preceding approximation seems inadequate in practice.
  #when implemented in scfmDriver, make use of correct derivation of p0 from pEscape based on L
  Temp <- data.frame("igLoc" = cellSample, "p0" = p0, "p" = pVec)

  return(Temp)
}

#' `scfmDriver`: `executeDesign`
#'
#' DESCRIPTION NEEDED
#'
#' @param L TODO
#' @param dT TODO
#' @param maxCells TODO
#'
#' @return TODO
#'
#' @export
#' @importFrom data.table getDTthreads setDTthreads
#' @importFrom terra rast setValues
#' @importFrom reproducible Cache
#' @importFrom SpaDES.tools adj spread2
executeDesign <- function(L, dT, maxCells) {
  ## extract elements of dT into a three column matrix where column 1,2,3 = igLoc, p0, p
  iter <- 0
  probRas <- L
  startTime <- Sys.time()

  .executeDesignInternal <- function(x, L, ProbRas, startTime) { ## L, P are rasters, passed by reference
    iter <<- iter + 1
    currentTime <- Sys.time()
    diffTime <- currentTime - startTime
    units(diffTime) <- "secs"
    timePer <- as.numeric(diffTime) / iter
    timeLeft <- (NROW(dT) - iter) * timePer
    timeLeft <- round(as.difftime(timeLeft, units = "mins")/60, 1)
    nrowDT <- NROW(dT)
    if (iter %% 200 == 0) {
      message("  ", iter, " of ", nrowDT, " total; estimated time remaining: ",
              format(timeLeft, units = "mins"))
    }

    threadsDT <- data.table::getDTthreads()
    data.table::setDTthreads(1)
    on.exit({data.table::setDTthreads(threadsDT)}, add = TRUE)

    i <- x[1]
    p0 <- x[2]
    p <- x[3]

    nbrs <- as.vector(SpaDES.tools::adj(x = L, i, pairs = FALSE, directions = 8))
    ## nbrs < nbrs[which(L[nbrs] == 1)] #or this?
    nbrs <- nbrs[L[nbrs] == 1] #only flammable neighbours please. also, verify NAs excluded.
    ## nbrs is a vector of flammable neighbours.
    nn <- length(nbrs)
    res <- c(nn, 0, 1)
    if (nn == 0)
      return(res) #really defaults
    ## P is still flammableMap.

    ProbVals <- as.vector(ProbRas)
    ProbVals[nbrs] <- p0
    ProbRas <- setValues(ProbRas, ProbVals)
    #Now it is 1, 0, p0, and NA
    spreadState0 <- SpaDES.tools::spread2(landscape = L,
                                          start = i,
                                          iterations = 1,
                                          maxSize = maxCells,
                                          spreadProb = ProbRas,
                                          asRaster = FALSE)

    tmp <- nrow(spreadState0)
    res[2:3] <- c(tmp - 1,tmp)
    if (tmp == 1) { # the fire did not spread
      return(res)
    }

    ProbRas <- setValues(ProbRas, as.vector(L) * p)
    spreadState1 <- SpaDES.tools::spread2(landscape = L,
                                          start = spreadState0,
                                          spreadProb = ProbRas,
                                          asRaster = FALSE,
                                          maxSize = maxCells)
    #calculate return data
    res[3] <- nrow(spreadState1)
    return(res)
  }

  ## [TM] 2019-02-15: Parallelizing isn't efficient here
  res <- Cache(apply, dT, 1, .executeDesignInternal, L, ProbRas = probRas, startTime = startTime)
  res <- data.frame("nNeighbours" = res[1, ], "initSpreadEvents" = res[2, ], "finalSize" = res[3, ])

  #cbind dT and res, then select the columns we need
  x <- cbind(dT, res)
  x <- x[x$finalSize > 1, ]
  return(x)
}

#' `scfmDriver`: `makeAndExecuteDesign`
#'
#' This is a wrapper around `makeDesign` and `executeDesign`.
#'
#' @param ... objects to pass through to `makeDesign` and `executeDesign`.
#'
#' @return output of `executeDesign` (TODO)
#'
#' @export
makeAndExecuteDesign <- function(...) {
  dots <- list(...)
  designTable <- makeDesign(indices = dots$indices, targetN = dots$targetN,
                            pmin = dots$pmin, pmax = dots$pmax,
                            pEscape = dots$pEscape)

  cD <- executeDesign(L = dots$L,
                      dT = designTable,
                      maxCells = dots$maxCells)

  return(cD)
}

#' `scfmDriver`: escape probability
#'
#' ````
#' 1 - (1-p0)**N = pEscape
#' 1 - pEscape = (1-p0)**N
#' (1 - pEscape)**1/N = 1 - p0
#' p0 = 1 - (1 - pEscape)**1/N
#' ```
#'
#' @param pEscape TODO
#' @param n TODO
#'
#' @return TODO
#'
#' @export
#' @rdname pEscape
hatP0 <- function(pEscape, n = 8) {
  1 - (1 - pEscape) ** (1 / n)
}

#' @param p0 TODO
#' @param w TODO
#' @param hatPE TODO
#'
#' @rdname pEscape
escapeProbDelta <- function(p0, w, hatPE) {
  ## TODO: a real clever boots would minimise the abs log odds ratio.
  abs(sum(w*(1 - (1 - p0) ** (0:8))) - hatPE)
}

#' `scfmDriver`: `calibrateFireRegimePolys`
#'
#' Calibrate fire regime polygons ... (TODO)
#'
#' @param polygonType the names of polygons, i.e. `PolyID`
#' @param targetN the number of fires to simulate during calibration
#' @param fireRegimePolys fire regime polygons
#' @param buffDist buffer distance for cells available to be burned outside of each regime polygon
#' @param pJmp default spread probability for degenerate polygons
#' @param pMin minimum spread probability
#' @param pMax maximum allowable spread probability
#' @param flammableMap a packed `SpatRaster` - see [terra::wrap()]
#' @param plotPath file name specifying an output directory to use for producing plots of the scam
#'                 fit for each polygon.
#' @param optimizer the numerical optimization method to use with scam fitting; see `?scam`.
#'
#' @return TODO
#'
#' @export
#' @importFrom grDevices dev.off png
#' @importFrom terra ncell unwrap
#' @importFrom reproducible Cache checkPath
#' @importFrom rlang eval_tidy
#' @importFrom scam scam
#' @importFrom stats as.formula optimise uniroot
#' @importFrom data.table melt.data.table
calibrateFireRegimePolys <- function(polygonType, targetN, fireRegimePolys,
                                     buffDist, pJmp, pMin, pMax, flammableMap = NULL,
                                     plotPath = NULL, optimizer = "bfgs") {
  #must be a file path when run in parallel as SpatRaster can't be serialized
  flammableMap <- terra::unwrap(flammableMap)
  fireRegimePoly <- fireRegimePolys[fireRegimePolys$PolyID == polygonType,]

  frp <- as.data.table(fireRegimePoly)#drop geometry
  frp <- unique(frp[, geometry := NULL])

  maxBurnCells <- as.integer(round(frp$emfs_ha / frp$cellSize)) ## will return NA if emfs is NA
  if (is.na(maxBurnCells)) {
    warning("maxBurnCells cannot be NA... there is a problem with scfmRegime")
    maxBurnCells = 1
  }

  message("generating buffered landscapes...")
  ## this function returns too much data to be worth caching (4 rasters per poly)
  if (is(fireRegimePolys, "quosure")) {
    fireRegimePolys <- eval_tidy(fireRegimePolys)
  }

  message("running genSimLand() ...")
  calibLand <- genSimLand(fireRegimePoly,
                          buffDist = buffDist, flammableMap = flammableMap)

  ## Need a vector of igniteable cells
  ## Item 1 = L, the flammable Map
  ## Item 2 = B (aka the landscape Index) this denotes buffer
  ## Item 3 = igLoc(index of igniteable cells) L[igloc] == 1 &&(B[igLoc]) == 1 (ie within core)
  index <- 1:ncell(calibLand$flammableMap)
  index[calibLand$flammableMap[] != 1 | is.na(calibLand$flammableMap[])] <- NA
  index[calibLand$landscapeIndex[] != 1 | is.na(calibLand$landscapeIndex[])] <- NA
  index <- index[!is.na(index)]
  if (length(index) == 0)
    stop("polygon has no flammable cells!")

  message(paste0("calibrating for polygon ", polygonType, " (Time: ", Sys.time(), ")"))

  ## NOTE: these functions have been wrapped to allow for simpler caching
  message("running makeAndExecuteDesign()...")
  cD <- Cache(makeAndExecuteDesign,
              indices = index,
              targetN = targetN,
              pmin = pMin, pmax = pMax,
              #TODO: change pEscape to use p0 which is calculated afterward (independently),
              #but naievely inside makeDesign (it assumes 8 neighbours)
              pEscape = ifelse(frp$pEscape == 0, 0.1, frp$pEscape),
              L = calibLand$flammableMap,
              maxCells = maxBurnCells,
              userTags = c("scfmDriver", "executeDesign", polygonType),
              omitArgs = c("indices"))

  count <- 0
  kcount <- 50
  scamFormula <- "finalSize ~ s(p, bs = 'micx', k = kcount)"

  message("fitting scam model...")
  stopifnot(optimizer %in% c("bfgs", "efs", "nlm", "nlm.fd", "optim")) ## see ?scam
  calibModel <- try({
    scam::scam(as.formula(scamFormula), data = cD, optimizer = optimizer)
  }, silent = TRUE)
  while (count < 5 & inherits(calibModel, "try-error")) {
    kcount <- kcount + 5
    count <- count + 1
    message("|_ failed! retrying scam fitting (attempt ", count, "/5) ...")
    calibModel <- try(scam::scam(as.formula(scamFormula), data = cD), silent = TRUE)
  }
  if (inherits(calibModel, "try-error")) {
    stop("could not calibrate fire model.")
  } else {
    message("|_ success!")
    plotPath <- checkPath(plotPath, create = TRUE)
    png(file.path(plotPath, sprintf("scfmDriver_scam_plot_Poly%s.png", polygonType)),
        height = 600, width = 800)
    plot(calibModel, main = paste("polygon", polygonType))
    dev.off()
  }
  xBar <- frp$xBar / frp$cellSize

  if (xBar > 0) {
    ## now for the inverse step.
    Res <- try(stats::uniroot(unirootFunction,
                              calibModel, xBar, # "..."
                              interval = c(min(cD$p), max(cD$p)),
                              extendInt = "no",
                              tol = 0.00001
    ), silent = TRUE)
    if (inherits(Res, "try-error")) {
      ## TODO: should pick the closest value (of min and max) if error is value not of opposite sign
      pJmp <- min(cD$p)
      message("the loess model may underestimate the spread probability for polygon ", polygonType)
    } else {
      pJmp <- Res$root
    }
  } else {
    calibModel <- "No Model"
    Res <- "No Uniroot result"
  }
  ## check convergence, and out of bounds errors etc

  nNbrs <- melt.data.table(frp, id.vars = "PolyID", measure.vars = patterns("nNbr"),
                           variable.name = "nNbr", value.name = "count")
  nNbrs <- nNbrs$count
  w <- nNbrs/sum(nNbrs)
  neighbours <- length(nNbrs) - 1 #because 0 is counted

  hatPE <- frp$pEscape
  browser()
  if (hatPE == 0) {
    # no fires in polygon zone escaped
    p0 <- 0
  } else if (hatPE == 1) {
    # all fires in polygon zone escaped
    p0 <- 1
  } else {
    message("running optimise() to determine p0...")
    res <- optimise(escapeProbDelta,
                    interval = c(hatP0(hatPE, neighbours),
                                 hatP0(hatPE, floor(sum(w * 0:neighbours)))),
                    tol = 1e-4,
                    w = w,
                    hatPE = hatPE)
    p0 <- res[["minimum"]]
    ## It is almost obvious that the true minimum must occur within the interval specified in the
    ## call to optimise, but I have not proved it, nor am I certain that the function being minimised is
    ## monotone.
  }
  ## don't forget to scale by number of years, as well, if your timestep is ever != 1yr
  rate <- fireRegimePoly$ignitionRate * frp$cellSize
  ## fireRegimeModel and this module must agree on an annual time step. How to test / enforce?
  pIgnition <- rate
  ## approximate Poisson arrivals as a Bernoulli process at cell level.
  ## for Poisson rate << 1, the expected values are the same, partially accounting
  ## for multiple arrivals within years. Formerly, I used a poorer approximation
  ## where 1-p = P[x==0 | lambda=rate] (Armstrong and Cumming 2003).
  driverResult <- list(
    PolyID = polygonType,
    pSpread = pJmp,
    p0 = p0,
    naiveP0 = hatP0(frp$pEscape, 8),
    pIgnition = pIgnition,
    maxBurnCells = maxBurnCells,
    # calibModel = calibModel,
    uniroot.Res = Res
  )
  return(driverResult)
}

#' `unirootFunction`
#'
#' @param x TODO
#' @param cM TODO
#' @param xBar TODO
#'
#' @return TODO
#'
#' @export
#' @importFrom stats predict
unirootFunction <- function(x, cM, xBar) {
  predict(cM, list("p" = x)) - xBar
}
