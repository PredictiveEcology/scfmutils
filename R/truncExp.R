## TODO: don't use `T` as a variable -- `T` is equivalent to `TRUE` in R

#' Various tools for a (right) Truncated Exponential Distribution
#'
#' \eqn{Z ~ TE(\theta, T)}: pdf is \eqn{\frac{\theta \exp(-z \theta)}{1 - \exp(-T\theta)}}
#'
#' where `T` is the truncation point or upper bound and \eqn{\theta} is the shape parameter
#' in this application, `x` are fire sizes `>= shift`, which is a lower bound
#' and `z = log(x / shift)` are the scaled log transformed sizes which seem to
#' fit a truncated exponential distribution fairly well.
#'
#' Originally written by Steve in 1999 in support of Cumming CJFR 2001.
#' Has been in use by BEACONs and was acquired from Pierre Vernier in May 17 2014.
#'
#' - `pTE()` is the distribution function;
#' - `dTE()` is the density function;
#' - `qTE()` is the quantile function;
#' - `rTE()` is the random generation function.
#'   In fire size applications `exp(rTE(n, theta, T)) * shift` will generate `n` random fire sizes;
#' - `qqTE()` produces a quantile-quantile plot of vector `x` against a `TE(theta, T)`;
#' - `ETE()` is TODO;
#' - `ETEx()` is TODO;
#' - `ExpBar()` is TODO;
#' - `HannonDayiha()` implements the estimator of Hannon and Dayiha (1999),
#'   ported from 1999 C language implementation by SGC June 2004.
#'
#' @references
#' Patrick M. Hannona & Ram C. Dahiyaa (1999) Estimation of parameters for the truncated
#' exponential distribution. Communications in Statistics - Theory and Methods 28(11): 2591-2612.
#' \doi{10.1080/03610929908832440}.
#'
#'
#' @param x TODO
#' @param theta TODO
#' @param T TODO
#'
#' @export
#' @rdname truncExp
pTE <- function(x, theta, T) {
  const <- 1 - exp(-T * theta)
  num <- 1 - exp(-x * theta)
  return(ifelse(x >= 0 & x <= T, num / const, NaN))
}

#' @export
#' @rdname truncExp
dTE <- function(x, theta, T) {
  const <- theta / (1 - exp(-T * theta))
  return(ifelse(x >= 0 & x <= T, const * exp(-x * theta), NaN))
}

#' @param p TODO
#'
#' @export
#' @rdname truncExp
qTE <- function(p, theta, T) {
  const <- 1 - exp(-T * theta)
  q <- p * const
  q <- log1p(-q)
  q <- q / -theta
  return(ifelse(p >= 0 & p <= 1, q, NaN))
}

#' @param n TODO
#'
#' @export
#' @importFrom stats runif
#' @rdname truncExp
rTE <- function(n, theta, T) {
  p <- runif(n)
  return(qTE(p, theta, T))
}

#' @param plot.it logical indicating whether to plot the resulting figure
#' @param xlab,ylab label text for x- and y-axes
#' @param ... additional arguments passed to `plot()`
#'
#' @export
#' @importFrom graphics abline plot
#' @importFrom stats ppoints
#' @rdname truncExp
qqTE <- function(x, theta = stop("missing theta"), T = stop("missing T"), plot.it = TRUE,
                 xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), ...) {
  sx <- sort(x)
  qs <- ppoints(sx)
  sy <- qTE(qs, theta, T)
  ux <- max(c(sx, sy))
  lx <- -0.05 * ux
  if (plot.it) {
    plot(sx, sy, xlab = xlab, ylab = ylab, xlim = c(lx, ux), ylim = c(lx, ux), ...)
    abline(0, 1) # may not want this always?
  }
  invisible(list(x = sx, y = sy))
}

#' @param a TODO
#'
#' @export
#' @rdname truncExp
ETE <- function(theta, T, a = 0) {
  # if Z ~ TE(theta,T) return E[Z|z>a]
  shape <- (T - a) * theta
  part1 <- exp(shape) - 1
  part2 <- 1 - shape * (part1^-1)
  return(a + part2 / theta)
}

#' @param shift TODO
#'
#' @export
#' @rdname truncExp
ETEx <- function(theta, T, a = 0, shift = 200) {
  # z = log(x/shift), x>=shift
  # x = exp(z)*shift
  # Z ~ TE(theta,T)
  # return E[X|z>a]
  eta <- 1 - theta
  n1 <- exp(T)^eta
  n2 <- exp(a)^eta
  num <- theta * exp(a * theta) * (n1 - n2)
  den <- (1 - exp(-(T - a) * theta)) * eta
  x <- num / den
  return(x * shift)
}

#' @param Z TODO
#'
#' @export
#' @rdname truncExp
ExpBar <- function(Z) {
  # max likelihood estimate of the shape parameter for
  # an unshifted exponential distribution
  return(sum(Z) / length(Z))
}

#' @param Tspec TODO
#'
#' @export
#' @rdname truncExp
HannonDayiha <- function(Z, Tspec = 0) {
  # Z is a vector of (log transformed and shifted) sizes
  # if Tspec nonzero, assume T=Tspec and just estimate theta
  # if a>0 use shifted subset Z[Z>=a]-a

  MAXSTEPS <- 100
  epsilon <- 1e-8

  mle_step <- function(theta, zbar, T) {
    x <- T / expm1(theta * T)
    x <- 1 / theta - x - zbar
    return(x)
  }

  estimate <- function(T0, theta, flag = 0) {
    if (flag == 0) {
      That <- Zmax
    } else {
      if (flag == 1) {
        x <- expm1(Zmax * theta) / length(Z)
        an <- log1p(x) / theta
        That <- Zmax + an
      } else {
        if (flag == 2) {
          That <- T0
        } else {
          stop(paste("Illegal: flag = ", flag))
        }
      }
    }
    lower <- zbar
    while (mle_step(lower, zbar, That) <= 0) {
      lower <- lower * 0.9
    }
    upper <- zbar
    while ((tmp <- mle_step(upper, zbar, That)) >= 0) {
      upper <- upper * 1.1
    }
    theta <- (upper + lower) / 2
    if (theta <= 0 || is.nan(theta)) {
      stop("theta estimate negative or NaN")
    }

    done <- FALSE
    Steps <- 0
    while (Steps < MAXSTEPS && abs(tmp <- mle_step(theta, zbar, That)) > epsilon) {
      if (tmp < 0) { # theta too big
        upper <- theta
        theta <- (theta + lower) / 2
      } else {
        lower <- theta
        theta <- (upper + theta) / 2
      }
      Steps <- Steps + 1
    }

    if (Steps == MAXSTEPS) {
      warning("Estimate: convergence failure")
      That <- -1
    }
    return(list(theta = theta, That = That))
  }

  sigvar <- function(theta, T, n) {
    xx <- exp(-1 * theta * T)
    x <- 1 - xx
    x <- x^-2.0
    x <- T^2 * xx * x
    x <- (theta^-2 - x)^-1
    x <- x / n
    x <- sqrt(x)
    return(x)
  }

  if (any(Z <= 0) || any(is.na(Z))) {
    stop("Negatives or NaNs")
  }

  nZ <- length(Z)
  Zmax <- max(Z)
  zbar <- ExpBar(Z)

  if (Tspec == 0) {
    step1 <- estimate(0, 0, 0)
    T0 <- step1$That
    theta <- step1$theta
    step2 <- estimate(T0, theta, 1)
    That <- step2$That
    theta <- step2$theta
  } else {
    step1 <- estimate(Tspec, 0, 2)
    That <- step1$That
    theta <- step1$theta
  }

  Tsd <- expm1(theta * That) / theta
  Tsd <- Tsd / nZ
  thsd <- sigvar(theta, That, nZ)
  return(list(theta = theta, thsd = thsd, That = That, Tsd = Tsd))
}
