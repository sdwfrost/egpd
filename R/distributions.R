## Distribution functions for EGPD, Discrete EGPD, and Zero-Inflated Discrete EGPD
## Adapted from degpd_functions.R in the degpd-and-zidegpd prototype

## Self-contained GPD CDF/PDF/quantile (replaces extraDistr dependency)

.pgpd_std <- function(q, sigma, xi) {
  ## GPD CDF with location=0, scale=sigma, shape=xi
  xi <- sign(xi) * pmax(abs(xi), 1e-8)
  z <- q / sigma
  p <- 1 - pmax(1 + xi * z, 0)^(-1/xi)
  p[q < 0] <- 0
  p
}

.dgpd_std <- function(x, sigma, xi, log = FALSE) {
  ## GPD density with location=0, scale=sigma, shape=xi
  xi <- sign(xi) * pmax(abs(xi), 1e-8)
  z <- x / sigma
  t <- 1 + xi * z
  t <- pmax(t, 0)
  if (log) {
    out <- -log(sigma) - (1/xi + 1) * log(t)
    out[x < 0] <- -Inf
    out[t <= 0] <- -Inf
  } else {
    out <- (1/sigma) * t^(-1/xi - 1)
    out[x < 0] <- 0
    out[t <= 0] <- 0
  }
  out
}

.qgpd_std <- function(p, sigma, xi) {
  ## GPD quantile with location=0, scale=sigma, shape=xi
  xi <- sign(xi) * pmax(abs(xi), 1e-8)
  sigma * ((1 - p)^(-xi) - 1) / xi
}

## Transformation functions G, g, G^{-1} for EGPD types 1-6

#' Transformation CDF for EGPD
#'
#' @param u values in (0,1)
#' @param type integer 1-6 specifying the transformation type
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#'
#' @return Transformed probabilities
#' @export
p.G <- function(u, type = 1, prob, kappa, delta) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1, 2, 3, 5, 6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (type == 0) {
    return(u)
  } else if (type == 1) {
    return(u^kappa)
  } else if (type == 2) {
    F.min <- 1 - pnorm(1, mean = 0, sd = 1/sqrt(kappa))
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
    p <- (pnorm(u, mean = 1, sd = 1/sqrt(kappa)) - F.min) / (F.max - F.min)
    return(p)
  } else if (type == 3) {
    lower <- 1/32
    upper <- 1/2
    stopifnot(lower <= upper, kappa > 0)
    aa <- rep(lower, length(u))
    bb <- rep(upper, length(u))
    normalize.factor <- pbeta(bb, kappa, kappa) - pbeta(aa, kappa, kappa)
    tt <- pbeta((upper - lower) * u + lower, kappa, kappa)
    tt <- tt / normalize.factor
    return(tt)
  } else if (type == 4) {
    return(1 - pbeta((1 - u)^delta, 1/delta, 2))
  } else if (type == 5) {
    return((1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2))
  } else if (type == 6) {
    return(prob * u^kappa + (1 - prob) * u^delta)
  }
}

#' Transformation density for EGPD
#'
#' @param u values in (0,1)
#' @param type integer 1-6 specifying the transformation type
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#' @param log logical: return log-density?
#'
#' @return Transformation density values
#' @export
d.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA, log = FALSE) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1, 2, 3, 5, 6) && identical(kappa, NA)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && identical(delta, NA)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && identical(prob, NA)) {
    stop("Argument `prob' missing.")
  }
  if (log == FALSE) {
    if (type == 0) {
      return(1)
    } else if (type == 1) {
      return(kappa * u^(kappa - 1))
    } else if (type == 2) {
      F.min <- 1 - pnorm(1, mean = 0, sd = 1/sqrt(kappa))
      F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
      d <- dnorm(u, mean = 1, sd = 1/sqrt(kappa))
      den <- (sqrt(kappa) * d) / (F.max - F.min)
      return(den)
    } else if (type == 3) {
      lower <- 1/32
      upper <- 1/2
      stopifnot(lower <= upper, kappa > 0)
      tt <- rep(0, length(u))
      normalize.factor <- pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa)
      tt[u >= lower & u <= upper] <- dbeta(u[u >= lower & u <= upper],
                                           kappa, kappa) / normalize.factor
      return(tt)
    } else if (type == 4) {
      return(dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 - u)^(delta - 1))
    } else if (type == 5) {
      return((kappa/2) * (1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2 - 1) * dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 - u)^(delta - 1))
    } else if (type == 6) {
      return(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1))
    }
  } else {
    if (type == 0) {
      return(0)
    } else if (type == 1) {
      return(log(kappa) + (kappa - 1) * log(u))
    } else if (type == 2) {
      F.min <- 1 - pnorm(1, mean = 0, sd = 1/sqrt(kappa))
      F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
      d <- dnorm(u, mean = 1, sd = 1/sqrt(kappa), log = TRUE)
      den <- log(sqrt(kappa)) + d - log(F.max - F.min)
      return(den)
    } else if (type == 3) {
      lower <- 1/32
      upper <- 1/2
      stopifnot(lower <= upper, kappa > 0)
      tt <- rep(-Inf, length(u))
      normalize.factor <- log(pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa))
      idx <- u >= lower & u <= upper
      tt[idx] <- dbeta(u[idx], kappa, kappa, log = TRUE) - normalize.factor
      return(tt)
    } else if (type == 4) {
      return(dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) + log(delta) + (delta - 1) * log(1 - u))
    } else if (type == 5) {
      return(log(kappa/2) + (kappa/2 - 1) * log(1 - pbeta((1 - u)^delta, 1/delta, 2)) + dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) +
               log(delta) + (delta - 1) * log(1 - u))
    } else if (type == 6) {
      return(log(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1)))
    }
  }
}

#' Inverse transformation (quantile) for EGPD
#'
#' @param u values in (0,1)
#' @param type integer 1-6
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#'
#' @return Inverse-transformed values
#' @export
q.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1, 2, 3, 5, 6) && identical(kappa, NA)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && identical(delta, NA)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && identical(prob, NA)) {
    stop("Argument `prob' missing.")
  }
  if (type == 0) {
    return(u)
  } else if (type == 1) {
    return(u^(1/kappa))
  } else if (type == 2) {
    F.min <- 1 - pnorm(1, mean = 0, sd = 1/sqrt(kappa))
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
    q <- (qnorm(u * (F.max - F.min) + F.min, mean = 0, sd = 1)) / sqrt(kappa) + 1
    return(q)
  } else if (type == 3) {
    lower <- 1/32
    upper <- 1/2
    tt <- u
    pin <- pbeta(lower, kappa, kappa) + u * (pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa))
    tt <- (qbeta(pin, kappa, kappa) - lower) / (upper - lower)
    return(tt)
  } else if (type == 4) {
    return(1 - qbeta(1 - u, 1/delta, 2)^(1/delta))
  } else if (type == 5) {
    return(1 - qbeta(1 - u^(2/kappa), 1/delta, 2)^(1/delta))
  } else if (type == 6) {
    dummy.func <- function(u, p, prob = NA, kappa = NA, delta = NA) {
      return(p.G(u = u, prob = prob, kappa = kappa, delta = delta, type = 6) - p)
    }
    find.root <- function(u, prob = NA, kappa = NA, delta = NA) {
      return(uniroot(dummy.func, interval = c(0, 1), p = u, prob = prob, kappa = kappa, delta = delta)$root)
    }
    return(sapply(u, FUN = find.root, prob = prob, kappa = kappa, delta = delta))
  }
}

#' Random generation from transformation G
#'
#' @param n number of samples
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#' @param type integer 1-6
#' @param unifsamp optional uniform samples
#' @param direct logical for type 6
#'
#' @return Random samples
#' @export
r.G <- function(n, prob = NA, kappa = NA, delta = NA, type = 1, unifsamp = NULL, direct = FALSE) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1, 2, 3, 5, 6) && identical(kappa, NA)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && identical(delta, NA)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && identical(prob, NA)) {
    stop("Argument `prob' missing.")
  }
  if (is.null(unifsamp)) {
    unifsamp <- runif(n)
  }
  if (type != 6 | (type == 6 & direct)) {
    return(q.G(unifsamp, prob = prob, kappa = kappa, delta = delta, type = type))
  } else if (type == 6 & !direct) {
    components <- sample(x = c(1, 2), size = n, replace = TRUE, prob = c(prob, 1 - prob))
    res <- numeric(n)
    res[components == 1] <- q.G(unifsamp[components == 1], prob = NA, kappa = kappa, delta = NA, type = 1)
    res[components == 2] <- q.G(unifsamp[components == 2], prob = NA, kappa = delta, delta = NA, type = 1)
    return(res)
  }
}

############################################
##    Continuous EGPD                     ##
############################################

#' CDF of the Extended GPD
#'
#' @param q quantiles
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return CDF values
#' @export
pegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(p.G(.pgpd_std(q, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type))
}

#' Density of the Extended GPD
#'
#' @param x values
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param log logical: return log-density?
#'
#' @return Density values
#' @export
degpd_density <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
  if (log == FALSE) {
    return(d.G(.pgpd_std(x, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type) *
             .dgpd_std(x, sigma = sigma, xi = xi))
  } else {
    return(d.G(.pgpd_std(x, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type, log = TRUE) +
             .dgpd_std(x, sigma = sigma, xi = xi, log = TRUE))
  }
}

#' Quantile function of the Extended GPD
#'
#' @param p probabilities
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return Quantile values
#' @export
qegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(.qgpd_std(q.G(p, prob = prob, kappa = kappa, delta = delta, type = type), sigma = sigma, xi = xi))
}

#' Random generation from the Extended GPD
#'
#' @param n number of samples
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param unifsamp optional uniform samples
#'
#' @return Random samples
#' @export
regpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  return(.qgpd_std(r.G(n, prob = prob, kappa = kappa, delta = delta, type = type, unifsamp = unifsamp), sigma = sigma, xi = xi))
}

############################################
##   Discrete Extended GPD               ##
############################################

#' Density of the Discrete Extended GPD
#'
#' @param x non-negative integer values
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return PMF values
#' @export
ddiscegpd <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(pegpd(x + 1, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type) -
         pegpd(x, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type))
}

#' CDF of the Discrete Extended GPD
#'
#' @param q non-negative integer quantiles
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return CDF values
#' @export
pdiscegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(pegpd(q + 1, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type))
}

#' Quantile function of the Discrete Extended GPD
#'
#' @param p probabilities
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return Quantile values (non-negative integers)
#' @export
qdiscegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  q <- ceiling(qegpd(p, prob = prob, kappa = kappa, delta = delta, type = type, sigma = sigma, xi = xi)) - 1
  q[q < 0] <- 0
  if (is.matrix(p)) {
    q <- matrix(q, ncol = ncol(p), nrow = nrow(p))
    colnames(q) <- colnames(p)
    rownames(q) <- rownames(p)
  }
  return(q)
}

#' Random generation from the Discrete Extended GPD
#'
#' @param n number of samples
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param unifsamp optional uniform samples
#'
#' @return Random non-negative integer samples
#' @export
rdiscegpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  return(floor(regpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type, unifsamp = unifsamp)))
}

############################################
##   Zero-Inflated Discrete Extended GPD  ##
############################################

#' Density of the Zero-Inflated Discrete Extended GPD
#'
#' @param x non-negative integer values
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return PMF values
#' @export
dzidiscegpd <- function(x, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  ly <- max(length(x), length(pi), length(sigma), length(xi))
  x <- rep(x, length = ly)
  pi <- rep(pi, length = ly)
  fy <- ifelse(x == 0,
    pi + (1 - pi) * ddiscegpd(0, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type),
    (1 - pi) * ddiscegpd(x, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type))
  fy
}

#' CDF of the Zero-Inflated Discrete Extended GPD
#'
#' @param q non-negative integer quantiles
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return CDF values
#' @export
pzidiscegpd <- function(q, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  ly <- max(length(q), length(pi), length(sigma), length(xi))
  q <- rep(q, length = ly)
  pi <- rep(pi, length = ly)
  cdf <- pdiscegpd(q, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type)
  cdf <- ifelse(q < 0, 0, pi + (1 - pi) * cdf)
  cdf
}

#' Quantile function of the Zero-Inflated Discrete Extended GPD
#'
#' @param p probabilities
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return Quantile values (non-negative integers)
#' @export
qzidiscegpd <- function(p, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  ly <- max(length(p), length(pi), length(sigma), length(xi))
  p <- rep(p, length = ly)
  pi <- rep(pi, length = ly)
  pnew <- ((p - pi) / (1 - pi)) - (1e-7)
  pnew <- ifelse(pnew > 0, pnew, 0)
  q <- qdiscegpd(pnew, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type)
  q
}

#' Random generation from the Zero-Inflated Discrete Extended GPD
#'
#' @param n number of samples
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param unifsamp optional uniform samples
#'
#' @return Random non-negative integer samples
#' @export
rzidiscegpd <- function(n, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  z <- rbinom(n, size = 1, prob = pi)
  y <- (1 - z) * rdiscegpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type, unifsamp = NULL)
  return(y)
}

############################################
##  Continuous Zero-Inflated EGPD         ##
############################################

#' Density of the Zero-Inflated Extended GPD
#'
#' Continuous zero-inflated EGPD density:
#' \eqn{h(x) = \pi I(x=0) + (1-\pi) f_{EGPD}(x)} for \eqn{x \ge 0}.
#'
#' @param x non-negative values
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param log logical: return log-density?
#'
#' @return Density values
#' @export
dziegpd <- function(x, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
  ly <- max(length(x), length(pi), length(sigma), length(xi))
  x <- rep(x, length = ly)
  pi <- rep(pi, length = ly)
  if (!log) {
    d <- ifelse(x == 0, pi,
                (1 - pi) * degpd_density(x, prob = prob, kappa = kappa,
                                         delta = delta, sigma = sigma,
                                         xi = xi, type = type))
    d[x < 0] <- 0
  } else {
    d <- ifelse(x == 0, log(pi),
                log(1 - pi) + degpd_density(x, prob = prob, kappa = kappa,
                                            delta = delta, sigma = sigma,
                                            xi = xi, type = type, log = TRUE))
    d[x < 0] <- -Inf
  }
  return(d)
}

#' CDF of the Zero-Inflated Extended GPD
#'
#' @param q quantiles
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return CDF values
#' @export
pziegpd <- function(q, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  ly <- max(length(q), length(pi), length(sigma), length(xi))
  q <- rep(q, length = ly)
  pi <- rep(pi, length = ly)
  cdf <- pegpd(q, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type)
  cdf <- ifelse(q < 0, 0, pi + (1 - pi) * cdf)
  cdf
}

#' Quantile function of the Zero-Inflated Extended GPD
#'
#' @param p probabilities
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return Quantile values
#' @export
qziegpd <- function(p, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  ly <- max(length(p), length(pi), length(sigma), length(xi))
  p <- rep(p, length = ly)
  pi <- rep(pi, length = ly)
  pnew <- ((p - pi) / (1 - pi)) - 1e-7
  pnew <- ifelse(pnew > 0, pnew, 0)
  qegpd(pnew, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type)
}

#' Random generation from the Zero-Inflated Extended GPD
#'
#' @param n number of samples
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return Random non-negative samples
#' @export
rziegpd <- function(n, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  z <- rbinom(n, size = 1, prob = pi)
  y <- (1 - z) * regpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type)
  return(y)
}
