## bamlss family constructors for EGPD models
##
## These allow fitting EGPD and zero-inflated EGPD models using the
## bamlss package's Bayesian distributional regression framework.
## Distribution functions are reused from the egpd package; no
## external distribution packages are needed.

## internal helper (mirrors bamlss:::parse.links)
.parse_links <- function(links, default.links, ...) {
  dots <- list(...)
  nl <- names(default.links)
  if (length(dots))
    links <- as.character(dots)
  if (is.null(names(links)))
    names(links) <- rep(nl, length.out = length(links))
  links <- as.list(links)
  for (j in nl) {
    if (is.null(links[[j]]))
      links[[j]] <- default.links[j]
  }
  links <- links[nl]
  links <- as.character(links)
  names(links) <- nl
  links
}

## m -> internal type mapping (same as families.R)
.m_to_type <- function(m) c(1L, 6L, 4L, 5L, 2L, 3L)[m]


#' bamlss family for continuous EGPD
#'
#' Creates a \code{family.bamlss} object for fitting continuous
#' Extended Generalized Pareto Distribution models with \code{bamlss()}.
#'
#' @param m integer 1--4 selecting the G transformation:
#'   \describe{
#'     \item{1}{Power: \eqn{G(u) = u^\kappa}}
#'     \item{2}{Mixture: \eqn{G(u) = p u^{\kappa} + (1-p) u^{\delta}}}
#'     \item{3}{Incomplete beta}
#'     \item{4}{Power-beta}
#'   }
#' @param ... arguments passed to link specification
#'
#' @return An object of class \code{family.bamlss}
#' @export
egpd_bamlss <- function(m = 1, ...) {
  if (!m %in% 1:6) stop("m must be 1, 2, 3, 4, 5, or 6.")
  dtype <- .m_to_type(m)

  if (m == 1) {
    params <- c("sigma", "xi", "kappa")
    default_links <- c(sigma = "log", xi = "log", kappa = "log")
  } else if (m == 2) {
    params <- c("sigma", "xi", "kappa", "delta", "prob")
    default_links <- c(sigma = "log", xi = "log", kappa = "log",
                       delta = "log", prob = "logit")
  } else if (m == 3) {
    params <- c("sigma", "xi", "delta")
    default_links <- c(sigma = "log", xi = "log", delta = "log")
  } else if (m == 4) {
    params <- c("sigma", "xi", "delta", "kappa")
    default_links <- c(sigma = "log", xi = "log", delta = "log", kappa = "log")
  } else if (m %in% c(5, 6)) {
    params <- c("sigma", "xi", "kappa")
    default_links <- c(sigma = "log", xi = "log", kappa = "log")
  }

  rval <- list(
    "family" = paste0("egpd", m),
    "names"  = params,
    "links"  = .parse_links(default_links, default_links, ...),

    "d" = function(y, par, log = FALSE) {
      degpd_density(y,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype, log = log)
    },

    "p" = function(y, par, ...) {
      pegpd(y,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "q" = function(p, par, ...) {
      qegpd(p,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "r" = function(n, par) {
      regpd(n,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    }
  )

  ## initialization
  if (m %in% c(1, 5, 6)) {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y))
    )
  } else if (m == 2) {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      delta = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      prob  = function(y, ...) rep(0.5, length(y))
    )
  } else if (m == 3) {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y))
    )
  } else {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y))
    )
  }

  class(rval) <- "family.bamlss"
  rval
}


#' bamlss family for zero-inflated continuous EGPD
#'
#' Creates a \code{family.bamlss} object for fitting continuous
#' zero-inflated Extended Generalized Pareto Distribution models
#' with \code{bamlss()}.
#'
#' @param m integer 1--4 selecting the G transformation (see
#'   \code{\link{egpd_bamlss}} for details)
#' @param ... arguments passed to link specification
#'
#' @return An object of class \code{family.bamlss}
#' @export
ziegpd_bamlss <- function(m = 1, ...) {
  if (!m %in% 1:6) stop("m must be 1, 2, 3, 4, 5, or 6.")
  dtype <- .m_to_type(m)

  if (m == 1) {
    params <- c("sigma", "xi", "kappa", "pi")
    default_links <- c(sigma = "log", xi = "log", kappa = "log", pi = "logit")
  } else if (m == 2) {
    params <- c("sigma", "xi", "kappa", "delta", "prob", "pi")
    default_links <- c(sigma = "log", xi = "log", kappa = "log",
                       delta = "log", prob = "logit", pi = "logit")
  } else if (m == 3) {
    params <- c("sigma", "xi", "delta", "pi")
    default_links <- c(sigma = "log", xi = "log", delta = "log", pi = "logit")
  } else if (m == 4) {
    params <- c("sigma", "xi", "delta", "kappa", "pi")
    default_links <- c(sigma = "log", xi = "log", delta = "log",
                       kappa = "log", pi = "logit")
  } else if (m %in% c(5, 6)) {
    params <- c("sigma", "xi", "kappa", "pi")
    default_links <- c(sigma = "log", xi = "log", kappa = "log", pi = "logit")
  }

  rval <- list(
    "family" = paste0("ziegpd", m),
    "names"  = params,
    "links"  = .parse_links(default_links, default_links, ...),

    "d" = function(y, par, log = FALSE) {
      dziegpd(y,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype, log = log)
    },

    "p" = function(y, par, ...) {
      pziegpd(y,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "q" = function(p, par, ...) {
      qziegpd(p,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "r" = function(n, par) {
      rziegpd(n,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    }
  )

  ## initialization
  npar <- length(params)
  pi_init <- function(y, ...) rep(max(mean(y == 0), 0.01), length(y))
  base <- list(
    sigma = function(y, ...) rep(max(sd(y[y > 0]), 1e-6), length(y)),
    xi    = function(y, ...) rep(0.3, length(y))
  )

  if (m %in% c(1, 5, 6)) {
    rval$initialize <- c(base, list(
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      pi = pi_init))
  } else if (m == 2) {
    rval$initialize <- c(base, list(
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      delta = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      prob  = function(y, ...) rep(0.5, length(y)),
      pi = pi_init))
  } else if (m == 3) {
    rval$initialize <- c(base, list(
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y)),
      pi = pi_init))
  } else {
    rval$initialize <- c(base, list(
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      pi = pi_init))
  }

  class(rval) <- "family.bamlss"
  rval
}


#' bamlss family for discrete EGPD (DEGPD)
#'
#' Creates a \code{family.bamlss} object for fitting discrete
#' Extended Generalized Pareto Distribution models with \code{bamlss()}.
#'
#' @param m integer 1--4 selecting the G transformation (see
#'   \code{\link{egpd_bamlss}} for details)
#' @param ... arguments passed to link specification
#'
#' @return An object of class \code{family.bamlss}
#' @export
degpd_bamlss <- function(m = 1, ...) {
  if (!m %in% 1:6) stop("m must be 1, 2, 3, 4, 5, or 6.")
  dtype <- .m_to_type(m)

  if (m == 1) {
    params <- c("sigma", "xi", "kappa")
    default_links <- c(sigma = "log", xi = "log", kappa = "log")
  } else if (m == 2) {
    params <- c("sigma", "xi", "kappa", "delta", "prob")
    default_links <- c(sigma = "log", xi = "log", kappa = "log",
                       delta = "log", prob = "logit")
  } else if (m == 3) {
    params <- c("sigma", "xi", "delta")
    default_links <- c(sigma = "log", xi = "log", delta = "log")
  } else if (m == 4) {
    params <- c("sigma", "xi", "delta", "kappa")
    default_links <- c(sigma = "log", xi = "log", delta = "log", kappa = "log")
  } else if (m %in% c(5, 6)) {
    params <- c("sigma", "xi", "kappa")
    default_links <- c(sigma = "log", xi = "log", kappa = "log")
  }

  rval <- list(
    "family" = paste0("degpd", m),
    "names"  = params,
    "links"  = .parse_links(default_links, default_links, ...),

    "d" = function(y, par, log = FALSE) {
      d <- ddiscegpd(y,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
      if (log) d <- ifelse(d > 0, log(d), -Inf)
      d
    },

    "p" = function(y, par, ...) {
      pdiscegpd(y,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "q" = function(p, par, ...) {
      qdiscegpd(p,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "r" = function(n, par) {
      rdiscegpd(n,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    }
  )

  ## initialization
  if (m %in% c(1, 5, 6)) {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y))
    )
  } else if (m == 2) {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      delta = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      prob  = function(y, ...) rep(0.5, length(y))
    )
  } else if (m == 3) {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y))
    )
  } else {
    rval$initialize <- list(
      sigma = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      xi    = function(y, ...) rep(0.3, length(y)),
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y))
    )
  }

  class(rval) <- "family.bamlss"
  rval
}


#' bamlss family for zero-inflated discrete EGPD (ZIDEGPD)
#'
#' Creates a \code{family.bamlss} object for fitting zero-inflated
#' discrete Extended Generalized Pareto Distribution models with
#' \code{bamlss()}.
#'
#' @param m integer 1--4 selecting the G transformation (see
#'   \code{\link{egpd_bamlss}} for details)
#' @param ... arguments passed to link specification
#'
#' @return An object of class \code{family.bamlss}
#' @export
zidegpd_bamlss <- function(m = 1, ...) {
  if (!m %in% 1:6) stop("m must be 1, 2, 3, 4, 5, or 6.")
  dtype <- .m_to_type(m)

  if (m == 1) {
    params <- c("sigma", "xi", "kappa", "pi")
    default_links <- c(sigma = "log", xi = "log", kappa = "log", pi = "logit")
  } else if (m == 2) {
    params <- c("sigma", "xi", "kappa", "delta", "prob", "pi")
    default_links <- c(sigma = "log", xi = "log", kappa = "log",
                       delta = "log", prob = "logit", pi = "logit")
  } else if (m == 3) {
    params <- c("sigma", "xi", "delta", "pi")
    default_links <- c(sigma = "log", xi = "log", delta = "log", pi = "logit")
  } else if (m == 4) {
    params <- c("sigma", "xi", "delta", "kappa", "pi")
    default_links <- c(sigma = "log", xi = "log", delta = "log",
                       kappa = "log", pi = "logit")
  } else if (m %in% c(5, 6)) {
    params <- c("sigma", "xi", "kappa", "pi")
    default_links <- c(sigma = "log", xi = "log", kappa = "log", pi = "logit")
  }

  rval <- list(
    "family" = paste0("zidegpd", m),
    "names"  = params,
    "links"  = .parse_links(default_links, default_links, ...),

    "d" = function(y, par, log = FALSE) {
      d <- dzidiscegpd(y,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
      if (log) d <- ifelse(d > 0, log(d), -Inf)
      d
    },

    "p" = function(y, par, ...) {
      pzidiscegpd(y,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "q" = function(p, par, ...) {
      qzidiscegpd(p,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    },

    "r" = function(n, par) {
      rzidiscegpd(n,
        pi    = par$pi,
        sigma = par$sigma, xi = par$xi,
        kappa = if (!is.null(par$kappa)) par$kappa else NA,
        delta = if (!is.null(par$delta)) par$delta else NA,
        prob  = if (!is.null(par$prob))  par$prob  else NA,
        type = dtype)
    }
  )

  ## initialization
  pi_init <- function(y, ...) rep(max(mean(y == 0), 0.01), length(y))
  base <- list(
    sigma = function(y, ...) rep(max(sd(y[y > 0]), 1e-6), length(y)),
    xi    = function(y, ...) rep(0.3, length(y))
  )

  if (m %in% c(1, 5, 6)) {
    rval$initialize <- c(base, list(
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      pi = pi_init))
  } else if (m == 2) {
    rval$initialize <- c(base, list(
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      delta = function(y, ...) rep(max(sd(y), 1e-6), length(y)),
      prob  = function(y, ...) rep(0.5, length(y)),
      pi = pi_init))
  } else if (m == 3) {
    rval$initialize <- c(base, list(
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y)),
      pi = pi_init))
  } else {
    rval$initialize <- c(base, list(
      delta = function(y, ...) rep(max(0.5 * mean(y), 1e-6), length(y)),
      kappa = function(y, ...) rep(max(mean(y), 1e-6), length(y)),
      pi = pi_init))
  }

  class(rval) <- "family.bamlss"
  rval
}
