#' Predictions from a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param newdata a data frame
#' @param type a character string: "link", "response", "lpmatrix", or "quantile"
#' @param prob a scalar or vector of probabilities for quantile estimation
#' @param se.fit logical: should standard errors be returned? Defaults to FALSE
#' @param marginal logical: should uncertainty integrate out smoothing parameter uncertainty? Defaults to TRUE
#' @param trace an integer controlling output verbosity
#' @param ... unused
#'
#' @return A data frame, list, or design matrix depending on \code{type}
#'
#' @export
predict.egpd <- function(object, newdata, type="link", prob=NULL, se.fit=FALSE, marginal=TRUE,
trace = 0, ...) {

## a few checks
family <- object$family

if (family == "egpd") {
  egpd_m <- object$likfns$m
  egpd_iG <- object$likfns$iG
}
if (family == "degpd") {
  degpd_m <- object$likfns$m
  degpd_iG <- object$likfns$iG
}
if (family == "zidegpd") {
  zidegpd_m <- object$likfns$m
  zidegpd_iG <- object$likfns$iG
}
if (family == "custom") {
  q_fn <- object$likfns$q
  unlink_fns <- object$likfns$unlink
}
if (family == "comppareto") {
  q_fn <- object$likfns$q
  p_fn <- object$likfns$p
}

if (!is.null(prob))
  type <- "quantile"
if (type == "quantile" & is.null(prob))
  stop("non-NULL `prob' required if `type = quantile'")

## end checks

## standard error setup
if (se.fit) {
  if (marginal) {
    V.type <- "Vc"
  } else {
    V.type <- "Vp"
  }
  conf.pars <- list(object$coefficients, object[[V.type]], object$idpars)
}

got.newdata <- !missing(newdata)

if (got.newdata) {
  pred.vars <- object$predictor.names
  missing.vars <- pred.vars[!(pred.vars %in% names(newdata))]
  if (length(missing.vars) > 0)
    stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `newdata'.", sep=""))
}

if (got.newdata) {
  ndat <- nrow(newdata)
} else {
  if (is.null(object$data))
    stop("Supply `egpd' with `removeData = FALSE' if not supplying `newdata'.")
  ndat <- nrow(object$data)
}

## X creation
X <- .X.egpd(object, newdata)
nX <- length(X)
nms <- names(object)[seq_len(nX)]

## offsets
if (got.newdata) {
  offsets <- attr(X, "offsets")
} else {
  offsets <- object$likdata$offsets
}

if (type == "lpmatrix") {
  attr(X, "offsets") <- offsets
  return(X)
} else {

out <- lapply(seq_len(nX), function(i) {
  eta <- X[[i]] %*% object[[i]]$coefficients
  if (length(offsets[[i]]) > 0) eta <- eta + offsets[[i]]
  eta
})
names(out) <- names(X)
out <- as.data.frame(lapply(out, function(x) x[,1]))

if (se.fit) {
  if (type != "quantile") {
    std.err <- lapply(seq_len(nX), function(i) sqrt(rowSums(X[[i]] * (X[[i]] %*% object[[i]][[V.type]]))))
    std.err <- as.data.frame(std.err)
    names(std.err) <- nms
  }
}

if (type %in% c("response", "quantile")) {

  if (family == "custom") {
    for (i in seq_along(nms)) {
      if (se.fit) {
        if (!is.null(attr(unlink_fns[[i]], "deriv")))
          std.err[, i] <- attr(unlink_fns[[i]], "deriv")(out[, i]) * std.err[, i]
      }
      if (!is.null(unlink_fns[[i]]))
        out[, i] <- unlink_fns[[i]](out[, i])
    }
  } else {
    unlink <- which(substr(nms, 1, 3) == "log")
    for (i in unlink) {
      if (substr(nms[i], 1, 5) == "logit") {
        ## logit link: pi = 1/(1+exp(-eta)), d(pi)/d(eta) = pi*(1-pi)
        out[, i] <- 1 / (1 + exp(-out[, i]))
        if (se.fit & type == "response")
          std.err[, i] <- out[, i] * (1 - out[, i]) * std.err[, i]
      } else {
        ## log link: theta = exp(eta), d(theta)/d(eta) = theta
        out[, i] <- exp(out[, i])
        if (se.fit & type == "response")
          std.err[, i] <- out[, i] * std.err[, i]
      }
    }
  }

  nms <- sub("^cloglog", "", nms)
  nms <- sub("^probit", "", nms)
  nms <- sub("^logit", "", nms)
  nms <- sub("^log", "", nms)
  names(out) <- nms

  if (se.fit & type == "response")
    names(std.err) <- nms

  ## Match the boundary guard applied inside the likelihood (src/gpig.cpp,
  ## eta_to_abc): the GPIG/ZIGPIG tail exponent a and down-weighting c are
  ## clamped to the same safe interior. In the likelihood a flat region beyond
  ## the clamp lets the optimiser leave eta unbounded, so the plain inverse link
  ## above can reconstruct an a/c marginally outside the bound; clamping here
  ## makes predict() (and the quantiles below) agree with the fitted model.
  if (family %in% c("gpig", "gpignat", "zigpig", "zigpignat")) {
    if ("a" %in% names(out)) out[["a"]] <- pmin(pmax(out[["a"]], 1e-3), 1 - 1e-6)
    if ("c" %in% names(out)) out[["c"]] <- pmin(pmax(out[["c"]], 1e-6), 1 - 1e-6)
  }
  ## The same guard for the count families of src/countfams.cpp. It matters most for the
  ## generalized Waring: rho -> Inf is its negative-binomial limit, and since rho enters the
  ## likelihood only through the clamp, the optimiser can push eta3 arbitrarily far with no
  ## change in fit. Unclamped, predict() then reports rho = Inf (or 1e220) while the fitted
  ## model actually used rho = 1e6 -- a boundary fit that looks like a converged one.
  if (family %in% c("gpois", "gwaring", "plnorm")) {
    bd <- cf_bounds_cpp()[[family]]        # single source of truth: the C++ clamps
    for (pn in c("mu", "lambda", "k", "rho", "sigma")) {
      lo <- paste0(pn, "_lo"); hi <- paste0(pn, "_hi")
      if (pn %in% names(out) && lo %in% names(bd))
        out[[pn]] <- pmin(pmax(out[[pn]], bd[[lo]]), bd[[hi]])
    }
  }

  ## issue #2: report the *bounded* tail index for DEGPD model 1 fitted with a
  ## finite xi.max, matching the link used in the likelihood. The unlink above
  ## produced exp(eta); map it through xi = xi.max*(1 - exp(-exp(eta)/xi.max)).
  ## Applied before the quantile block so quantiles use the same bounded xi.
  if (isTRUE(object$likdata$bounded.xi) && "shape" %in% names(out)) {
    a <- object$likdata$xi.max
    out[["shape"]] <- a * (1 - exp(-out[["shape"]] / a))
  }

  if (type == "quantile") {

    pars <- out
    nprob <- length(prob)
    out <- matrix(NA, ndat, nprob)

    for (j in seq_len(nprob)) {
      pj <- prob[j]

      if (family == "egpd") {
        if (egpd_m %in% c(1, 3, 5, 6)) {
          pj <- egpd_iG(pj, pars[, 3])
        } else {
          if (egpd_m == 2) {
            # Reparameterization: col 4 is exp(ldkappa) = kappa2 - kappa1; need kappa2 = kappa1 + dkappa
            kappa2 <- pars[, 3] + pars[, 4]
            pj <- egpd_iG(pj, pars[, 3], kappa2, pars[, 5])
          } else {
            pj <- egpd_iG(pj, pars[, 3], pars[, 4])
          }
        }
      }

      if (family == "degpd") {
        if (degpd_m %in% c(1, 3, 5, 6)) {
          pj <- degpd_iG(pj, pars[, 3])
        } else {
          if (degpd_m == 2) {
            # Reparameterization: col 4 is exp(ldkappa) = kappa2 - kappa1; need kappa2 = kappa1 + dkappa
            kappa2 <- pars[, 3] + pars[, 4]
            pj <- degpd_iG(pj, pars[, 3], kappa2, pars[, 5])
          } else {
            pj <- degpd_iG(pj, pars[, 3], pars[, 4])
          }
        }
      }

      if (family == "zidegpd") {
        if (zidegpd_m %in% c(1, 3, 5, 6)) {
          pj <- ((pj - pars[, 4]) / (1 - pars[, 4])) - (1e-7)
          pj <- ifelse(pj > 0, pj, 0)
          pj <- zidegpd_iG(pj, pars[, 3])
        } else {
          if (zidegpd_m == 2) {
            pj <- ((pj - pars[, 6]) / (1 - pars[, 6])) - (1e-7)
            pj <- ifelse(pj > 0, pj, 0)
            # Reparameterization: col 4 is exp(ldkappa) = kappa2 - kappa1; need kappa2 = kappa1 + dkappa
            kappa2 <- pars[, 3] + pars[, 4]
            pj <- zidegpd_iG(pj, pars[, 3], kappa2, pars[, 5])
          } else {
            pj <- ((pj - pars[, 5]) / (1 - pars[, 5])) - (1e-7)
            pj <- ifelse(pj > 0, pj, 0)
            if (zidegpd_m == 4) {
              pj <- zidegpd_iG(pj, pars[, 4], pars[, 3])
            } else {
              pj <- zidegpd_iG(pj, pars[, 3], pars[, 4])
            }
          }
        }
      }

      if (family == "custom") {
        if (length(nms) == 1) {
          out[, j] <- q_fn(pj, pars[,1])
        } else if (length(nms) == 2) {
          out[, j] <- q_fn(pj, pars[,1], pars[,2])
        } else if (length(nms) == 3) {
          out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3])
        } else {
          out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3], pars[,4])
        }
      } else if (family == "comppareto") {
        if (length(nms) == 3) {
          out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3])
        } else {
          out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3], pars[,4])
        }
      } else {
        if (family %in% c("gpd", "egpd")) {
          out[, j] <- .qgpd(pj, 0, pars[,1], pars[,2])
        } else if (family == "degpd") {
          out[, j] <- ceiling(.qgpd(pj, 0, pars[,1], pars[,2])) - 1
          out[, j][out[, j] < 0] <- 0
        } else if (family == "zidegpd") {
          out[, j] <- ceiling(.qgpd(pj, 0, pars[,1], pars[,2])) - 1
          out[, j][out[, j] < 0] <- 0
        } else {
          stop("invalid family for quantile prediction")
        }
      }
    }

    if (se.fit) {
      if (family %in% c("egpd", "degpd", "zidegpd"))
        stop("Standard errors not yet available for this family.")
    }

    out <- as.data.frame(out)
    names(out) <- paste("q", round(prob, 3), sep=":")

  } ## end quantile

  ## Convenience alias: the GPD tail index is returned in the column named
  ## "shape" (from the logshape/shape link). Many users reach for `xi`; expose
  ## it as an alias so `pars[["xi"]]` does not silently return NULL (see #5).
  if (type == "response" && is.data.frame(out) &&
      "shape" %in% names(out) && !("xi" %in% names(out))) {
    out$xi <- out[["shape"]]
  }

  if (se.fit) {
    out <- list(fitted = out, se.fit = std.err)
  }

} ## end response/quantile

if (type == "link" & se.fit) {
  out <- list(fitted = out, se.fit = std.err)
}

return(out)

} ## end not lpmatrix

}


#' Randomized quantile residuals for a fitted \code{egpd} model
#'
#' Computes randomized quantile residuals (Dunn & Smyth, 1996) for
#' discrete or continuous EGPD models. For discrete models, a uniform
#' random variate is drawn between the lower and upper CDF bounds at
#' each observation and transformed to the normal scale via
#' \code{qnorm}. For continuous models, the CDF value is transformed
#' directly.
#'
#' @param object a fitted \code{egpd} object
#' @param seed optional random seed for reproducibility
#' @param ... unused
#'
#' @return A numeric vector of randomized quantile residuals on the
#'   standard normal scale
#'
#' @export
rqresid <- function(object, ...) UseMethod("rqresid")

#' @rdname rqresid
#' @export
rqresid.egpd <- function(object, seed = NULL, ...) {

  if (!is.null(seed)) set.seed(seed)
  if (is.null(object$data))
    stop("rqresid requires the original data; refit with removeData = FALSE", call. = FALSE)

  family <- object$family
  if (family == "comppareto") {
    p_fn <- object$likfns$p
    m <- NULL
    dtype <- NULL
  } else {
    m <- object$likfns$m
    type_map <- c(1L, 6L, 4L, 5L, 2L, 3L)
    dtype <- type_map[m]
  }

  ## response vector
  y <- object$data[[object$response.name]]
  n <- length(y)

  ## predicted parameters on response scale
  pars <- predict(object, type = "response")

  ## Count families (PIG/GPIG/Bell and their zero-inflated variants) use their
  ## own mean-convention CDFs rather than the EGPD sigma/xi/kappa structure.
  ## Randomized quantile residual: u ~ Unif(F(y-1), F(y)), r = qnorm(u).
  count_fams <- c("pig", "zipig", "gpig", "zigpig", "bell", "zibell")
  if (family %in% count_fams) {
    cdf <- switch(family,
      pig    = function(q) ppig(q, mu = pars$mu, sigma = pars$sigma),
      zipig  = function(q) pzipig(q, mu = pars$mu, sigma = pars$sigma, pi = pars$pi),
      gpig   = function(q) pGPIG(q, mu = pars$mu, sigma = pars$a, nu = pars$c),
      zigpig = function(q) pZIGPIG(q, mu = pars$mu, sigma = pars$a, nu = pars$c,
                                   tau = pars$pi),
      bell   = function(q) pBELL(q, mu = pars$mu),
      zibell = function(q) pZIBELL(q, mu = pars$mu, sigma = pars$pi))
    p_upper <- cdf(y)
    p_lower <- cdf(y - 1)
    p_lower[y == 0] <- 0
    r <- qnorm(runif(n, p_lower, p_upper))
    r[is.infinite(r)] <- NA
    return(r)
  }

  if (family != "comppareto") {
    sigma <- pars[[1]]
    xi    <- pars[[2]]

    ## build CDF arguments depending on model type
    cdf_args <- list(sigma = sigma, xi = xi, type = dtype)
    if (m == 1) {
      cdf_args$kappa <- pars[[3]]
    } else if (m == 2) {
      cdf_args$kappa <- pars[[3]]                  # kappa1
      cdf_args$delta <- pars[[3]] + pars[[4]]      # kappa2 = kappa1 + dkappa
      cdf_args$prob  <- pars[[5]]
    } else if (m == 3) {
      cdf_args$delta <- pars[[3]]
    } else if (m == 4) {
      cdf_args$delta <- pars[[4]]
      cdf_args$kappa <- pars[[3]]
    } else if (m %in% c(5, 6)) {
      cdf_args$kappa <- pars[[3]]
    }
  }

  if (family == "egpd") {
    ## continuous: residual = qnorm(F(y))
    cdf_args$q <- y
    Fu <- do.call(pegpd, cdf_args)
    r  <- qnorm(Fu)
  } else if (family == "degpd") {
    ## discrete: randomize between F(y-1) and F(y)
    cdf_args$q <- y
    p_upper <- do.call(pdiscegpd, cdf_args)
    cdf_args$q <- y - 1
    p_lower <- do.call(pdiscegpd, cdf_args)
    p_lower[y == 0] <- 0
    u <- runif(n, p_lower, p_upper)
    r <- qnorm(u)
  } else if (family == "zidegpd") {
    ## zero-inflated discrete
    if (m %in% c(1, 3, 5, 6)) {
      pi_val <- pars[[4]]
    } else if (m == 2) {
      pi_val <- pars[[6]]
    } else {
      pi_val <- pars[[5]]
    }
    cdf_args$pi <- pi_val
    cdf_args$q  <- y
    p_upper <- do.call(pzidiscegpd, cdf_args)
    cdf_args$q <- y - 1
    p_lower <- do.call(pzidiscegpd, cdf_args)
    p_lower[y == 0] <- 0
    u <- runif(n, p_lower, p_upper)
    r <- qnorm(u)
  } else if (family == "comppareto") {
    if (length(pars) == 3L) {
      Fu <- p_fn(y, pars[[1]], pars[[2]], pars[[3]])
    } else {
      Fu <- p_fn(y, pars[[1]], pars[[2]], pars[[3]], pars[[4]])
    }
    r <- qnorm(Fu)
  } else {
    stop("rqresid not implemented for family '", family, "'")
  }

  r[is.infinite(r)] <- NA
  return(r)
}
