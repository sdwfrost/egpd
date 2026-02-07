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

  nms <- gsub("cloglog", "", nms)
  nms <- gsub("probit", "", nms)
  nms <- gsub("logit", "", nms)
  nms <- gsub("log", "", nms)
  names(out) <- nms

  if (se.fit & type == "response")
    names(std.err) <- nms

  if (type == "quantile") {

    pars <- out
    nprob <- length(prob)
    out <- matrix(NA, ndat, nprob)

    for (j in seq_len(nprob)) {
      pj <- prob[j]

      if (family == "egpd") {
        if (egpd_m %in% c(1, 3)) {
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
        if (degpd_m %in% c(1, 3)) {
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
        if (zidegpd_m %in% c(1, 3)) {
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
            pj <- zidegpd_iG(pj, pars[, 3], pars[, 4])
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

  family <- object$family
  m <- object$likfns$m
  type_map <- c(1L, 6L, 4L, 5L)
  dtype <- type_map[m]

  ## response vector
  y <- object$data[[object$response.name]]
  n <- length(y)

  ## predicted parameters on response scale
  pars <- predict(object, type = "response")

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
    cdf_args$delta <- pars[[3]]
    cdf_args$kappa <- pars[[4]]
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
    if (m %in% c(1, 3)) {
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
  } else {
    stop("rqresid not implemented for family '", family, "'")
  }

  r[is.infinite(r)] <- NA
  return(r)
}
