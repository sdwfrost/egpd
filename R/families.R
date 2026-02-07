.setup.family.egpd <- function(family, egpd, degpd, zidegpd, formula, likfns) {
  if (family == "egpd") {
    if (is.null(egpd$m))
      egpd$m <- 1
    if (egpd$m == 1) {
      lik.fns <- .egpd1fns
      npar <- 3
      nms <- c("lpsi", "xi", "lkappa")
      nms2 <- c('logscale', 'shape', 'logkappa')
      attr(family, "type") <- 1
    } else {
      if (egpd$m == 2) {
        lik.fns <- .egpd2fns
        npar <- 5
        nms <- c("lpsi", "xi", "lkappa1", "ldkappa", "logitp")
        nms2 <- c('logscale', 'shape', 'logkappa1', 'logdkappa', 'logitp')
        attr(family, "type") <- 6
      } else {
        if (egpd$m == 3) {
          lik.fns <- .egpd3fns
          npar <- 3
          nms <- c("lpsi", "xi", "ldelta")
          nms2 <- c('logscale', 'shape', 'logdelta')
          attr(family, "type") <- 4
        } else {
          lik.fns <- .egpd4fns
          npar <- 4
          nms <- c("lpsi", "xi", "ldelta", "lkappa")
          nms2 <- c('logscale', 'shape', 'logdelta', 'logkappa')
          attr(family, "type") <- 5
        }
      }
    }
  } else if (family == "degpd") {
    if (is.null(degpd$m))
      degpd$m <- 1
    if (degpd$m == 1) {
      lik.fns <- .degpd1fns
      npar <- 3
      nms <- c("lsigma", "lxi", "lkappa")
      nms2 <- c('logscale', 'logshape', 'logkappa')
      attr(family, "type") <- 1
    } else {
      if (degpd$m == 2) {
        lik.fns <- .degpd2fns
        npar <- 5
        nms <- c("lsigma", "lxi", "lkappa1", "ldkappa", "logitp")
        nms2 <- c('logscale', 'logshape', 'logkappa1', 'logdkappa', 'logitp')
        attr(family, "type") <- 6
      } else {
        if (degpd$m == 3) {
          lik.fns <- .degpd3fns
          npar <- 3
          nms <- c("lsigma", "lxi", "ldelta")
          nms2 <- c('logscale', 'logshape', 'logdelta')
          attr(family, "type") <- 4
        } else {
          lik.fns <- .degpd4fns
          npar <- 4
          nms <- c("lsigma", "lxi", "ldelta", "lkappa")
          nms2 <- c('logscale', 'logshape', 'logdelta', 'logkappa')
          attr(family, "type") <- 5
        }
      }
    }
  } else if (family == "zidegpd") {
    if (is.null(zidegpd$m))
      zidegpd$m <- 1
    if (zidegpd$m == 1) {
      lik.fns <- .zidegpd1fns
      npar <- 4
      nms <- c("lsigma", "lxi", "lkappa", "logitpi")
      nms2 <- c('logscale', 'logshape', 'logkappa', 'logitpi')
      attr(family, "type") <- 1
    } else {
      if (zidegpd$m == 2) {
        lik.fns <- .zidegpd2fns
        npar <- 6
        nms <- c("lsigma", "lxi", "lkappa1", "ldkappa", "logitp", "logitpi")
        nms2 <- c('logscale', 'logshape', 'logkappa1', 'logdkappa', 'logitp', 'logitpi')
        attr(family, "type") <- 6
      } else {
        if (zidegpd$m == 3) {
          lik.fns <- .zidegpd3fns
          npar <- 4
          nms <- c("lsigma", "lxi", "ldelta", "logitpi")
          nms2 <- c('logscale', 'logshape', 'logdelta', 'logitpi')
          attr(family, "type") <- 4
        } else {
          lik.fns <- .zidegpd4fns
          npar <- 5
          nms <- c("lsigma", "lxi", "lkappa", "ldelta", "logitpi")
          nms2 <- c('logscale', 'logshape', 'logkappa', 'logdelta', 'logitpi')
          attr(family, "type") <- 5
        }
      }
    }
  } else {
    if (length(likfns)) {
      lik.fns <- likfns
      family <- "custom"
      npar <- length(formula)
      if (is.null(names(formula))) {
        nms <- paste("par", seq_along(formula), sep = "_")
      } else {
        nms <- names(formula)
      }
      nms2 <- nms
    } else {
      stop(paste("Family '", family, "' not supported. Use 'egpd', 'degpd', or 'zidegpd'.", sep=""))
    }
  }
  out <- list(npar=npar, npar2=npar, lik.fns=lik.fns, nms=nms, family=family, nms2 = nms2)
}
