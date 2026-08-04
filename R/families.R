## Is an identity link on xi requested? Available for every DEGPD carrier: model 1
## through the symbolically-derived src/degpd1_identity.cpp, models 2-6 by handing
## xi straight to the existing C++ and converting the returned log-link derivatives
## in .identity_xi_chain() (see the note in R/degpd_lik.R).
##
## This validates rather than ignores. Before, an unrecognised link was silently
## dropped, so degpd.args = list(m = 3, link = "identity") quietly fitted the LOG
## link and returned a strictly positive xi -- the caller had asked for the ability
## to represent a bounded tail and got no indication they had not received it.
.degpd_identity_link <- function(degpd) {
  if (is.null(degpd$link)) return(FALSE)
  if (!is.character(degpd$link) || length(degpd$link) != 1L ||
      !degpd$link %in% c("log", "identity"))
    stop("degpd.args$link must be \"log\" or \"identity\".", call. = FALSE)
  identical(degpd$link, "identity")
}

.setup.family.egpd <- function(family, egpd, degpd, zidegpd, comppareto, formula, likfns) {
  if (family == "egpd") {
    if (is.null(egpd$m))
      egpd$m <- 1
    if (egpd$m == 1) {
      lik.fns <- .egpd1fns
      npar <- 3
      nms <- c("lpsi", "xi", "lkappa")
      nms2 <- c('logscale', 'shape', 'logkappa')
      attr(family, "type") <- 1
    } else if (egpd$m == 2) {
      lik.fns <- .egpd2fns
      npar <- 5
      nms <- c("lpsi", "xi", "lkappa1", "ldkappa", "logitp")
      nms2 <- c('logscale', 'shape', 'logkappa1', 'logdkappa', 'logitp')
      attr(family, "type") <- 6
    } else if (egpd$m == 3) {
      lik.fns <- .egpd3fns
      npar <- 3
      nms <- c("lpsi", "xi", "ldelta")
      nms2 <- c('logscale', 'shape', 'logdelta')
      attr(family, "type") <- 4
    } else if (egpd$m == 4) {
      lik.fns <- .egpd4fns
      npar <- 4
      nms <- c("lpsi", "xi", "ldelta", "lkappa")
      nms2 <- c('logscale', 'shape', 'logdelta', 'logkappa')
      attr(family, "type") <- 5
    } else if (egpd$m == 5) {
      lik.fns <- .egpd5fns
      npar <- 3
      nms <- c("lpsi", "xi", "lkappa")
      nms2 <- c('logscale', 'shape', 'logkappa')
      attr(family, "type") <- 2
    } else if (egpd$m == 6) {
      lik.fns <- .egpd6fns
      npar <- 3
      nms <- c("lpsi", "xi", "lkappa")
      nms2 <- c('logscale', 'shape', 'logkappa')
      attr(family, "type") <- 3
    } else {
      stop("egpd$m must be 1, 2, 3, 4, 5, or 6.")
    }
  } else if (family == "degpd") {
    if (is.null(degpd$m))
      degpd$m <- 1
    if (degpd$m == 1) {
      if (.degpd_identity_link(degpd)) {
        lik.fns <- .degpd1idfns                       # identity link on xi (xi may be < 0)
        npar <- 3
        nms <- c("lsigma", "xi", "lkappa")
        nms2 <- c('logscale', 'shape', 'logkappa')
      } else {
        lik.fns <- .degpd1fns
        npar <- 3
        nms <- c("lsigma", "lxi", "lkappa")
        nms2 <- c('logscale', 'logshape', 'logkappa')
      }
      attr(family, "type") <- 1
    } else if (degpd$m == 2) {
      idl <- .degpd_identity_link(degpd)
      lik.fns <- if (idl) .degpd2idfns else .degpd2fns
      npar <- 5
      nms <- c("lsigma", if (idl) "xi" else "lxi", "lkappa1", "ldkappa", "logitp")
      nms2 <- c('logscale', if (idl) 'shape' else 'logshape', 'logkappa1', 'logdkappa', 'logitp')
      attr(family, "type") <- 6
    } else if (degpd$m == 3) {
      idl <- .degpd_identity_link(degpd)
      lik.fns <- if (idl) .degpd3idfns else .degpd3fns
      npar <- 3
      nms <- c("lsigma", if (idl) "xi" else "lxi", "ldelta")
      nms2 <- c('logscale', if (idl) 'shape' else 'logshape', 'logdelta')
      attr(family, "type") <- 4
    } else if (degpd$m == 4) {
      idl <- .degpd_identity_link(degpd)
      lik.fns <- if (idl) .degpd4idfns else .degpd4fns
      npar <- 4
      nms <- c("lsigma", if (idl) "xi" else "lxi", "ldelta", "lkappa")
      nms2 <- c('logscale', if (idl) 'shape' else 'logshape', 'logdelta', 'logkappa')
      attr(family, "type") <- 5
    } else if (degpd$m == 5) {
      idl <- .degpd_identity_link(degpd)
      lik.fns <- if (idl) .degpd5idfns else .degpd5fns
      npar <- 3
      nms <- c("lsigma", if (idl) "xi" else "lxi", "lkappa")
      nms2 <- c('logscale', if (idl) 'shape' else 'logshape', 'logkappa')
      attr(family, "type") <- 2
    } else if (degpd$m == 6) {
      idl <- .degpd_identity_link(degpd)
      lik.fns <- if (idl) .degpd6idfns else .degpd6fns
      npar <- 3
      nms <- c("lsigma", if (idl) "xi" else "lxi", "lkappa")
      nms2 <- c('logscale', if (idl) 'shape' else 'logshape', 'logkappa')
      attr(family, "type") <- 3
    } else {
      stop("degpd$m must be 1, 2, 3, 4, 5, or 6.")
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
    } else if (zidegpd$m == 2) {
      lik.fns <- .zidegpd2fns
      npar <- 6
      nms <- c("lsigma", "lxi", "lkappa1", "ldkappa", "logitp", "logitpi")
      nms2 <- c('logscale', 'logshape', 'logkappa1', 'logdkappa', 'logitp', 'logitpi')
      attr(family, "type") <- 6
    } else if (zidegpd$m == 3) {
      lik.fns <- .zidegpd3fns
      npar <- 4
      nms <- c("lsigma", "lxi", "ldelta", "logitpi")
      nms2 <- c('logscale', 'logshape', 'logdelta', 'logitpi')
      attr(family, "type") <- 4
    } else if (zidegpd$m == 4) {
      lik.fns <- .zidegpd4fns
      npar <- 5
      nms <- c("lsigma", "lxi", "lkappa", "ldelta", "logitpi")
      nms2 <- c('logscale', 'logshape', 'logkappa', 'logdelta', 'logitpi')
      attr(family, "type") <- 5
    } else if (zidegpd$m == 5) {
      lik.fns <- .zidegpd5fns
      npar <- 4
      nms <- c("lsigma", "lxi", "lkappa", "logitpi")
      nms2 <- c('logscale', 'logshape', 'logkappa', 'logitpi')
      attr(family, "type") <- 2
    } else if (zidegpd$m == 6) {
      lik.fns <- .zidegpd6fns
      npar <- 4
      nms <- c("lsigma", "lxi", "lkappa", "logitpi")
      nms2 <- c('logscale', 'logshape', 'logkappa', 'logitpi')
      attr(family, "type") <- 3
    } else {
      stop("zidegpd$m must be 1, 2, 3, 4, 5, or 6.")
    }
  } else if (family == "pig") {
    ## Poisson-inverse Gaussian (mixed-Poisson count model; not an EGPD).
    ## mu = mean (log link), sigma = dispersion (log link).
    lik.fns <- .pig1fns
    npar <- 2
    nms <- c("lmu", "lsigma")
    nms2 <- c('logmu', 'logsigma')
    attr(family, "type") <- 1
  } else if (family == "zipig") {
    ## Zero-inflated PIG: adds a logit-link zero-inflation probability pi.
    lik.fns <- .zipig1fns
    npar <- 3
    nms <- c("lmu", "lsigma", "logitpi")
    nms2 <- c('logmu', 'logsigma', 'logitpi')
    attr(family, "type") <- 1
  } else if (family == "gpig") {
    ## Generalised PIG (Zhu & Joe 2009), mean parameterisation:
    ## mu = mean (log link), a = tail exponent (logit), c = down-weight (logit).
    lik.fns <- .gpig1fns
    npar <- 3
    nms <- c("lmu", "logita", "logitc")
    nms2 <- c('logmu', 'logita', 'logitc')
    attr(family, "type") <- 1
  } else if (family == "gpignat") {
    ## Generalised PIG, native (a, b, c) parameterisation (paper-faithful):
    ## a (logit), b = level (log), c = down-weight (logit).
    lik.fns <- .gpignat1fns
    npar <- 3
    nms <- c("logita", "lb", "logitc")
    nms2 <- c('logita', 'logb', 'logitc')
    attr(family, "type") <- 1
  } else if (family == "zigpig") {
    ## Zero-inflated GPIG, mean parameterisation: adds logit-link pi.
    lik.fns <- .zigpig1fns
    npar <- 4
    nms <- c("lmu", "logita", "logitc", "logitpi")
    nms2 <- c('logmu', 'logita', 'logitc', 'logitpi')
    attr(family, "type") <- 1
  } else if (family == "zigpignat") {
    ## Zero-inflated GPIG, native parameterisation: adds logit-link pi.
    lik.fns <- .zigpignat1fns
    npar <- 4
    nms <- c("logita", "lb", "logitc", "logitpi")
    nms2 <- c('logita', 'logb', 'logitc', 'logitpi')
    attr(family, "type") <- 1
  } else if (family == "bell") {
    ## Bell (Castellares et al. 2018), mean parameterisation:
    ## mu = mean (log link); internally theta = W0(mu).
    lik.fns <- .bell1fns
    npar <- 1
    nms <- c("lmu")
    nms2 <- c('logmu')
    attr(family, "type") <- 1
  } else if (family == "bellnat") {
    ## Bell, native theta parameterisation (paper Definition 1): theta (log link).
    lik.fns <- .bellnat1fns
    npar <- 1
    nms <- c("ltheta")
    nms2 <- c('logtheta')
    attr(family, "type") <- 1
  } else if (family == "zibell") {
    ## Zero-inflated Bell, mean parameterisation: adds logit-link pi.
    lik.fns <- .zibell1fns
    npar <- 2
    nms <- c("lmu", "logitpi")
    nms2 <- c('logmu', 'logitpi')
    attr(family, "type") <- 1
  } else if (family == "zibellnat") {
    ## Zero-inflated Bell, native parameterisation: adds logit-link pi.
    lik.fns <- .zibellnat1fns
    npar <- 2
    nms <- c("ltheta", "logitpi")
    nms2 <- c('logtheta', 'logitpi')
    attr(family, "type") <- 1
  } else if (family == "cmp") {
    ## Conway-Maxwell-Poisson, mean parameterisation: mu = mean (log), nu = dispersion (log).
    ## nu < 1 overdispersed, nu = 1 Poisson, nu > 1 underdispersed. Tail is lighter than
    ## geometric for any nu > 0, so xi = 0: a dispersion model, not a heavy-tail model.
    lik.fns <- .cmp1fns
    npar <- 2
    nms <- c("lmu", "lnu")
    nms2 <- c('logmu', 'lognu')
    attr(family, "type") <- 1
  } else if (family == "prl") {
    ## Poisson Ramos-Louzada (Alkhairy 2023), mean parameterisation: mu = mean (log link).
    ## One parameter: tau = (mu + sqrt(mu^2-4mu))/2, so mu >= 4 and DI ~ mu are structural.
    lik.fns <- .prl1fns
    npar <- 1
    nms <- c("lmu")
    nms2 <- c('logmu')
    attr(family, "type") <- 1
  } else if (family == "gpois") {
    ## Generalized Poisson (Consul & Jain), mean parameterisation:
    ## mu = mean (log link), lambda = dispersion (logit); Var = mu/(1-lambda)^2.
    lik.fns <- .gpois1fns
    npar <- 2
    nms <- c("lmu", "logitlambda")
    nms2 <- c('logmu', 'logitlambda')
    attr(family, "type") <- 1
  } else if (family == "gwaring") {
    ## Generalized Waring, mean parameterisation: mu = mean (log), k (log), rho (log).
    ## a = mu(rho-1)/k. Tail P(Y=y) ~ y^-(rho+1): the only heavy-tailed member here.
    lik.fns <- .gwaring1fns
    npar <- 3
    nms <- c("lmu", "lk", "lrho")
    nms2 <- c('logmu', 'logk', 'logrho')
    attr(family, "type") <- 1
  } else if (family == "ztgwaring") {
    ## Zero-truncated generalized Waring: the positive block of a hurdle model
    ## (see fit_hurdle). Same parameters as "gwaring"; support is y >= 1.
    lik.fns <- .ztgwaring1fns
    npar <- 3
    nms <- c("lmu", "lk", "lrho")
    nms2 <- c('logmu', 'logk', 'logrho')
    attr(family, "type") <- 1
  } else if (family == "ztdegpd") {
    ## Zero-truncated DEGPD model 1: the positive block of a hurdle model.
    ## Same parameters and link choice as "degpd"; support is y >= 1.
    if (is.null(degpd$m)) degpd$m <- 1
    if (!identical(as.integer(degpd$m), 1L))
      stop("ztdegpd currently supports only degpd.args = list(m = 1).")
    if (!is.null(degpd$link) && identical(degpd$link, "identity")) {
      lik.fns <- .ztdegpd1idfns
      nms <- c("lsigma", "xi", "lkappa")
      nms2 <- c('logscale', 'shape', 'logkappa')
    } else {
      lik.fns <- .ztdegpd1fns
      nms <- c("lsigma", "lxi", "lkappa")
      nms2 <- c('logscale', 'logshape', 'logkappa')
    }
    npar <- 3
    attr(family, "type") <- 1
  } else if (family == "plnorm") {
    ## Poisson-lognormal, mean parameterisation: mu = mean (log), sigma (log).
    ## muz = log(mu) - sigma^2/2; pmf by adaptive Gauss-Hermite quadrature.
    lik.fns <- .plnorm1fns
    npar <- 2
    nms <- c("lmu", "lsigma")
    nms2 <- c('logmu', 'logsigma')
    attr(family, "type") <- 1
  } else if (family == "comppareto") {
    comp.info <- .comppareto_family_setup(comppareto)
    lik.fns <- comp.info$lik.fns
    npar <- comp.info$npar
    nms <- comp.info$nms
    nms2 <- comp.info$nms2
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
      stop(paste("Family '", family, "' not supported. Use 'egpd', 'degpd', 'zidegpd', 'pig', 'gpig', 'bell', 'gpois', 'gwaring', 'ztgwaring', 'ztdegpd', 'plnorm', 'prl', 'cmp', or 'comppareto'.", sep=""))
    }
  }
  out <- list(npar=npar, npar2=npar, lik.fns=lik.fns, nms=nms, family=family, nms2 = nms2)
}
