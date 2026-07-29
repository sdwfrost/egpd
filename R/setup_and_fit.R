## Setup and fitting helper functions
## Adapted from evgam by Ben Youngman

############ .setup.formulae ##########################

.setup.formulae <- function(formula, npar, npar2, data, trace, nms) {
if (inherits(formula, "formula"))
  formula <- lapply(seq_len(npar), function(i) formula)
if (npar == 1) {
  if (!(length(formula) %in% c(npar, 1)))
    stop("length(formula) for this family should be 1")
} else {
  if (!(length(formula) %in% c(npar, 1)))
    stop(paste("length(formula) for this family should be", npar, "(or 1 if all parameters are to have the same formula)"))
}
pred.vars <- unique(unlist(lapply(lapply(formula, mgcv::interpret.gam), "[[", "fake.names")))
off.vars <- pred.vars[grepl("^offset\\(", pred.vars)]
pred.vars <- pred.vars[!grepl("^offset\\(", pred.vars)]
if (length(off.vars) > 0)
  pred.vars <- unique(c(pred.vars, all.vars(reformulate(off.vars))))
if (!all(pred.vars %in% names(data))) {
  missing.vars <- pred.vars[!(pred.vars %in% names(data))]
  stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `data'.", sep=""))
}
terms.list <- lapply(formula, terms.formula, specials=c("s", "te", "ti"))
got.specials <- sapply(lapply(terms.list, function(x) unlist(attr(x, "specials"))), any)
termlabels.list <- lapply(terms.list, attr, "term.labels")
got.intercept <- sapply(terms.list, attr, "intercept") == 1
## Per-element offset term strings (e.g. "offset(z)"). terms.formula() strips
## offsets out of term.labels, so the reformulate() calls below used to drop
## them silently for every non-first parameter -- meaning `lxi = ~ offset(z)`
## had no effect on the fit. We capture them here and reattach them below.
offterms.list <- lapply(terms.list, function(tt) {
  oi <- attr(tt, "offset")
  if (is.null(oi)) return(character(0))
  vars <- attr(tt, "variables")
  vapply(oi, function(k) deparse(vars[[k + 1L]]), character(1))
})
for (i in seq_along(termlabels.list)) {
  if (length(termlabels.list[[i]]) == 0) {
    if (got.intercept[i]) {
      termlabels.list[[i]] <- "1"
    } else if (length(offterms.list[[i]]) > 0) {
      ## offset-only predictor, e.g. ~ offset(z) - 1: a fixed/known linear
      ## predictor with no free coefficients (0-column design). The parameter
      ## is pinned to its offset. Leave term labels empty; the offset is
      ## reattached and the dropped intercept is respected below.
    } else {
      stop(paste("formula element", i, "incorrectly specified"))
    }
  }
}
got.response <- sapply(terms.list, attr, "response") == 1
if (!got.response[1]) {
  stop("formula has no response")
} else {
  response.name <- as.character(formula[[1]])[2]
}
if (any(!got.response)) {
  for (i in which(!got.response)) {
    formula[[i]] <- reformulate(termlabels=c(termlabels.list[[i]], offterms.list[[i]]),
                                response=response.name, intercept=got.intercept[i])
  }
}
stripped.formula <- lapply(seq_along(termlabels.list), function(i)
  reformulate(termlabels=c(termlabels.list[[i]], offterms.list[[i]]),
              intercept=got.intercept[i]))
censored <- FALSE
if (substr(response.name, 1, 5) == "cens(") {
  response.name <- substr(response.name, 6, nchar(response.name) - 1)
  response.name <- gsub(" ", "", response.name)
  response.name <- strsplit(response.name, ",")[[1]]
  if (length(response.name) > 2)
    stop("Censored response can only contain two variables.")
  rr <- response.name[2]
  formula <- lapply(termlabels.list, function(x) reformulate(termlabels=x, response=rr))
  censored <- TRUE
}
attr(formula, "response.name") <- response.name
pred.vars <- pred.vars[!(pred.vars %in% response.name)]
attr(formula, "predictor.names") <- pred.vars
attr(formula, "stripped") <- stripped.formula
attr(formula, "censored") <- censored
attr(formula, "smooths") <- got.specials
attr(formula, "npar") <- npar
for (i in seq_along(formula)) {
  attr(formula[[i]], "intercept") <- got.intercept[i]
  attr(formula[[i]], "smooth") <- got.specials[i]
}
names(formula) <- nms
formula
}

############ .predictable.gam ##########################

.predictable.gam <- function(G, formula) {
keep <- c("dev.extra", "pterms", "nsdf", "X", "terms", "mf", "smooth", "sp", "term.names", "offset")
G <- G[keep]
## Offset-only predictor (e.g. ~ offset(z) - 1): mgcv returns a 0-column design,
## which would create zero-width Hessian blocks (Armadillo submat out-of-bounds in
## gradHess). Inject a single all-zero column so the parameter is fixed to its
## offset (X %*% beta == 0) without any C++ change; its coefficient is inert (zero
## gradient/Hessian) and is absorbed by the existing pinv/perturb path.
if (ncol(G$X) == 0) {
  G$X <- matrix(0, nrow = nrow(G$X), ncol = 1)
  colnames(G$X) <- "(fixed)"
  G$term.names <- "(fixed)"
}
G$nb <- ncol(G$X)
G$coefficients <- numeric(G$nb)
old <- c("mf", "pP", "cl")
new <- c("model", "paraPen", "call")
is.in <- !is.na(match(old, names(G)))
if (any(is.in)) names(G)[match(old[is.in], names(G))] <- new[is.in]
G$formula <- formula
class(G) <- "gamlist"
G
}

############ .X.egpd ##########################

.X.egpd <- function(object, newdata) {
object <- object[sapply(object, inherits, what="gamlist")]
if (missing(newdata)) {
  X <- lapply(object, function(x) x$X)
  offsets <- lapply(object, function(x) if (is.null(x$offset)) numeric(0) else x$offset)
} else {
  for (i in seq_along(object))
    class(object[[i]]) <- "gam"
  X <- lapply(object, mgcv::predict.gam, newdata=newdata, type="lpmatrix")
  ## mirror the offset-only zero-column injection from .predictable.gam so that
  ## newdata prediction matches the fitted design (preserve the model offset).
  X <- lapply(X, function(x) {
    if (ncol(x) == 0) {
      off <- attr(x, "model.offset")
      z <- matrix(0, nrow = nrow(x), ncol = 1)
      if (!is.null(off)) attr(z, "model.offset") <- off
      z
    } else x
  })
  offsets <- lapply(X, function(x) {
    off <- attr(x, "model.offset")
    if (is.null(off)) numeric(0) else off
  })
}
for (i in seq_along(X))
    colnames(X[[i]]) <- object[[i]]$term.names
names(X) <- names(object)
attr(X, "offsets") <- offsets
X
}

############ .setup.data ##########################

.setup.data <- function(data, responsename, formula, family, nms, removeData,
knots, maxdata, maxspline, compact, sargs, outer, trace, gamma) {

## data
for (i in seq_along(responsename)) {
  dm <- as.matrix(data[,responsename[i]])
  data <- data[rowSums(!is.na(dm)) == ncol(dm), , drop = FALSE]
}

if (nrow(data) > maxdata) {
    id <- sort(sample(nrow(data), maxdata))
    data <- data[id,]
    if (trace >= 0)
      message("`data' truncated to `maxdata' rows. Re-supply `data' to, e.g., `predict.egpd'")
}

if  (compact) {
data.undup <- as.list(data[,unique(unlist(lapply(formula, function(y) unlist(lapply(mgcv::interpret.gam(y)$smooth.spec, function(x) x$term))))), drop=FALSE])
data.undup <- lapply(data.undup, function(x) as.integer(as.factor(x)))
if (length(data.undup) > 1) for (i in 2:length(data.undup)) data.undup[[1]] <- paste(data.undup[[1]], data.undup[[i]], sep=":")
data.undup <- data.undup[[1]]
gc()
unq.id <- which(!duplicated(data.undup))
data.unq <- data.undup[unq.id]
dup.id <- match(data.undup, data.unq)
}

subsampling <- FALSE

## gams
gams <- list()

for (i in seq_along(formula)) {
  if (nrow(data) > maxspline) {
    id <- sample(nrow(data), maxspline)
    gams[[i]] <- mgcv::gam(formula[[i]], data=data[id,], fit=FALSE, knots=knots, method="REML")
  } else {
    gams[[i]] <- mgcv::gam(formula[[i]], data=data, fit=FALSE, knots=knots, method="REML")
  }
  gams[[i]] <- .predictable.gam(gams[[i]], formula[[i]])
}

gc()

## likelihood
lik.data <- list()
lik.data$control <- list()
lik.data$outer <- outer
lik.data$control$outer <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-2, stepmax=3)
lik.data$control$inner <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-4, stepmax=1e2)
lik.data$y <- as.matrix(data[,responsename, drop=FALSE])
lik.data$Mp <- sum(unlist(sapply(gams, function(y) c(1, sapply(y$smooth, function(x) x$null.space.dim)))))
lik.data$const <- .5 * lik.data$Mp * log(2 * pi)
lik.data$nobs <- nrow(lik.data$y)
if (attr(formula, "censored")) {
  lik.data$censored <- TRUE
  if (any(lik.data$y[,2] < lik.data$y[,1]))
    stop("For censored response need right >= left in `cens(left, right)'")
  lik.data$cens.id <- lik.data$y[,2] > lik.data$y[,1]
  if (trace >= 0 & sum(lik.data$cens.id) == 0) {
    message("No response data appear to be censored. Switching to uncensored likelihood.")
    lik.data$censored <- FALSE
  }
} else {
  lik.data$censored <- FALSE
}
lik.data$sandwich <- !is.null(sargs$id)
if (lik.data$sandwich)
  lik.data$sandwich.split <- data[,sargs$id]

if (!compact) {
  if (nrow(data) > maxspline) {
    lik.data$X <- .X.egpd(gams, data)
  } else {
    lik.data$X <- .X.egpd(gams)
  }
  lik.data$dupid <- 0
  lik.data$duplicate <- 0
  lik.data$offsets <- attr(lik.data$X, "offsets")
} else {
  lik.data$X <- .X.egpd(gams, data[unq.id,])
  lik.data$dupid <- dup.id - 1
  lik.data$duplicate <- 1
  ## offsets are per-observation (full-length), expand from unique rows
  lik.data$offsets <- lapply(attr(lik.data$X, "offsets"), function(off) {
    if (length(off) == 0) {
      off
    } else if (length(off) == 1) {
      rep(off, length(dup.id))
    } else {
      off[dup.id]
    }
  })
}
for (i in seq_along(gams)) {
  if (removeData)
    gams[[i]]$y <- NULL
}
if (length(lik.data$X) == 1 & length(nms) > 1) {
  for (i in 2:length(nms)) {
    lik.data$X[[i]] <- lik.data$X[[1]]
    gams[[i]] <- gams[[1]]
    lik.data$offsets[[i]] <- lik.data$offsets[[1]]
  }
}
nbk <- sapply(lik.data$X, ncol)
lik.data$nb <- sum(nbk)
lik.data$idpars <- rep(seq_along(lik.data$X), nbk)
lik.data$LAid <- lik.data$idpars > 0
lik.data$subsampling <- subsampling
gotsmooth <- which(sapply(gams, function(x) length(x$sp)) > 0)
lik.data$k <- 1 / gamma
if (is.null(sargs$id)) {
  lik.data$adjust <- 0
} else {
  if (gamma != 1)
    stop("Can't have gamma != 1 and sandwich adjustment.")
  if (is.null(sargs$method))
    sargs$method <- "magnitude"
  if (sargs$method == "curvature") {
    if (trace > 0)
      message(paste("Sandwich adjustment method: curvature"))
    lik.data$adjust <- 2
  } else {
    if (trace > 0)
      message(paste("Sandwich adjustment method: magnitude"))
    lik.data$adjust <- 1
  }
}
if (is.null(sargs$force))
  sargs$force <- FALSE
lik.data$force <- sargs$force
lik.data$npar <- attr(formula, "npar")
list(lik.data=lik.data, gotsmooth=gotsmooth, data=data, gams=gams, sandwich=lik.data$adjust > 0)
}

############ .setup.inner.inits ##########################

.setup.inner.inits <- function(inits, likdata, likfns, npar, family) {

likdata0 <- likdata
likdata0$X <- lapply(seq_along(likdata$X), function(i) matrix(1, nrow=nrow(likdata$y), ncol=1))
likdata0$dupid <- 0
likdata0$duplicate <- 0
likdata0$S <- diag(0, npar)
likdata0$idpars <- seq_len(npar)

if (is.null(inits)) {
  if (!is.null(likfns$inits)) {
    inits <- likfns$inits(likdata)
  } else {
  if (family %in% c("egpd", "degpd")) {
    inits <- numeric(npar)
    ## The second linear predictor is lxi under the default log link but xi itself
    ## under the identity link (.degpd1idfns, flagged by identity.xi). Seeding both
    ## with 0.05 therefore starts the identity-link fit from xi = 0.05 rather than
    ## xi = exp(0.05) = 1.051 -- a tail 20x lighter. On series whose largest counts
    ## sit far out relative to sigma = mean(y), the DEGPD log-density then underflows
    ## at those counts, log p = -Inf, and the fit aborts with "Some gradient
    ## non-finite" before the optimiser can take a step. Seed the same xi under both.
    xi0 <- if (isTRUE(likfns$identity.xi)) exp(.05) else .05
    inits[1:2] <- c(log(mean(likdata$y[,1])), xi0)
    if (attr(family, "type") == 6)
      inits <- c(inits[1:2], 0, 0, 0)
  } else {
    if (family == "zidegpd") {
      inits <- numeric(npar)
      inits[1:2] <- c(log(mean(likdata$y[,1])), .05)
      if (attr(family, "type") == 6)
        inits <- c(inits[1:2], 0, 0, 0, 0)
    } else if (family %in% c("pig", "zipig")) {
      ## PIG/ZIPIG: seed log(mu) at the log mean and log(sigma) at a moment
      ## estimate of the dispersion (var = mu + sigma*mu^2). zipig's logitpi
      ## is seeded near the observed zero proportion.
      y1 <- likdata$y[, 1]
      inits <- numeric(npar)
      inits[1] <- log(max(mean(y1), 1e-3))
      s0 <- (stats::var(y1) - mean(y1)) / max(mean(y1)^2, 1e-8)
      inits[2] <- log(min(max(s0, 0.1), 10))
      if (family == "zipig")
        inits[3] <- stats::qlogis(min(max(mean(y1 == 0), 1e-3), 0.5))
    } else if (family %in% c("gpig", "gpignat", "zigpig", "zigpignat")) {
      ## GPIG/ZIGPIG (Zhu & Joe 2009). Seed a = 1/2 (the PIG special case),
      ## derive c from the empirical dispersion index D via D-1=(1-a)c/(1-c),
      ## and mu (or the implied b) from the sample mean.
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3)
      vv <- stats::var(y1)
      D <- if (is.finite(vv) && ybar > 0) max(vv / ybar, 1.05) else 2
      a0 <- 0.5
      cc <- min(max((D - 1) / ((1 - a0) + (D - 1)), 0.05), 0.95)
      inits <- numeric(npar)
      if (family %in% c("gpig", "zigpig")) {        # mean: logmu, logita, logitc
        inits[1] <- log(ybar)
        inits[2] <- stats::qlogis(a0)
        inits[3] <- stats::qlogis(cc)
      } else {                                      # native: logita, logb, logitc
        b0 <- ybar * (1 - cc)^(1 - a0) / (a0 * cc)
        inits[1] <- stats::qlogis(a0)
        inits[2] <- log(max(b0, 1e-3))
        inits[3] <- stats::qlogis(cc)
      }
      if (family %in% c("zigpig", "zigpignat"))
        inits[4] <- stats::qlogis(min(max(mean(y1 == 0), 1e-3), 0.5))
    } else if (family == "cmp") {
      ## CMP: Var ~ mu/nu asymptotically, so seed nu at 1/DI (clamped into the valid range).
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3); vv <- stats::var(y1)
      DI <- if (is.finite(vv) && vv > 0) max(vv / ybar, 1e-3) else 1
      inits <- c(log(ybar), log(min(max(1 / DI, 1e-5), 50)))
    } else if (family == "prl") {
      ## PRL cannot represent a mean below 4 (mu = tau^2/(tau-1) is minimised at tau = 2).
      y1 <- likdata$y[, 1]
      inits <- log(max(mean(y1), 4 + 1e-3))
    } else if (family == "gpois") {
      ## Generalized Poisson: Var/mean = 1/(1-lambda)^2, so lambda = 1 - sqrt(mean/var).
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3); vv <- stats::var(y1)
      lam0 <- if (is.finite(vv) && vv > ybar) 1 - sqrt(ybar / vv) else 0.1
      lam0 <- min(max(lam0, 0.01), 0.98)
      inits <- c(log(ybar), stats::qlogis(lam0))
    } else if (family == "gwaring") {
      ## Generalized Waring: seed the mean at ybar, k = 1 (the Waring special case), and
      ## rho = 2.5, just inside the finite-variance region rho > 2.
      y1 <- likdata$y[, 1]
      inits <- c(log(max(mean(y1), 1e-3)), 0, log(2.5))
    } else if (family == "plnorm") {
      ## Poisson-lognormal: matching the first two moments gives
      ## sigma^2 = log(1 + (var - mean)/mean^2).
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3); vv <- stats::var(y1)
      s2 <- if (is.finite(vv) && vv > ybar) log(1 + (vv - ybar) / ybar^2) else 0.25
      inits <- c(log(ybar), 0.5 * log(min(max(s2, 1e-4), 25)))
    } else if (family %in% c("bell", "bellnat", "zibell", "zibellnat")) {
      ## Bell/ZIBell (Castellares et al. 2018). Native theta is seeded at the
      ## MLE theta.hat = W0(ybar); the mean parameterisation at log(ybar).
      ## ZI variants seed logitpi near the observed zero proportion.
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3)
      inits <- numeric(npar)
      if (family %in% c("bell", "zibell")) {          # mean: log mu
        inits[1] <- log(ybar)
      } else {                                        # native: log theta
        inits[1] <- log(max(bell_W0_cpp(ybar), 1e-3))
      }
      if (family %in% c("zibell", "zibellnat"))
        inits[2] <- stats::qlogis(min(max(mean(y1 == 0), 1e-3), 0.5))
    } else {
      inits <- numeric(npar)
    }
  }
  }
  likdata0$CH <- diag(length(inits))
  likdata0$compmode <- numeric(length(inits))
  beta0 <- .newton_step_inner(inits, .nllh.nopen, .search.nopen, likdata=likdata0, likfns=likfns, control=likdata$control$inner)$par
} else {
  if (is.list(inits)) {
    betamat <- expand.grid(inits)
    betanllh <- numeric(nrow(betamat))
    for (i in seq_len(nrow(betamat))) {
      beta0 <- unlist(betamat[i,])
      betanllh[i] <- likfns$d0(beta0, likdata0)
    }
    beta0 <- betamat[which.min(betanllh),]
  } else {
    beta0 <- inits
  }
}
beta0 <- unlist(lapply(seq_len(npar), function(i) {
  nci <- ncol(likdata$X[[i]])
  if (nci == 0) numeric(0) else c(beta0[i], rep(0, nci - 1))  # 0-column (offset-only) block contributes no coefficients
}))
compmode <- 0 * beta0
CH <- diag(compmode + 1)
k <- likdata$k
likdata[c("k", "CH", "compmode")] <- list(k, CH, compmode)
diagH <- diag(.gH.nopen(beta0, likdata=likdata, likfns=likfns)[[2]])
if (likdata$sandwich) {
  beta0 <- .newton_step(beta0, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)$par
  H <- .gH.nopen(beta0, likdata=likdata, likfns=likfns, sandwich=TRUE)
  J <- split(as.data.frame(t(H[[1]])), likdata$sandwich.split)
  J <- sapply(J, colSums)
  J <- tcrossprod(J)
  H <- H[[2]]
  diagH <- diag(H)
  cholH <- try(chol(H), silent=TRUE)
  if (inherits(cholH, "try-error")) {
    if (!likdata$force) {
      stop("Hessian of unpenalised MLE not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
    } else {
      if (trace >= 0)
        message("Hessian perturbed to be positive definite for sandwich adjustment.")
      iH <- pinv(H)
    }
  } else {
    iH <- chol2inv(cholH)
  }
  if (likdata$adjust == 2) {
    cholJ <- try(chol(J), silent=TRUE)
    if (!inherits(cholJ, "try-error")) {
      HA <- crossprod(backsolve(cholJ, H, transpose=TRUE))
    } else {
      iHA <- tcrossprod(crossprod(iH, J), iH)
      choliHA <- try(chol(iHA), silent=TRUE)
      if (inherits(choliHA, "try-error")) {
        if (!likdata$force) {
          stop("Sandwich variance not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
        } else {
          if (trace >= 0)
            message("Sandwich variance perturbed to be positive definite.")
          HA <- pinv(iHA)
        }
      } else {
        HA <- chol2inv(choliHA)
      }
    }
    sH <- svd(H)
    M <- sqrt(sH$d) * t(sH$v)
    sHA <- svd(HA)
    MA <- sqrt(sHA$d) * t(sHA$v)
    CH <- solve(M, MA)
    compmode <- beta0
  } else {
    k <- 1 / mean(diag(crossprod(iH, J)))
  }
}
attr(beta0, "k") <- k
attr(beta0, "CH") <- CH
attr(beta0, "compmode") <- compmode
attr(beta0, "diagH") <- diagH
beta0
}

############ .guess ##########################

.guess <- function(x, d, s) {
okay <- s != 0
val <- d / (d + exp(x) * s)
mean(val[okay]) - .4
}

############ .sandwich ##########################

.sandwich <- function(likdata, beta) {
likdata$k <- attr(beta, "k")
likdata$CH <- attr(beta, "CH")
likdata$compmode <- attr(beta, "compmode")
bigX <- do.call(cbind, likdata$X)
CHX <- bigX %*% likdata$CH
CHX <- lapply(unique(likdata$idpars), function(i) CHX[,likdata$idpars == i])
likdata$CHX <- CHX
likdata
}

############ .outer ##########################

.outer <- function(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace) {

attr(rho0, "beta") <- beta

if (outer == "fixed") {

  fit.reml <- .reml0_fixed(rho0, likfns=likfns, likdata=likdata, Sdata=Sdata)

} else {

if (is.null(likfns$d340) & outer != "fd")
  outer <- "fd"

if (outer == "newton") {
  fit.reml <- .newton_step_inner(rho0, .reml0, .search.reml, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
} else {
  if (outer == "fd") {
    fit.reml <- .BFGS(rho0, .reml0, .reml1.fd, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  } else {
    fit.reml <- .BFGS(rho0, .reml0, .reml1, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  }
  rho1 <- fit.reml$par
  attr(rho1, "beta") <- fit.reml$beta
  if (correctV) {
    fit.reml$Hessian <- try(.reml12(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)[[2]], silent=TRUE)
    if (inherits(fit.reml$Hessian, "try-error"))
      fit.reml$Hessian <- try(.reml2.fd(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata), silent=TRUE)
    if (inherits(fit.reml$Hessian, "try-error"))
      fit.reml$Hessian <- .reml2.fdfd(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)
  }
}

if (correctV)
  fit.reml$invHessian <- .solve_egpd(fit.reml$Hessian)

fit.reml$trace <- trace

if (trace == 1) {
  report <- "\n Final max(|grad|))"
  likdata$S <- .makeS(Sdata, exp(fit.reml$par))
  report <- c(report, paste("   Inner:", signif(max(abs(.gH.pen(fit.reml$beta, likdata, likfns)[[1]])), 3)))
  report <- c(report, paste("   Outer:", signif(max(abs(fit.reml$gradient)), 3)))
  report <- c(report, "", "")
  cat(paste(report, collapse="\n"))
}

}

fit.reml

}

############ .outer.nosmooth ##########################

.outer.nosmooth <- function(beta, likfns, likdata, control, trace) {

fit.inner <- .newton_step(beta, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)

fit.inner$beta <- fit.inner$par
fit.inner

}

############ .VpVc ##########################

.VpVc <- function(fitreml, likfns, likdata, Sdata, correctV, sandwich, smooths, trace) {
lsp <- fitreml$par
H0 <- .gH.nopen(fitreml$beta, likdata, likfns)[[2]]
if (smooths) {
  sp <- exp(lsp)
  H <- H0 + likdata$S
} else {
  H <- H0
}
## Coefficients pinned by an offset-only predictor (`fixed=`, see .predictable.gam)
## carry an all-zero design column, so their row and column of H are exactly zero and
## H is singular. That direction is a nuisance -- orthogonal to every estimated
## coefficient -- but left alone it corrupts the whole covariance matrix through
## pinv() and the smoothing-parameter correction below; on one test fit it produced a
## negative variance for the shape parameter, so confint() returned NA. Giving those
## directions a unit diagonal makes H block-diagonal and invertible without altering
## the inverse of the free block; their rows and columns are zeroed out afterwards,
## which also leaves them exactly zero effective degrees of freedom in .edf().
inert <- which(rowSums(abs(H)) == 0 & colSums(abs(H)) == 0)
if (length(inert)) diag(H)[inert] <- 1

cholH <- try(chol(H), silent=TRUE)
if (inherits(cholH, "try-error") & trace >= 0)
  message("Final Hessian of negative penalized log-likelihood not numerically positive definite.")
Vc <- Vp <- pinv(H)
if (smooths) {
if (correctV) {
## Vc corrects Vp for uncertainty in the smoothing parameters. That correction only
## exists if they were estimated: when the caller supplies `sp` they are fixed
## quantities carrying no uncertainty, fitreml then has no invHessian to propagate,
## and the correction term is exactly zero, so Vc = Vp. Previously this branch ran
## regardless and failed on the missing invHessian ("requires numeric/complex matrix
## arguments"), which made any `sp = ` fit unusable at the default correctV = TRUE.
Vrho <- fitreml$invHessian
if (is.null(Vrho) || !is.matrix(Vrho) || nrow(Vrho) != length(lsp)) {
  Vrho <- 0
} else {
  cholVp <- try(chol(Vp), silent=TRUE)
  if (inherits(cholVp, "try-error")) {
      cholVp <- attr(.perturb(Vp), "chol")
  }
  attr(lsp, "beta") <- fitreml$beta
  spSl <- Map("*", attr(Sdata, "Sl"), exp(lsp))
  dbeta <- .d1beta(lsp, fitreml$beta, spSl, .Hdata(H))$d1
  Vbetarho <- tcrossprod(dbeta %*% Vrho, dbeta)
  VR <- matrix(0, nrow=likdata$nb, ncol=likdata$nb)
  Vc <- .perturb(Vp + Vbetarho + VR)
}
} else {
  Vrho <- 0
}
} else {
  Vrho <- 0
}
## a pinned coefficient has no sampling variability and no covariance with anything
if (length(inert)) {
  Vp[inert, ] <- 0; Vp[, inert] <- 0
  Vc[inert, ] <- 0; Vc[, inert] <- 0
}
list(Vp=Vp, Vc=Vc, Vlsp=Vrho, H0=H0, H=H)
}

############ .edf ##########################

.edf <- function(beta, likfns, likdata, VpVc, sandwich) {
diag(crossprod(VpVc$Vp, VpVc$H0))
}

############ .swap ##########################

.swap <- function(fitreml, gams, likdata, VpVc, gotsmooth, edf, smooths) {
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) {
  spl <- split(exp(fitreml$par), unlist(sapply(seq_along(gams), function(x) rep(x, length(gams[[x]]$sp)))))
  sp <- replace(lapply(seq_along(gams), function(x) NULL), gotsmooth, spl)
}
for (i in seq_along(gams)) {
  idi <- likdata$idpars == i
  gams[[i]]$coefficients <- fitreml$beta[idi]
  names(gams[[i]]$coefficients) <- gams[[i]]$term.names
  gams[[i]]$Vp <- Vp[idi, idi, drop = FALSE]
  gams[[i]]$Vc <- Vc[idi, idi, drop = FALSE]
  if (i %in% gotsmooth) gams[[i]]$sp <- sp[[i]]
  gams[[i]]$edf <- edf[idi]
}
gams
}

############ .finalise ##########################

.finalise <- function(gams, data, likfns, likdata, Sdata, fitreml, VpVc, family, gotsmooth,
formula, responsenm, removeData, edf, linkNames) {
names(gams) <- linkNames
smooths <- length(gotsmooth) > 0
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) gams$sp <- exp(fitreml$par)
gams$nobs <- likdata$nobs
gams$convergence <- if (is.null(fitreml$convergence)) NA_integer_ else as.integer(fitreml$convergence)
gams$logLik <- NA_real_
## NB: optimizers return a *double* 0 on success, so the old
## identical(gams$convergence, 0L) test was always FALSE and logLik/AIC/BIC
## were never populated. Compare numerically instead.
converged <- isTRUE(gams$convergence == 0L)
gams$df <- sum(edf)
if (converged) {
## The unpenalised log-likelihood at the fitted coefficients, and nothing else.
## This used to have likdata$const = .5 * Mp * log(2 * pi) subtracted from it,
## which is the REML criterion's normalising constant and has no place in a
## log-likelihood: it is not a Laplace correction (that would *add* the term and
## carry .5*log|S|+ - .5*log|H| with it), and .nllh.nopen is already used bare
## as the log-likelihood by the profile intervals in confint(). Because Mp counts
## one per linear predictor plus each smooth's null-space dimension, the offset
## grew with model complexity, so it biased AIC/BIC against the more flexible
## model rather than cancelling in a comparison.
gams$logLik <- -.nllh.nopen(fitreml$beta, likdata, likfns)
gams$AIC <- -2 * gams$logLik + 2 * gams$df
gams$BIC <- -2 * gams$logLik + log(likdata$nobs) * gams$df
} else {
gams$AIC <- gams$BIC <- 1e20
}
attr(gams, "df") <- sum(edf)
gams$simulate <- list(mu=fitreml$beta, Sigma=Vp)
gams$family <- family
gams$idpars <- likdata$idpars
nms <- names(gams)[seq_along(formula)]
if (family != "custom") {
  logits <- substr(nms, 1, 5) == "logit"
  if (any(logits))
    nms[logits] <- gsub("logit", "", nms[logits])
  logs <- substr(nms, 1, 3) == "log"
  if (any(logs))
    nms[logs] <- gsub("log", "", nms[logs])
  probits <- substr(nms, 1, 6) == "probit"
  if (any(probits))
    nms[probits] <- gsub("probit", "", nms[probits])
}
gams$predictor.names <- attr(formula, "predictor.names")
formula <- attr(formula, "stripped")
names(formula) <- nms
gams$call <- formula
gams$response.name <- responsenm
gams$gotsmooth <- gotsmooth
if (!removeData) {
    gams$data <- data
}
gams$Vc <- Vc
gams$Vp <- Vp
gams$Vlsp <- VpVc$Vlsp
gams$negREML <- if (is.null(fitreml$objective)) NA_real_ else fitreml$objective
gams$coefficients <- as.vector(fitreml$beta)
for (i in seq_along(likdata$X)) {
  gams[[i]]$X <- likdata$X[[i]]
  if (likdata$duplicate == 1)
    gams[[i]]$X <- gams[[i]]$X[likdata$dupid + 1, , drop = FALSE]
  gams[[i]]$fitted <- as.vector(gams[[i]]$X %*% gams[[i]]$coefficients)
  if (length(likdata$offsets[[i]]) > 0)
    gams[[i]]$fitted <- gams[[i]]$fitted + likdata$offsets[[i]]
  names(gams[[i]]$coefficients) <- colnames(gams[[i]]$X)
}
gams$likdata <- likdata
gams$likfns <- likfns
if (smooths) gams$Sdata <- Sdata
gams$formula <- formula
gams$compacted <- likdata$duplicate == 1
if (gams$compacted) gams$compactid <- likdata$dupid + 1
smooth.terms <- unique(lapply(lapply(gams[gotsmooth], function(x) x$smooth), function(y) lapply(y, function(z) z$term)))
smooth.terms <- unique(unlist(smooth.terms, recursive=FALSE))
gams$plotdata <- lapply(smooth.terms, function(x) unique(data[,x, drop=FALSE]))
names(gams$coefficients) <- unlist(lapply(seq_along(likdata$X), function(i) paste(names(gams)[i], names(gams[[i]]$coefficients), sep = "_")))
gams$ngam <- length(formula)
for (i in seq_along(gams[nms])[-gotsmooth])
  gams[[i]]$smooth <- NULL
class(gams) <- "egpd"
return(gams)
}

############ .setup.fixed.egpd ##########################

## Sugar for holding distribution parameters at known values.
##
## `fixed` is a named list of RESPONSE-scale values, e.g. list(kappa = 1). Each is
## translated into the offset-only predictor `~ offset(<col>) - 1`, which pins that
## parameter without any free coefficient (see .predictable.gam). The translation
## exists because the offset has to be supplied on the LINK scale, which is the part
## users get wrong: parameters whose internal name starts with "l" (lsigma, lkappa,
## lxi, ...) use a log link, so the offset is log(value); the rest (e.g. `xi` under
## degpd.args link = "identity") are pinned on their natural scale.
##
## Names are matched against the family's parameter names both with and without the
## leading link letter, so fixed = list(kappa = 1) and list(lkappa = 1) both work and
## both mean "kappa = 1".
.setup.fixed.egpd <- function(fixed, formula, data, nms) {
  if (is.null(names(fixed)) || any(!nzchar(names(fixed))))
    stop("`fixed` must be a named list, e.g. fixed = list(kappa = 1)", call. = FALSE)
  if (is.null(nms))
    stop("this family does not expose named parameters, so `fixed` cannot be used",
         call. = FALSE)
  if (!is.list(formula))
    stop("`fixed` requires the multi-parameter list form of `formula`", call. = FALSE)

  bare <- sub("^l", "", nms)          # lsigma -> sigma, lkappa -> kappa, xi -> xi
  n    <- nrow(data)
  offsets <- list(); values <- list()

  for (nm in names(fixed)) {
    v <- fixed[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      stop(sprintf("fixed$%s must be a single finite number", nm), call. = FALSE)

    i <- match(nm, nms)
    if (is.na(i)) i <- match(nm, bare)
    if (is.na(i))
      stop(sprintf("'%s' is not a parameter of this family; available: %s",
                   nm, paste(nms, collapse = ", ")), call. = FALSE)

    logscale <- startsWith(nms[i], "l")
    if (logscale && v <= 0)
      stop(sprintf("fixed$%s must be positive (it is fitted on the log scale)", nm),
           call. = FALSE)
    eta <- if (logscale) log(v) else v

    col <- paste0(".fixedpar_", nms[i])
    data[[col]] <- rep(eta, n)
    formula[[i]] <- stats::as.formula(paste0("~ offset(", col, ") - 1"))
    offsets[[col]] <- eta            # link scale, for predict() on newdata
    values[[nms[i]]] <- v            # response scale, for the user
  }

  list(formula = formula, data = data, offsets = offsets, values = values)
}
