## Internal helpers: Newton step functions, p-value tests, summary helpers
## Adapted from evgam by Ben Youngman
## P-value testing code from mgcv by Simon Wood

## p-values for smooths
##
## To avoid use of mgcv::: this code is directly taken from mgcv 1.8-14
## by Professor Simon Wood
##
## Thank you, Simon.
##
.smoothTest <- function(b,X,V,eps=.Machine$double.eps^.5) {
  qrx <- qr(X)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  k <- length(ed$values)
  f <- t(ed$vectors[,1:k])%*%R%*%b
  t <- sum(f^2)
  k <- ncol(X)
  lambda <- as.numeric(ed$values[1:k])
  pval <- .liu2(t,lambda)
  list(stat=t,pval=pval)
}

.liu2 <- function(x, lambda, h = rep(1,length(lambda)),lower.tail=FALSE) {
  if (length(h) != length(lambda)) stop("lambda and h should have the same length!")
  lh <- lambda*h
  muQ <- sum(lh)
  lh <- lh*lambda
  c2 <- sum(lh)
  lh <- lh*lambda
  c3 <- sum(lh)
  s1 <- c3/c2^1.5
  s2 <- sum(lh*lambda)/c2^2
  sigQ <- sqrt(2*c2)
  t <- (x-muQ)/sigQ
  if (s1^2>s2) {
    a <- 1/(s1-sqrt(s1^2-s2))
    delta <- s1*a^3-a^2
    l <- a^2-2*delta
  } else {
    a <- 1/s1
    delta <- 0
    l <- c2^3/c3^2
  }
  muX <- l+delta
  sigX <- sqrt(2)*a
  return(pchisq(t*sigX+muX,df=l,ncp=delta,lower.tail=lower.tail))
}

.simf <- function(x,a,df,nq=50) {
  p <- (1:nq-.5)/nq
  q <- qchisq(p,df)
  x <- x*q/df
  pr <- sum(.liu2(x,a))
  pr/nq
}

.testStat <- function(p,X,V,rank=NULL,type=0,res.df= -1) {
  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  siv <- sign(ed$vectors[1,]);siv[siv==0] <- 1
  ed$vectors <- sweep(ed$vectors,2,siv,"*")
  k <- max(0,floor(rank))
  nu <- abs(rank - k)
  if (type < -.5) {
    res <- .smoothTest(p,X,V)
    res$rank <- rank
    return(res)
  } else  if (type==1) {
    if (rank > k + .05||k==0) k <- k + 1
    nu <- 0;rank <- k
  } else if (type==2) {
    nu <- 0;rank <- k <- max(1,round(rank))
    warning("p-values may give low power in some circumstances")
  } else if (type==3) {
    nu <- 0; rank <- k <- max(1,ceiling(rank))
    warning("p-values un-reliable")
  } else if (type==4) {
    rank <- k <- max(sum(ed$values>1e-3*max(ed$values)),1)
    nu <- 0
    warning("p-values may give very low power")
  }
  if (nu>0) k1 <- k+1 else k1 <- k
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]
  if (nu>0&&k>0) {
     if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
     b12 <- .5*nu*(1-nu)
     if (b12<0) b12 <- 0
     b12 <- sqrt(b12)
     B <- matrix(c(1,b12,b12,nu),2,2)
     ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
     B <- ev%*%B%*%ev
     eb <- eigen(B,symmetric=TRUE)
     rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
     vec1 <- vec
     vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
     vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    vec1 <- vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else
            t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
  d <- t(vec)%*%(R%*%p)
  d <- sum(d^2)
  d1 <- t(vec1)%*%(R%*%p)
  d1 <- sum(d1^2)
  rank1 <- rank
  if (nu>0) {
     if (k1==1) rank1 <- val <- 1 else {
       val <- rep(1,k1)
       rp <- nu+1
       val[k] <- (rp + sqrt(rp*(2-rp)))/2
       val[k1] <- (rp - val[k])
     }
     if (res.df <= 0) pval <- (.liu2(d,val) + .liu2(d1,val))/2 else
     pval <- (.simf(d,val,res.df) + .simf(d1,val,res.df))/2
  } else { pval <- 2 }
  if (pval > .5) {
    if (res.df <= 0) pval <- (pchisq(d,df=rank1,lower.tail=FALSE)+pchisq(d1,df=rank1,lower.tail=FALSE))/2 else
    pval <- (pf(d/rank1,rank1,res.df,lower.tail=FALSE)+pf(d1/rank1,rank1,res.df,lower.tail=FALSE))/2
  }
  list(stat=d,pval=min(1,pval),rank=rank)
}


## For t2() terms egpd keeps the coefficients, their covariance and the design in the
## basis the smooths were fitted in (see .predictable.gam). That basis mixes the
## intercept into the smooth, so summary() maps them to mgcv's prediction basis, as
## summary.gam reports them: X = X_fit P^-1, beta = P beta_fit, V = P V_fit P'.
.summary_prediction_basis <- function(obj) {
P <- obj[["P"]]
if (is.null(P)) return(obj)
X <- obj$X
if (rcond(P) > sqrt(.Machine$double.eps)) {
  obj$X <- X %*% solve(P)
} else {
  ## A fully penalised random effect can make the prediction-to-fit map
  ## rank deficient even with full penalised information. Use the actual
  ## training prediction design, not an arbitrary generalized inverse.
  xp <- obj[["summary.X"]]
  if (is.null(xp)) {
    if (!is.data.frame(obj$model) || nrow(obj$model) != nrow(X))
      stop("Training prediction design unavailable for this smooth summary; refit the model.", call. = FALSE)
    g <- obj
    class(g) <- "gam"
    xp <- .gam_predict_rows(g, obj$model, "lpmatrix")
  }
  obj$X <- xp
}
dimnames(obj$X) <- dimnames(X)
.gam_prediction_basis(obj)
}

.boundary_smooth <- function(sm) length(sm$S) > 0L && isTRUE(sm$null.space.dim == 0L)

## Known-scale version of Wood's random-effect test (mgcv's reTest/recov).
## Work in the joint fit/penalty basis: separate predictor designs do not
## contain the likelihood cross-information between distributional parameters.
.re_test_setup <- function(object) {
  ld <- object$likdata
  if (!isTRUE(object$convergence == 0L)) stop("the fit did not converge")
  if (isTRUE(ld$sandwich)) stop("no cluster-adjusted null score law is available")
  H <- object[["re.info"]]
  if (is.null(H)) stop("joint likelihood information was not retained; refit the model")
  if (!isTRUE(object$hessian.rank[1L] == object$hessian.rank[2L]))
    stop("the penalised information is numerically singular")
  free <- which(!ld$idpars %in% which(vapply(ld$X, .is_pinned_block, logical(1))))
  d <- .coef_scale(ld)[free]
  dd <- outer(d, d)
  H <- H[free, free, drop = FALSE] / dd
  e <- eigen(H, symmetric = TRUE)
  tol <- nrow(H) * .Machine$double.eps * max(abs(e$values))
  if (any(e$values < -tol)) stop("the joint likelihood information is indefinite")
  pos <- e$values > tol
  if (!any(pos)) stop("the joint likelihood information has zero rank")
  random <- unlist(lapply(object$gotsmooth, function(i) {
    block <- which(ld$idpars == i)
    unlist(lapply(object[[i]]$smooth, function(sm)
      if (isTRUE(sm$random) && .boundary_smooth(sm)) block[sm$first.para:sm$last.para]))
  }))
  list(R = sqrt(e$values[pos]) * t(e$vectors[, pos, drop = FALSE]),
       S = ld$gamma * ld$S[free, free, drop = FALSE] / dd,
       V = object$Vp[free, free, drop = FALSE] * dd, b = coef(object)[free] * d,
       random = match(random, free), free = free)
}

.re_test <- function(re, id) {
  if (inherits(re, "error")) stop(re)
  id <- match(id, re$free)
  R <- re$R
  other <- setdiff(re$random, id)
  if (length(other)) {
    ## With H = R'R, integrating the other random effects adds
    ## H[,o] S[o,o]^-1 H[o,] to the null score covariance.
    L <- chol(re$S[other, other, drop = FALSE])
    R <- rbind(R, forwardsolve(t(L), crossprod(R[, other, drop = FALSE], R)))
  }
  ## Profile the remaining nuisance coefficients, including ordinary smooths
  ## in any predictor. Check what the data identify BEFORE adding the tested
  ## term's penalty, which can be arbitrarily large even for an aliased term.
  n <- setdiff(seq_along(re$b), c(other, id))
  Z <- R[, id, drop = FALSE]
  K <- ZZ <- crossprod(Z)
  if (length(n)) {
    X <- R[, n, drop = FALSE]
    L <- chol(crossprod(X) + re$S[n, n, drop = FALSE])
    K <- ZZ - crossprod(forwardsolve(t(L), crossprod(X, Z)))
  }
  if (max(K) <= ncol(R) * .Machine$double.eps^.8 * max(ZZ))
    stop("no estimable random-effect contrast remains after accounting for nuisance terms")
  Rm <- chol(K + re$S[id, id, drop = FALSE])
  stat <- sum((Rm %*% re$b[id])^2)
  V <- re$V[, id, drop = FALSE]
  B <- R %*% V
  ## Null law of ||Rm b_hat||^2: independent chi-squares with these weights,
  ## not a chi-square on the edf, nor a universal 50:50 boundary mixture.
  w <- svd(B %*% t(Rm), nu = 0L, nv = 0L)$d^2
  w <- w[w > max(w) * .Machine$double.eps^.8]
  if (!length(w)) stop("the null covariance has zero rank")
  p <- if (length(w) == 1L) pchisq(stat / w, 1, lower.tail = FALSE) else
    mgcv::psum.chisq(stat / max(w), w / max(w))
  if (!is.finite(p) || p < 0 || p > 1) stop("the weighted chi-square tail could not be evaluated")
  c(stat, p)
}

.giveSmoothTest <- function(object, i, j, re) {
  obj <- object[[i]]
  sm <- obj$smooth[[j]]
  id <- sm$first.para:sm$last.para
  edf <- min(length(id), sum(obj$edf[id]))
  if (.boundary_smooth(sm)) {
    z <- tryCatch(.re_test(re, which(object$idpars == i)[id]), error = function(e) {
      warning(sprintf("Boundary test unavailable for %s in %s: %s.",
                      sm$label, names(object)[i], conditionMessage(e)), call. = FALSE)
      c(NA_real_, NA_real_)
    })
    out <- data.frame(edf, length(id), z[1L], z[2L])
  } else {
    obj <- .summary_prediction_basis(obj)
    stat <- .testStat(coef(obj)[id], obj$X[, id, drop = FALSE], obj$Vp[id, id, drop = FALSE], edf)
    out <- data.frame(stat$rank, length(id), stat$stat, stat$pval)
  }
  names(out) <- c("edf", "max.df", "Chi.sq", "p-value")
  rownames(out) <- sm$label
  out
}

.giveParametricTest <- function(obj) {
obj <- .summary_prediction_basis(obj)
np <- length(coef(obj))
id <- !logical(np)
for (j in seq_along(obj$smooth)) {
if (!is.null(obj$smooth[[j]]$first.para)) {
idsmooth <- obj$smooth[[j]]$first.para:obj$smooth[[j]]$last.para
id[idsmooth] <- FALSE
}
}
id2 <- which(id)
beta <- coef(obj)[id]
ese <- sqrt(obj$Vp[cbind(id2, id2)])
## Wald z statistics (known-scale likelihood), two-sided as in mgcv's summary.gam
tval <- beta / ese
pval <- 2 * pnorm(-abs(tval))
## A parameter held by `fixed` has one inert, all-zero "(fixed)" design column whose
## coefficient is not estimated: report NA rather than 0 with 0/0 statistics.
if (.is_pinned_block(obj$X))
  beta[] <- ese[] <- tval[] <- pval[] <- NA_real_
out <- try(data.frame(beta, ese, tval, pval))
rownames(out) <- colnames(obj$X[,id, drop=FALSE])
out
}

.tidySmoothTable <- function(tab) {
for (i in 1:3) tab[,i] <- round(tab[,i], 2)
tab[,4] <- signif(tab[,4], 3)
tab[,4] <- replace(unlist(tab[,4]), unlist(tab[,4]) < 2e-16, "<2e-16")
tab
}

## The one-line note of print() and summary() on a fit whose Hessian behind the standard
## errors is singular: `h` is fit$hessian.rank, c(rank, npar) (.hessian_check), and
## `component` how the note refers to it; NULL when it has full rank, or was not assessed
## (NA, or a fit from before it was recorded).
.hessian_note <- function(h, component = "`hessian.rank`") {
if (length(h) != 2L || anyNA(h) || h[[1L]] >= h[[2L]])
  return(NULL)
sprintf(paste0("Note: the Hessian behind the standard errors is singular (numerical rank ",
               "%d of %d estimated coefficients; %s, see ?egpd): they do not ",
               "reflect the uncertainty along the unidentified direction(s)."),
        as.integer(h[[1L]]), as.integer(h[[2L]]), component)
}

.tidyParametricTable <- function(tab) {
for (i in 1:3) tab[,i] <- round(tab[,i], 2)
tab[,4] <- signif(tab[,4], 3)
tab[,4] <- replace(unlist(tab[,4]), unlist(tab[,4]) < 2e-16, "<2e-16")
tab
}

.smooth.summary.egpd <- function(object, re = tryCatch(.re_test_setup(object), error = identity)) {
## `re` is evaluated once, only if a boundary test is requested.
out <- lapply(object$gotsmooth, function(i) do.call(rbind,
  lapply(seq_along(object[[i]]$smooth), function(j) .giveSmoothTest(object, i, j, re))))
nms <- names(object)[object$gotsmooth]
nms[nms == "mu"] <- "location"
nms[nms == "lpsi"] <- "logscale"
nms[nms == "xi"] <- "shape"
names(out) <- nms
out
}

.parametric.summary.egpd <- function(object) {
idgamlist <- sapply(object, class) == "gamlist"
npar <- sum(idgamlist)
out <- lapply(which(idgamlist), function(i) .giveParametricTest(object[[i]]))
for (i in seq_along(out)) names(out[[i]]) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
nms <- names(object)[idgamlist]
nms[nms == "mu"] <- "location"
nms[nms == "lpsi"] <- "logscale"
nms[nms == "xi"] <- "shape"
names(out) <- nms
out
}

## Newton functions

.zerosinvec <- function(x, id) {
out <- numeric(length(id))
out[id] <- x
out
}

.zerosinmat <- function(x, id) {
out <- matrix(0, length(id), length(id))
out[id, id] <- x
out
}

.newton_step_inner <- function(pars, fn, sfn, ..., control, trace=0, kept=NULL, alpha0=1) {

steptol <- control$steptol
itlim <- control$itlim
fntol <- control$fntol
gradtol <- control$gradtol
stepmax <- control$stepmax

pars0 <- pars
if (is.null(kept))
  kept <- !logical(length(pars))

it <- 1
okay <- TRUE
f0 <- fn(pars, ...)
if (isTRUE(attr(f0, "inner.convergence") != 0L))
  stop("Inner coefficient optimization did not converge: ",
       attr(f0, "inner.report"), call.=FALSE)
step1 <- NULL

while (okay) {
  if (it > 1)
    g0 <- g
  if (!is.null(step1)) {
    step0 <- step1
    g <- attr(step0, "gradient")
  } else {
    attr(pars, "beta") <- attr(f0, "beta")
    step0 <- sfn(pars, ..., kept=kept)
    g <- attr(step0, "gradient")
  }
  if (trace)
    .itreport(f0, g, it - 1)
  if (mean(abs(g)) < gradtol) {
    report <- c("gradient tolerance reached")
    break
  }

  ## Limit the step componentwise, unless that turns it into an ascent direction. A
  ## nearly singular Hessian (e.g. the outer-product Hessian of the zero-truncated
  ## families along a flat direction) asks for huge steps in some coordinates, and
  ## clipping only those can reverse the sign of g'step; then no step length lowers
  ## the objective and the search stopped on "step tolerance reached" far from the
  ## optimum. Rescaling the whole step keeps the descent direction.
  decrement <- sum(g * step0)
  clipped <- sign(step0) * pmin(abs(step0), stepmax)
  if (sum(g * clipped) <= 0 && decrement > 0)
    clipped <- step0 * (stepmax / max(abs(step0)))
  ## When the search direction is already the gradient there is nothing to fall back on.
  sd.tried <- isTRUE(max(abs(as.vector(step0) - as.vector(g))) == 0)
  step0 <- clipped
  alpha <- alpha0
  report <- NULL
  ls <- TRUE
  while(ls & is.null(report)) {
    step <- alpha * step0
    stepokay <- mean(abs(step)) > steptol
    if (!stepokay) {
      ## The Newton direction predicted a real decrease (half the Newton decrement
      ## exceeds fntol) that no step along it achieves, so the quadratic model is
      ## wrong here. Try steepest descent once before stopping, rather than report a
      ## point with a large gradient as converged.
      if (!sd.tried && .5 * decrement > fntol) {
        sd.tried <- TRUE
        step0[] <- sign(g) * pmin(abs(g), stepmax)
        alpha <- 1
        next
      }
      report <- c("step tolerance reached")
    } else {
      theta1 <- pars - step
      f1 <- fn(theta1, ...)
      d <- f1 - f0
      if (isTRUE(attr(f1, "inner.convergence") != 0L)) d <- Inf
      if (!is.finite(d))
        d <- 10
      if (d < 0) {
        attr(theta1, "beta") <- attr(f1, "beta")
        ## A trial point whose search direction fails is rejected, except when the
        ## analytic third derivatives are unavailable there: that condition must
        ## reach .outer(), which retries with finite differences (see .BFGS).
        step1 <- tryCatch(sfn(theta1, ..., kept=kept), error = function(e) {
          if (inherits(e, "egpd_d34_unavailable")) stop(e)
          structure(conditionMessage(e), class="try-error", condition=e)
        })
      if (inherits(step1, "try-error"))
        d <- 1
      if (any(!is.finite(attr(step1, "gradient"))) | any(!is.finite(attr(step1, "PP"))))
        d <- 1
      }
      if (d < 0) {
        f0 <- f1
        pars <- theta1
        step0 <- step1
        ls <- FALSE
        ## A negligible decrease ends the search only if the step was not predicted
        ## to do far better: a tiny gain from a step for which the local quadratic
        ## model promised a large one (a poor Hessian, cut back by the line search)
        ## is a stall, not convergence, so keep iterating.
        if (-d < fntol && !isTRUE(.5 * decrement > 1e3 * fntol))
          report <- "function tolerance reached"
      } else {
        alpha <- .5 * alpha
      }
    }
  }
  if (!is.null(report))
    break
  if (mean(abs(attr(step0, "gradient"))) < gradtol) {
    report <- "gradient tolerance reached"
    break
  }
  if (it >= itlim) {
    report <- c("iteration limit reached")
    okay <- FALSE
  } else {
    it <- it + 1
  }
}

convergence <- as.integer(identical(report, "iteration limit reached"))
## The objectives map invalid parameter values to 1e20 and zero their gradients, so
## a start in an invalid region satisfies the gradient test at once. Never report
## that as convergence.
if (.egpd_invalid_objective(f0)) {
  convergence <- 2L
  report <- "objective not finite at the returned parameters"
}
convergence <- .inner_stall(f0, convergence, report)
report <- attr(convergence, "report")
convergence <- as.vector(convergence)
if (trace)
  cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(pars=as.vector(pars), objective=f0)
out$gradient <- attr(step0, "gradient")
## g'H^{-1}g at the returned point (the Newton decrement; twice the decrease the
## local quadratic model still predicts). When the search stopped inside a line
## search, step0 has been clipped, so use the decrement computed before clipping.
out$decrement <- if (identical(report, "step tolerance reached")) decrement else
  sum(attr(step0, "gradient") * step0)
out$Hessian <- attr(step0, "Hessian")
if (!is.null(attr(step0, "PP")))
  kept <- .new.kept(attr(step0, "PP"), kept)
out$kept <- kept
out$cholHessian <- attr(step0, "cholH")
out$diagHessian <- attr(step0, "diagH")
out$rankHessian <- attr(step0, "rank")
out$convergence <- convergence
out$report <- report
out$iterations <- it
out$gradconv <- convergence == 0L && substr(report, 1, 4) == "grad"
if (!is.null(attr(pars, "beta")))
  out$beta <- attr(pars, "beta")
out
}

.rank2drop <- function(x) {
nc <- ncol(x)
R <- qr(x)
r <- R$rank
if (r < nc) {
  drop <- R$pivot[-seq_len(r)]
} else {
  drop <- integer(0)
}
drop
}

.new.kept <- function(x, kept) {
x <- x[kept, kept, drop=FALSE]
nc <- ncol(x)
R <- qr(x)
r <- R$rank
if (r < nc) {
  drop <- R$pivot[-seq_len(r)]
} else {
  drop <- integer(0)
}
if (length(drop) > 0)
  kept[which(kept)[drop]] <- FALSE
kept
}

.newton_step <- function(pars, fn, sfn, ..., control, trace=0, alpha0=1) {

nkept <- length(pars)
fit0 <- .newton_step_inner(pars, fn, sfn, ..., control=control, trace=trace, alpha0=alpha0)
nkept1 <- sum(fit0$kept)

while(nkept > nkept1) {
  nkept <- nkept1
  fit0 <- .newton_step_inner(fit0$par, fn, sfn, ..., control=control, trace=trace, kept=fit0$kept, alpha0=alpha0)
  nkept1 <- sum(fit0$kept)
}

fit0

}

## Final coefficient fit (no smoothing parameters to estimate, or fixed ones), with
## the retry the profile REML objective uses: when Newton stops short of the gradient
## criterion, take a few short gradient steps and resume.
##
## The other stopping rules do not establish convergence. A line search can stop on its
## step or function tolerance far from the optimum when the Hessian is a poor model
## there, which shows as a Newton step that still predicts a decrease of more than
## `dectol` (code 3); and with the outer-product (BHHH) Hessian of several families
## Newton converges only linearly and can run out of iterations next to the optimum
## (code 1). That Hessian can also overstate the remaining decrease by orders of
## magnitude along a poorly determined direction. So before such a fit is reported as
## not converged, a short BFGS search from the returned point, using only the objective
## and its gradient (.polish_inner), checks it:
## - if the search meets the gradient criterion, the fit has converged there;
## - if it gains `dectol` or more, Newton resumes from the better point, once;
## - otherwise a fit that ran out of iterations stays code 1 (at the better point),
##   and a stalled fit is accepted -- the predicted decrease was not realisable --
##   unless the searches keep gaining `dectol` or more (a ridge along which the
##   objective keeps decreasing), or the objective fails to fall along its own gradient
##   although the gradient and the curvature along it promise more than `dectol`
##   (.sd_pred; a likelihood that is numerically unreliable there, or a gradient that
##   does not match it). Those stay code 3.
## Fits that converge by the gradient criterion, or stop with a negligible predicted
## decrease, never run the check -- unless they stopped with a non-negligible gradient
## where the objective is not convex along it (.inner_saddle): there .perturb() makes
## the Newton step, and so the predicted decrease, negligible however far the point is
## from a minimum (e.g. on the unbounded likelihood of an EGPD shape below -1).
.fit_inner <- function(beta, fn, sfn, likdata, likfns, dectol=1e-3) {
fit <- .fit_inner_newton(beta, fn, sfn, likdata, likfns, dectol)
if (!is.null(likfns$finish_inner))
  fit <- likfns$finish_inner(fit, fn, sfn, likdata, likfns)
fit
}

.fit_inner_newton <- function(beta, fn, sfn, likdata, likfns, dectol) {
ctl <- likdata$control$inner
newton <- function(b) {
  fit <- .newton_step(b, fn, sfn, likdata=likdata, likfns=likfns, control=ctl)
  if (!fit$gradconv && fit$convergence != 2L) {
    ctl10 <- ctl
    ctl10$itlim <- 10
    fit <- .newton_step(fit$par, fn, sfn, likdata=likdata, likfns=likfns, control=ctl10,
                        newton=FALSE, alpha0=.05)
    fit <- .newton_step(fit$par, fn, sfn, likdata=likdata, likfns=likfns, control=ctl)
  }
  fit
}
saddle <- function(fit)
  .inner_saddle(fit, sfn, likdata, likfns, gradtol=ctl$gradtol)
unresolved <- function(fit)
  !fit$gradconv && (fit$convergence == 1L ||
                    (fit$convergence == 0L && isTRUE(.5 * fit$decrement > dectol)) ||
                    !is.null(saddle(fit)))
polish <- function(fit)
  .polish_inner(fit, fn, sfn, likdata, likfns, gradtol=ctl$gradtol,
                maxit=min(100, ctl$itlim))
fit <- newton(beta)
if (!unresolved(fit))
  return(fit)
pol <- polish(fit)
if (!pol$gradconv && pol$gain >= dectol) {
  fit <- newton(pol$par)
  if (!unresolved(fit)) {
    fit$report <- paste(fit$report, "(Newton resumed after a gradient search)")
    return(fit)
  }
  pol <- polish(fit)
}
stopped <- sprintf("%s, with max |gradient| %.3g and a further decrease of %.3g predicted by the Newton step",
                   fit$report, max(abs(fit$gradient)), .5 * fit$decrement)
if (pol$gradconv) {
  fit <- .refresh_inner(fit, pol$par, fn, sfn, likdata, likfns)
  fit$convergence <- 0L
  fit$gradconv <- TRUE
  fit$report <- sprintf("gradient tolerance reached by a gradient search after Newton stopped (%s)",
                        stopped)
  return(fit)
}
better <- if (pol$gain > 0) .refresh_inner(fit, pol$par, fn, sfn, likdata, likfns) else fit
if (fit$convergence == 1L)
  return(better)
sd <- .sd_pred(fit, sfn, likdata, likfns)
## a stop found unresolved by .inner_saddle alone (its Newton decrement is negligible)
why <- if (!isTRUE(.5 * fit$decrement > dectol)) saddle(fit)
if (pol$gain >= dectol) {
  fit <- better
  fit$convergence <- 3L
  fit$report <- sprintf("%s; gradient searches keep decreasing the objective (by %.3g in the last)%s",
                        stopped, pol$gain, if (is.null(why)) "" else paste0("; ", why))
} else if (!is.null(why)) {
  fit$convergence <- 3L
  fit$report <- sprintf("%s; %s", stopped, why)
} else if (sd > dectol && pol$gain < .01 * sd) {
  fit$convergence <- 3L
  fit$report <- sprintf(paste0("%s; a gradient search gained %.3g where the gradient and the ",
                               "curvature along it predict %s, so the objective and its gradient ",
                               "disagree (the likelihood may be numerically unreliable here)"),
                        stopped, pol$gain,
                        if (is.finite(sd)) sprintf("%.3g", sd) else "an unbounded decrease")
} else {
  fit <- better
  fit$report <- sprintf("%s, which was not realisable (a gradient search gained %.3g)",
                        stopped, pol$gain)
}
fit
}

## A coefficient search that stopped (code 0, short of the gradient criterion) with a
## non-negligible gradient where the objective is not convex along that gradient
## (.sd_pred finds no positive curvature there, or cannot evaluate it). The Hessian is
## then indefinite, and .perturb() adds a ridge that makes the Newton step, and the
## decrease it predicts, negligible, so the stop looks like convergence. It is no
## minimum: typically the search ran into the edge of the support of a likelihood that
## is unbounded there -- the continuous EGPD with a shape below -1, whose density grows
## without bound at the upper end of its support, so that every further step leaves the
## support. Returns a description, or NULL. A count family on its clamp, where the
## likelihood is flat, has a zero gradient there and is not affected.
.inner_saddle <- function(fit, sfn, likdata, likfns, gradtol) {
if (isTRUE(fit$gradconv) || !isTRUE(fit$convergence == 0L))
  return(NULL)
if (!isTRUE(mean(abs(fit$gradient)) >= gradtol) || .sd_pred(fit, sfn, likdata, likfns) < Inf)
  return(NULL)
paste0("the objective is not convex along its gradient there: a saddle point, or a ",
       "likelihood that is unbounded (as for an EGPD shape below -1 with an observation ",
       "at the upper end of the support)")
}

## A short BFGS search from a fit's coefficients, over the coefficients it kept, using
## only the objective and its analytic gradient. Returns the lowest point evaluated, its
## gain over the fit's objective, and whether the gradient criterion holds there; a
## failure part-way keeps what was found so far.
.polish_inner <- function(fit, fn, sfn, likdata, likfns, gradtol, maxit=100) {
par0 <- as.vector(fit$pars)
kept <- fit$kept
if (is.null(kept))
  kept <- !logical(length(par0))
idx <- which(kept)
f0 <- as.numeric(fit$objective)
best <- list(par=par0, value=f0)
full <- function(b) {
  p <- par0
  p[idx] <- b
  p
}
f <- function(b) {
  p <- full(b)
  v <- as.numeric(fn(p, likdata=likdata, likfns=likfns))
  if (!is.finite(v) || .egpd_invalid_objective(v))
    return(Inf)
  if (v < best$value)
    best <<- list(par=p, value=v)
  v
}
g <- function(b) {
  s <- sfn(full(b), likfns=likfns, likdata=likdata, kept=kept, newton=FALSE)
  out <- as.vector(attr(s, "gradient"))[idx]
  if (any(!is.finite(out)))
    stop("non-finite gradient")
  out
}
if (length(idx) && is.finite(f0))
  try(stats::optim(par0[idx], f, g, method="BFGS",
                   control=list(maxit=maxit, reltol=1e-12)), silent=TRUE)
gradconv <- FALSE
if (best$value < f0) {
  gb <- try(g(best$par[idx]), silent=TRUE)
  gradconv <- !inherits(gb, "try-error") && sum(abs(gb)) / length(par0) < gradtol
}
list(par=best$par, value=best$value, gain=f0 - best$value, gradconv=gradconv)
}

## The decrease a steepest-descent step should achieve from a fit's point, |g|^2 / (2 c),
## with c the curvature of the objective along its gradient g, from a difference of the
## analytic gradient (not the Hessian, which may be an approximation). Inf when that
## curvature is not positive or cannot be computed.
.sd_pred <- function(fit, sfn, likdata, likfns, h=1e-6) {
kept <- fit$kept
if (is.null(kept))
  kept <- !logical(length(fit$pars))
g0 <- as.vector(fit$gradient)
g0[!kept] <- 0
ng <- sqrt(sum(g0^2))
if (!is.finite(ng))
  return(Inf)
if (ng == 0)
  return(0)
u <- g0 / ng
s <- try(sfn(as.vector(fit$pars) + h * u, likfns=likfns, likdata=likdata, kept=kept,
             newton=FALSE), silent=TRUE)
if (inherits(s, "try-error"))
  return(Inf)
g1 <- as.vector(attr(s, "gradient"))
g1[!kept] <- 0
curv <- sum(u * (g1 - g0)) / h
if (!is.finite(curv) || curv <= 0)
  return(Inf)
ng^2 / (2 * curv)
}

## The Newton quantities of a fit at another point, without moving from it.
.refresh_inner <- function(fit, par, fn, sfn, likdata, likfns) {
s <- try(sfn(par, likfns=likfns, likdata=likdata, kept=fit$kept), silent=TRUE)
if (inherits(s, "try-error"))
  return(fit)
fit$pars <- as.vector(par)
fit$objective <- fn(par, likdata=likdata, likfns=likfns)
fit$gradient <- attr(s, "gradient")
fit$decrement <- sum(attr(s, "gradient") * s)
fit$Hessian <- attr(s, "Hessian")
fit$cholHessian <- attr(s, "cholH")
fit$diagHessian <- attr(s, "diagH")
fit$rankHessian <- attr(s, "rank")
fit
}

.itreport <- function(f, g, it) {
    report <- paste("\n Outer iteration ", it, ":", sep="")
    rep1 <- paste("  Outer max(|grad|):", signif(max(abs(g)), 3))
    rep2 <- paste("  Inner max(|grad|): ", signif(max(abs(attr(f, "gradient"))), 3), ".", sep="")
    report <- c(report, paste(rep1, rep2, sep="; "))
    cat(paste(report, collapse="\n"))
}

.BFGS <- function(pars, fn, gfn, ..., control, trace=0) {

steptol <- control$steptol
itlim <- control$itlim
fntol <- control$fntol
gradtol <- control$gradtol
stepmax <- control$stepmax

it <- 1
okay <- TRUE
f0 <- fn(pars, ...)
if (isTRUE(attr(f0, "inner.convergence") != 0L))
  stop("Inner coefficient optimization did not converge: ",
       attr(f0, "inner.report"), call.=FALSE)
g1 <- NULL
I <- iH <- H <- diag(length(pars))

while (okay) {
if (it > 1) g0 <- g
if (!is.null(g1)) {
g <- g1
} else {
attr(pars, "beta") <- attr(f0, "beta")
g <- gfn(pars, ...)
}
if (trace) .itreport(f0, g, it - 1)
if (mean(abs(g)) < gradtol) {
report <- c("gradient tolerance reached")
break
}
step0 <- crossprod(iH, g)
step0 <- sign(step0) * pmin(abs(step0), stepmax)
alpha <- 1
report <- NULL
ls <- TRUE
while(ls & is.null(report)) {
step <- alpha * step0
stepokay <- any(abs(step) > steptol)
if (!stepokay) {
report <- c("step tolerance reached")
} else {
theta1 <- pars - step
f1 <- fn(theta1, ...)
d <- f1 - f0
if (isTRUE(attr(f1, "inner.convergence") != 0L)) d <- Inf
if (!is.finite(d)) d <- Inf
if (d < 0) {
attr(theta1, "beta") <- attr(f1, "beta")
## The profile gradients refit the inner coefficients at shifted smoothing
## parameters and stop() when such a refit fails. A trial point whose gradient
## cannot be evaluated is rejected like one that does not decrease the objective,
## rather than aborting the whole fit -- except when the analytic third derivatives
## are unavailable there (egpd_d34_unavailable): .outer() must see that condition to
## retry with finite differences, or the search halts at the edge of the region
## where they are unavailable and reports convergence.
g1 <- tryCatch(gfn(theta1, ...), error = function(e)
  if (inherits(e, "egpd_d34_unavailable")) stop(e) else NULL)
if (is.null(g1) || any(!is.finite(g1))) {
  d <- 1
  g1 <- NULL
}
}
if (d < 0) {
yk <- g1 - g
denom <- sum(- yk * step)
if (is.finite(denom) && denom > 0) {
  t1 <- I - tcrossprod(- step, yk) / denom
  t2 <- I - tcrossprod(yk, - step) / denom
  t3 <- tcrossprod(- step) / denom
  iH <- t1 %*% iH %*% t2 + t3
  if (any(!is.finite(iH))) iH <- I
} else {
  iH <- I
}
f0 <- f1
pars <- theta1
g <- g1
ls <- FALSE
if (-d < fntol) report <- "function tolerance reached"
} else {
alpha <- .5 * alpha
}
}
}
if (!is.null(report)) break
if (mean(abs(g)) < gradtol) {
report <- "gradient tolerance reached"
break
}
if (it >= itlim) {
report <- c("iteration limit reached")
okay <- FALSE
} else {
it <- it + 1
}
}
convergence <- as.integer(identical(report, "iteration limit reached"))
if (.egpd_invalid_objective(f0)) {
  convergence <- 2L
  report <- "objective not finite at the returned parameters"
}
convergence <- .inner_stall(f0, convergence, report)
report <- attr(convergence, "report")
convergence <- as.vector(convergence)
if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(par=as.vector(pars), objective=f0)
out$gradient <- g
out$convergence <- convergence
out$report <- report
out$iterations <- it
if (!is.null(attr(pars, "beta"))) out$beta <- attr(pars, "beta")
out
}

## Objectives return 1e20 (or a non-finite value) at invalid parameters.
.egpd_invalid_objective <- function(f) {
f <- as.numeric(f)
length(f) != 1L || !is.finite(f) || f >= 1e19
}

## The profile REML objective (.reml0) marks, in attribute "inner.stall", a coefficient
## fit that stopped on a saddle point or an unbounded likelihood (.inner_saddle). The
## outer search is free to pass through such points, but if it ends on one its result
## is no fit, whatever the outer stopping rule says: code 3. Returns the code with the
## report as attribute "report".
.inner_stall <- function(f0, convergence, report) {
stall <- attr(f0, "inner.stall")
if (identical(convergence, 0L) && !is.null(stall)) {
  convergence <- 3L
  report <- paste0(paste(report, collapse=" "), "; but at the returned smoothing parameters ",
                   stall)
}
structure(convergence, report=report)
}
