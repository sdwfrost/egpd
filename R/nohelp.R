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


.giveSmoothTest <- function(obj, j) {
id <- obj$smooth[[j]]$first.para:obj$smooth[[j]]$last.para
p <- coef(obj)[id]
X <- obj$X[,id]
V <- obj$Vp[id, id]
edf1 <- sum(obj$edf[id])
edf1 <- min(length(id), edf1)
stat <- .testStat(p, X, V, edf1)
out <- data.frame(stat[[3]], length(id), stat[[1]], stat[[2]])
rownames(out) <- obj$smooth[[j]]$label
out
}

.giveParametricTest <- function(obj) {
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
tval <- beta / ese
pval <- pnorm(-abs(tval))
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

.tidyParametricTable <- function(tab) {
for (i in 1:3) tab[,i] <- round(tab[,i], 2)
tab[,4] <- signif(tab[,4], 3)
tab[,4] <- replace(unlist(tab[,4]), unlist(tab[,4]) < 2e-16, "<2e-16")
tab
}

.smooth.summary.egpd <- function(object) {
idgamlist <- sapply(object, class) == "gamlist"
npar <- sum(idgamlist)
out <- lapply(object$gotsmooth, function(i) do.call(rbind, lapply(seq_along(object[[i]]$smooth), function(j) .giveSmoothTest(object[[i]], j))))
for (i in seq_along(out)) names(out[[i]]) <- c("edf", "max.df", "Chi.sq", "Pr(>|t|)")
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
for (i in seq_along(out)) names(out[[i]]) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
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

  step0 <- sign(step0) * pmin(abs(step0), stepmax)
  alpha <- alpha0
  report <- NULL
  ls <- TRUE
  while(ls & is.null(report)) {
    step <- alpha * step0
    stepokay <- mean(abs(step)) > steptol
    if (!stepokay) {
      report <- c("step tolerance reached")
    } else {
      theta1 <- pars - step
      f1 <- fn(theta1, ...)
      d <- f1 - f0
      if (!is.finite(d))
        d <- 10
      if (d < 0) {
        attr(theta1, "beta") <- attr(f1, "beta")
        step1 <- try(sfn(theta1, ..., kept=kept), silent=TRUE)
      if (inherits(step1, "try-error"))
        d <- 1
      if (any(!is.finite(attr(step1, "gradient"))) | any(!is.finite(attr(step1, "PP"))))
        d <- 1
      }
      if (d < 0) {
        f0 <- f1
        pars <- theta1
        ls <- FALSE
      } else {
        if (d < fntol)
          report <- c("function tolerance reached")
        alpha <- .5 * alpha
      }
    }
  }
  if (!is.null(report))
    break
  it <- it + 1
  if (it == itlim) {
    report <- c("iteration limit reached")
    okay <- FALSE
  }
}

if (trace)
  cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(pars=as.vector(pars), objective=f0)
out$gradient <- attr(step0, "gradient")
out$Hessian <- attr(step0, "Hessian")
if (!is.null(attr(step0, "PP")))
  kept <- .new.kept(attr(step0, "PP"), kept)
out$kept <- kept
out$cholHessian <- attr(step0, "cholH")
out$diagHessian <- attr(step0, "diagH")
out$rankHessian <- attr(step0, "rank")
out$convergence <- 0
out$report <- report
out$iterations <- it
out$gradconv <- substr(report, 1, 4) == "grad"
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
stepokay <- all(abs(step) > steptol)
if (!stepokay) {
report <- c("step tolerance reached")
} else {
theta1 <- pars - step
f1 <- fn(theta1, ...)
d <- f1 - f0
if (d < 0) {
attr(theta1, "beta") <- attr(f1, "beta")
g1 <- gfn(theta1, ...)
if (any(!is.finite(g1))) d <- 1
yk <- g1 - g
denom <- sum(- yk * step)
t1 <- I - tcrossprod(- step, yk) / denom
t2 <- I - tcrossprod(yk, - step) / denom
t3 <- tcrossprod(- step) / denom
iH <- t1 %*% iH %*% t2 + t3
if (any(!is.finite(iH))) d <- 1
}
if (d < 0) {
f0 <- f1
pars <- theta1
ls <- FALSE
} else {
if (d < fntol) {
report <- c("function tolerance reached")
}
alpha <- .5 * alpha
}
}
}
if (!is.null(report)) break
it <- it + 1
if (it == itlim) {
report <- c("iteration limit reached")
okay <- FALSE
}
}
if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(par=as.vector(pars), objective=f0)
out$gradient <- g
out$convergence <- 0
out$report <- report
out$iterations <- it
if (!is.null(attr(pars, "beta"))) out$beta <- attr(pars, "beta")
out
}
