#!/usr/bin/env Rscript
# Extensive verification of the identity-link bounded-support DEGPD m=1
# implementation (src/degpd1_identity.cpp). Four independent checks.
suppressMessages(library(egpd))
set.seed(42)
ok <- TRUE
chk <- function(name, pass, extra="") { pass <- isTRUE(pass); cat(sprintf("[%s] %s %s\n", ifelse(pass,"PASS","FAIL"), name, extra)); if(!pass) ok<<-FALSE }

## reference R implementation of the model (independent of the C++) ----
Fdisc <- function(t, sigma, xi, kappa) (1 - (1 + xi*t/sigma)^(-1/xi))^kappa   # F(t)=H(t)^kappa
pmf   <- function(k, sigma, xi, kappa) Fdisc(k+1, sigma, xi, kappa) - Fdisc(k, sigma, xi, kappa)

## raw C++ callers (intercept-only design) ----
mk <- function(n) list(X=matrix(1,n,1), dupid=integer(0), dcate=0L,
                       off=list(numeric(0),numeric(0),numeric(0)))
nll0 <- function(y, lsig, xi, lkap) {
  d <- mk(length(y))
  egpd:::degpd1id0(list(lsig, xi, lkap), d$X,d$X,d$X, as.numeric(y), d$dupid, d$dcate, d$off, Inf)
}
gh   <- function(y, lsig, xi, lkap) {
  d <- mk(length(y))
  m <- egpd:::degpd1id12(list(lsig, xi, lkap), d$X,d$X,d$X, as.numeric(y), d$dupid, d$dcate, d$off, Inf)
  colSums(m)
}

## ---- TEST 1: analytic grad/Hessian vs finite differences ----
## (FD reference is unreliable in the |xi|<0.05 catastrophic-cancellation band,
##  where the formula matches exact sympy but loses float precision; the xi~0
##  region is covered exactly by TEST 1b below.)
maxg <- 0; maxh <- 0
for (rep in 1:50) {
  xi <- sample(c(-1,1),1) * runif(1,0.05,0.45); lsig <- runif(1,-0.5,1.5); lkap <- runif(1,-0.5,1)
  sig <- exp(lsig); kap <- exp(lkap)
  endpoint <- if (xi<0) sig*(-1/xi) else Inf
  ymax <- if (is.finite(endpoint)) max(1, floor(0.6*endpoint)) else 12   # stay well interior
  y <- sample(0:ymax, 30, replace=TRUE)
  pv <- pmf(y,sig,xi,kap); y <- y[is.finite(pv) & pv > 0]
  if (length(y) < 5) next
  g <- gh(y, lsig, xi, lkap)               # [g_s,g_x,g_k, hss,hsx,hsk,hxx,hxk,hkk]
  h <- 1e-5
  # FD gradient
  fdg <- c(
    (nll0(y,lsig+h,xi,lkap)-nll0(y,lsig-h,xi,lkap))/(2*h),
    (nll0(y,lsig,xi+h,lkap)-nll0(y,lsig,xi-h,lkap))/(2*h),
    (nll0(y,lsig,xi,lkap+h)-nll0(y,lsig,xi,lkap-h))/(2*h))
  maxg <- max(maxg, abs(fdg - g[1:3]))
  # FD Hessian (diagonal + a cross term)
  f <- function(a,b,c) nll0(y,a,b,c)
  hss <- (f(lsig+h,xi,lkap)-2*f(lsig,xi,lkap)+f(lsig-h,xi,lkap))/h^2
  hxx <- (f(lsig,xi+h,lkap)-2*f(lsig,xi,lkap)+f(lsig,xi-h,lkap))/h^2
  hkk <- (f(lsig,xi,lkap+h)-2*f(lsig,xi,lkap)+f(lsig,xi,lkap-h))/h^2
  hsx <- (f(lsig+h,xi+h,lkap)-f(lsig+h,xi-h,lkap)-f(lsig-h,xi+h,lkap)+f(lsig-h,xi-h,lkap))/(4*h^2)
  maxh <- max(maxh, abs(c(hss,hxx,hkk,hsx) - g[c(4,7,9,5)]))
}
chk("grad vs FD (|xi|>0.05)",  maxg < 1e-4, sprintf("(max abs err %.2e)", maxg))
chk("Hessian vs FD sanity (noisy)", maxh < 1e3, sprintf("(max abs err %.2e; FD 2nd-diff is noisy, see 1b)", maxh))

## ---- TEST 1b: exact-value regression vs high-precision SymPy (grad AND Hessian) ----
## (ls=0.8, lk=0.2; 10-sig-fig reference from SymPy at rational inputs)
refg <- list("3_-0.3"  = c(-0.9638081388, 0.2231744119, -0.8341773583),
             "5_0.2"   = c(-0.9190876232,-0.05040498104,-0.8212996840),
             "8_-0.15" = c(-6.467901819, -12.47460562,  -0.9952867043))
refh <- list("3_-0.3"  = c(3.928396593,3.174867172,0.4886081773,1.839673372,0.4544155908,0.1641182419),
             "5_0.2"   = c(1.389734953,1.627632807,0.3148211907,1.265261142,0.3344168954,0.1784254203),
             "8_-0.15" = c(17.09650638,56.29183438,0.03986485342,153.5321665,0.09456468769,0.004711560557))
meg <- 0; meh <- 0
for (key in names(refg)) { p <- as.numeric(strsplit(key,"_")[[1]])
  d <- mk(1); row <- egpd:::degpd1id12(list(0.8,p[2],0.2), d$X,d$X,d$X, p[1], d$dupid,d$dcate,d$off, Inf)[1,]
  meg <- max(meg, abs(row[1:3] - refg[[key]])); meh <- max(meh, abs(row[4:9] - refh[[key]])) }
chk("grad == SymPy-exact (3 fixed pts)",    meg < 1e-6, sprintf("(max abs err %.2e)", meg))
chk("Hessian == SymPy-exact (3 fixed pts)", meh < 1e-6, sprintf("(max abs err %.2e)", meh))

## ---- TEST 2: nll(identity, xi>0) == nll(log-link degpd1d0) ----
d <- mk(8); y <- 0:7; lsig <- 0.3; xi <- 0.2; lkap <- 0.1
a <- nll0(y, lsig, xi, lkap)
b <- egpd:::degpd1d0(list(lsig, log(xi), lkap), d$X,d$X,d$X, as.numeric(y), d$dupid, d$dcate, d$off, Inf)
chk("identity nll == log-link nll (xi>0)", abs(a-b) < 1e-9, sprintf("(diff %.2e)", abs(a-b)))

## ---- TEST 3: pmf sums to 1 for xi<0 (bounded), via exp(-per-obs nll) ----
sig <- 5; xi <- -0.3; kap <- 1.4; ep <- floor(sig*(-1/xi))
ks <- 0:(ep+2)
ps <- sapply(ks, function(k) if (1 + xi*k/sig > 0) exp(-nll0(k, log(sig), xi, log(kap))) else 0)  # k in support
chk("pmf sums to 1 (xi<0)", abs(sum(ps)-1) < 1e-6, sprintf("(sum=%.6f, endpoint=%d)", sum(ps), ep))

## ---- TEST 3b: endpoint-cell grad/Hessian vs FD (the P=1-F(y) branch) ----
y0 <- ep                                  # the cell containing the non-integer endpoint
g <- gh(y0, log(sig), xi, log(kap)); h <- 1e-5
fdg <- c((nll0(y0,log(sig)+h,xi,log(kap))-nll0(y0,log(sig)-h,xi,log(kap)))/(2*h),
         (nll0(y0,log(sig),xi+h,log(kap))-nll0(y0,log(sig),xi-h,log(kap)))/(2*h),
         (nll0(y0,log(sig),xi,log(kap)+h)-nll0(y0,log(sig),xi,log(kap)-h))/(2*h))
chk("endpoint-cell grad vs FD", max(abs(fdg-g[1:3])) < 1e-4, sprintf("(err %.2e)", max(abs(fdg-g[1:3]))))

## ---- TEST 4: GAM recovers a KNOWN xi<0; log-link floors at 0 ----
sim_one <- function(n, sigma, xi, kappa) {
  ep <- if (xi<0) floor(sigma*(-1/xi)) else 60
  ks <- 0:(ep+1); p <- pmf(ks, sigma, xi, kappa); p[!is.finite(p)|p<0] <- 0; p <- p/sum(p)
  sample(ks, n, replace=TRUE, prob=p)
}
y <- sim_one(4000, sigma=6, xi=-0.30, kappa=1.5)
df <- data.frame(cases=y)
fit_id <- egpd(list(lsigma=cases~1, lxi=~1, lkappa=~1), data=df,
               family="degpd", degpd.args=list(m=1, link="identity"))
fit_lg <- egpd(list(lsigma=cases~1, lxi=~1, lkappa=~1), data=df,
               family="degpd", degpd.args=list(m=1))
xi_id <- as.data.frame(predict(fit_id, newdata=df[1,,drop=FALSE], type="response"))$shape[1]
xi_lg <- as.data.frame(predict(fit_lg, newdata=df[1,,drop=FALSE], type="response"))$shape[1]
chk("GAM identity recovers xi<0 (truth -0.30)", xi_id < -0.15 && xi_id > -0.45, sprintf("(xi_hat=%.3f)", xi_id))
chk("GAM log-link floors at >=0 on same data",  xi_lg >= -1e-3, sprintf("(xi_hat=%.4f)", xi_lg))

cat(sprintf("\n==== %s ====\n", ifelse(ok,"ALL CHECKS PASSED","SOME CHECKS FAILED")))
