#!/usr/bin/env Rscript
# Recovery study for the identity-link DEGPD m=1 GAM. Simulates discrete-EGPD counts with a
# known shape (including xi<0, bounded tails) and shows the identity link recovers the true
# shape across its sign, while the default log-link DEGPD GAM floors at 0 for xi<=0.
# Self-contained: depends only on egpd + simulated data (no external files).
suppressMessages(library(egpd)); set.seed(1)

Fd  <- function(t,s,xi,k) { H <- if (abs(xi)<1e-8) 1-exp(-t/s) else 1-(1+xi*t/s)^(-1/xi); H^k }
pmf <- function(k,s,xi,kp) Fd(k+1,s,xi,kp)-Fd(k,s,xi,kp)
sim <- function(n,s,xi,kp){ ep <- if(xi<0) floor(s*(-1/xi)) else 80
  ks<-0:(ep+1); p<-pmf(ks,s,xi,kp); p[!is.finite(p)|p<0]<-0; p<-p/sum(p); sample(ks,n,TRUE,p) }
fit_xi <- function(y, link){
  df<-data.frame(cases=y)
  args<-if(link=="identity") list(m=1,link="identity") else list(m=1)
  f<-tryCatch(egpd(list(lsigma=cases~1,lxi=~1,lkappa=~1),data=df,family="degpd",degpd.args=args),error=function(e)NULL)
  if(is.null(f)) return(NA_real_)
  as.data.frame(predict(f,newdata=df[1,,drop=FALSE],type="response"))$shape[1]
}

cat("=== Recovery: true xi vs mean xi_hat (n=3000, 15 reps), sigma=7 kappa=1.3 ===\n")
cat(sprintf("%6s | %-22s | %-22s\n","true","identity-link xi_hat","log-link xi_hat"))
for (xi0 in c(-0.30,-0.20,-0.10,0.00,0.10,0.20)) {
  id<-lg<-numeric(0)
  for (r in 1:15) { y<-sim(3000,7,xi0,1.3); id<-c(id,fit_xi(y,"identity")); lg<-c(lg,fit_xi(y,"log")) }
  cat(sprintf("%+6.2f | mean %+.3f sd %.3f (n=%d) | mean %+.3f sd %.3f\n",
      xi0, mean(id,na.rm=T), sd(id,na.rm=T), sum(!is.na(id)), mean(lg,na.rm=T), sd(lg,na.rm=T)))
}
cat("Done.\n")
