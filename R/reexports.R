## Re-exports of the Poisson-inverse Gaussian family constructors from
## gamlss.dist, so they are available directly from egpd alongside the
## package's own gamlss families (e.g. DEGPD1). gamlss.dist already provides
## complete gamlss.family objects for PIG and ZIPIG.

#' @importFrom gamlss.dist PIG
#' @export
gamlss.dist::PIG

#' @importFrom gamlss.dist ZIPIG
#' @export
gamlss.dist::ZIPIG
