## Analytical derivatives for gamlss families using the Deriv package
##
## Log-density/PMF functions are symbolically differentiated on first use
## (lazy initialisation).  The results are cached in .deriv_cache and
## dispatched via .egpd_ad3() / .egpd_ad4(), which mirror the numerical
## .egpd_nd3() / .egpd_nd4().

## ---------------------------------------------------------------
## Continuous EGPD log-density functions (plain R — no Deriv calls)
## ---------------------------------------------------------------

## Type 1 / model 1 (EGPD1): G(u) = u^kappa
.logdens_egpd1 <- function(y, mu, sigma, nu) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  log(nu) + (nu - 1) * log(u) - log(mu) - (1 / sigma + 1) * log(v)
}

## Type 4 / model 3 (EGPD3): g4 simplifies to avoid dbeta/pbeta
.logdens_egpd3 <- function(y, mu, sigma, nu) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  log(nu) + log(1 - (1 - u)^nu) - lbeta(1 / nu, 2) -
    log(mu) - (1 / sigma + 1) * log(v)
}

## Type 2 / model 5 (EGPD5): truncated normal
.logdens_egpd5 <- function(y, mu, sigma, nu) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  C <- pnorm(sqrt(nu), 0, 1) - 0.5
  log(nu) - 0.5 * log(2 * pi) - 0.5 * nu * (u - 1)^2 -
    log(C) - log(mu) - (1 / sigma + 1) * log(v)
}

## Type 3 / model 6 (EGPD6): truncated beta (main part without log(C))
## d.G type 3 evaluates dbeta(u, κ, κ)/C with u = GPD CDF directly
## (no rescaling to [1/32, 0.5] — the density is 0 outside that range).
.logdens_egpd6_main <- function(y, mu, sigma, nu) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  (nu - 1) * log(u) + (nu - 1) * log(1 - u) - lbeta(nu, nu) -
    log(mu) - (1 / sigma + 1) * log(v)
}

## Normalising constant for type 3
.logC_egpd6 <- function(nu) {
  log(pbeta(0.5, nu, nu) - pbeta(1 / 32, nu, nu))
}

## Numerical d/dnu of log(C) for type 3
.dlogC_dnu_egpd6 <- function(nu) {
  h <- pmax(abs(nu) * 1e-6, 1e-8)
  (.logC_egpd6(nu + h) - .logC_egpd6(nu - h)) / (2 * h)
}

## Type 5 / model 4 (EGPD4): power-beta (Deriv-compatible part)
.logdens_egpd4_deriv <- function(y, mu, sigma, nu, tau) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  log(tau / 2) + log(nu) + log(1 - (1 - u)^nu) - lbeta(1 / nu, 2) -
    log(mu) - (1 / sigma + 1) * log(v)
}

## G4 CDF term for type 5: (tau/2 - 1) * log(G4(u; nu))
.logG4_term <- function(y, mu, sigma, nu, tau) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  G4 <- 1 - pbeta((1 - u)^nu, 1 / nu, 2)
  (tau / 2 - 1) * log(G4)
}

## d/dtau of G4 term = (1/2)*log(G4) (hand-coded)
.ad_dldt_egpd4_G4_fn <- function(y, mu, sigma, nu, tau) {
  v <- 1 + sigma * y / mu
  u <- 1 - v^(-1 / sigma)
  G4 <- 1 - pbeta((1 - u)^nu, 1 / nu, 2)
  0.5 * log(G4)
}

## ---------------------------------------------------------------
## Discrete EGPD log-PMF functions (plain R — no Deriv calls)
## ---------------------------------------------------------------

## DEGPD1 (type 1): P = u1^nu - u0^nu
.logpmf_degpd1 <- function(y, mu, sigma, nu) {
  v1 <- 1 + sigma * (y + 1) / mu
  u1 <- 1 - v1^(-1 / sigma)
  v0 <- 1 + sigma * y / mu
  u0 <- 1 - v0^(-1 / sigma)
  log(u1^nu - u0^nu)
}

## DEGPD3 (type 4): G4(u) = 1 - pbeta((1-u)^d, 1/d, 2)
.logpmf_degpd3 <- function(y, mu, sigma, nu) {
  v1 <- 1 + sigma * (y + 1) / mu
  u1 <- 1 - v1^(-1 / sigma)
  v0 <- 1 + sigma * y / mu
  u0 <- 1 - v0^(-1 / sigma)
  G4_1 <- 1 - pbeta((1 - u1)^nu, 1 / nu, 2)
  G4_0 <- 1 - pbeta((1 - u0)^nu, 1 / nu, 2)
  log(G4_1 - G4_0)
}

## DEGPD5 (type 2): truncated normal CDF
.logpmf_degpd5 <- function(y, mu, sigma, nu) {
  v1 <- 1 + sigma * (y + 1) / mu
  u1 <- 1 - v1^(-1 / sigma)
  v0 <- 1 + sigma * y / mu
  u0 <- 1 - v0^(-1 / sigma)
  sd_val <- 1 / sqrt(nu)
  Fmin <- 1 - pnorm(1, 0, sd_val)
  Fmax <- pnorm(1, 1, sd_val)
  C <- Fmax - Fmin
  G2_1 <- (pnorm(u1, 1, sd_val) - Fmin) / C
  G2_0 <- (pnorm(u0, 1, sd_val) - Fmin) / C
  log(G2_1 - G2_0)
}

## DEGPD6 (type 3): truncated beta CDF
.logpmf_degpd6_main <- function(y, mu, sigma, nu) {
  v1 <- 1 + sigma * (y + 1) / mu
  u1 <- 1 - v1^(-1 / sigma)
  v0 <- 1 + sigma * y / mu
  u0 <- 1 - v0^(-1 / sigma)
  u_s1 <- (0.5 - 1 / 32) * u1 + 1 / 32
  u_s0 <- (0.5 - 1 / 32) * u0 + 1 / 32
  C_beta <- pbeta(0.5, nu, nu) - pbeta(1 / 32, nu, nu)
  G3_1 <- (pbeta(u_s1, nu, nu) - pbeta(1 / 32, nu, nu)) / C_beta
  G3_0 <- (pbeta(u_s0, nu, nu) - pbeta(1 / 32, nu, nu)) / C_beta
  log(G3_1 - G3_0)
}

## DEGPD4 (type 5, 4-param)
.logpmf_degpd4 <- function(y, mu, sigma, nu, tau) {
  v1 <- 1 + sigma * (y + 1) / mu
  u1 <- 1 - v1^(-1 / sigma)
  v0 <- 1 + sigma * y / mu
  u0 <- 1 - v0^(-1 / sigma)
  G4_1 <- 1 - pbeta((1 - u1)^nu, 1 / nu, 2)
  G4_0 <- 1 - pbeta((1 - u0)^nu, 1 / nu, 2)
  G5_1 <- G4_1^(tau / 2)
  G5_0 <- G4_0^(tau / 2)
  log(G5_1 - G5_0)
}

## ---------------------------------------------------------------
## Lazy-init cache for Deriv-generated derivative functions
## ---------------------------------------------------------------

.deriv_cache <- new.env(parent = emptyenv())
.deriv_cache$initialized <- FALSE

.ensure_deriv_cache <- function() {
  if (.deriv_cache$initialized) return(invisible(NULL))

  ## Register custom Deriv rules for pbeta / dbeta
  ## Must capture the drule environment first — Deriv::drule[["x"]] <- val
  ## fails because R treats it as a namespace assignment.
  dr <- Deriv::drule
  dr[["pbeta"]] <- alist(
    q      = dbeta(q, shape1, shape2),
    shape1 = stop("pbeta shape1 deriv not supported"),
    shape2 = stop("pbeta shape2 deriv not supported")
  )
  dr[["dbeta"]] <- alist(
    x      = dbeta(x, shape1, shape2) * ((shape1 - 1) / x - (shape2 - 1) / (1 - x)),
    shape1 = dbeta(x, shape1, shape2) * (log(x)      - digamma(shape1) + digamma(shape1 + shape2)),
    shape2 = dbeta(x, shape1, shape2) * (log(1 - x)  - digamma(shape2) + digamma(shape1 + shape2))
  )

  ## --- Continuous EGPD derivatives ---

  ## EGPD1: fully analytical
  .deriv_cache$dldm_egpd1 <- Deriv::Deriv(.logdens_egpd1, "mu")
  .deriv_cache$dldd_egpd1 <- Deriv::Deriv(.logdens_egpd1, "sigma")
  .deriv_cache$dldv_egpd1 <- Deriv::Deriv(.logdens_egpd1, "nu")

  ## EGPD3: fully analytical
  .deriv_cache$dldm_egpd3 <- Deriv::Deriv(.logdens_egpd3, "mu")
  .deriv_cache$dldd_egpd3 <- Deriv::Deriv(.logdens_egpd3, "sigma")
  .deriv_cache$dldv_egpd3 <- Deriv::Deriv(.logdens_egpd3, "nu")

  ## EGPD5: fully analytical (pnorm handled by Deriv)
  .deriv_cache$dldm_egpd5 <- Deriv::Deriv(.logdens_egpd5, "mu")
  .deriv_cache$dldd_egpd5 <- Deriv::Deriv(.logdens_egpd5, "sigma")
  .deriv_cache$dldv_egpd5 <- Deriv::Deriv(.logdens_egpd5, "nu")

  ## EGPD6: main part analytical; dldv gets numerical correction
  .deriv_cache$dldm_egpd6 <- Deriv::Deriv(.logdens_egpd6_main, "mu")
  .deriv_cache$dldd_egpd6 <- Deriv::Deriv(.logdens_egpd6_main, "sigma")
  .deriv_cache$dldv_egpd6_main <- Deriv::Deriv(.logdens_egpd6_main, "nu")

  ## EGPD4: Deriv part (g4 density + GPD)
  .deriv_cache$dldm_egpd4_d <- Deriv::Deriv(.logdens_egpd4_deriv, "mu")
  .deriv_cache$dldd_egpd4_d <- Deriv::Deriv(.logdens_egpd4_deriv, "sigma")
  .deriv_cache$dldv_egpd4_d <- Deriv::Deriv(.logdens_egpd4_deriv, "nu")
  .deriv_cache$dldt_egpd4_d <- Deriv::Deriv(.logdens_egpd4_deriv, "tau")

  ## EGPD4: G4 term — d/dmu, d/dsigma analytical (pbeta q-rule)
  .deriv_cache$dldm_egpd4_G4 <- Deriv::Deriv(.logG4_term, "mu")
  .deriv_cache$dldd_egpd4_G4 <- Deriv::Deriv(.logG4_term, "sigma")

  ## --- Discrete EGPD derivatives ---

  ## DEGPD1: fully analytical
  .deriv_cache$dldm_degpd1 <- Deriv::Deriv(.logpmf_degpd1, "mu")
  .deriv_cache$dldd_degpd1 <- Deriv::Deriv(.logpmf_degpd1, "sigma")
  .deriv_cache$dldv_degpd1 <- Deriv::Deriv(.logpmf_degpd1, "nu")

  ## DEGPD3: d/d(mu,sigma) analytical; d/dnu numerical fallback
  .deriv_cache$dldm_degpd3 <- Deriv::Deriv(.logpmf_degpd3, "mu")
  .deriv_cache$dldd_degpd3 <- Deriv::Deriv(.logpmf_degpd3, "sigma")

  ## DEGPD5: fully analytical
  .deriv_cache$dldm_degpd5 <- Deriv::Deriv(.logpmf_degpd5, "mu")
  .deriv_cache$dldd_degpd5 <- Deriv::Deriv(.logpmf_degpd5, "sigma")
  .deriv_cache$dldv_degpd5 <- Deriv::Deriv(.logpmf_degpd5, "nu")

  ## DEGPD6: d/d(mu,sigma) analytical; d/dnu numerical fallback
  .deriv_cache$dldm_degpd6 <- Deriv::Deriv(.logpmf_degpd6_main, "mu")
  .deriv_cache$dldd_degpd6 <- Deriv::Deriv(.logpmf_degpd6_main, "sigma")

  ## DEGPD4: d/d(mu,sigma) analytical; d/d(nu,tau) numerical fallback
  .deriv_cache$dldm_degpd4 <- Deriv::Deriv(.logpmf_degpd4, "mu")
  .deriv_cache$dldd_degpd4 <- Deriv::Deriv(.logpmf_degpd4, "sigma")

  .deriv_cache$initialized <- TRUE
  invisible(NULL)
}

## ---------------------------------------------------------------
## Dispatcher: 3-parameter analytical derivatives
## ---------------------------------------------------------------

#' @keywords internal
.egpd_ad3 <- function(y, mu, sigma, nu, family, idx) {
  .ensure_deriv_cache()

  dl <- switch(
    paste0(family, "_", idx),

    ## --- Continuous EGPD ---
    "egpd1_1" = .deriv_cache$dldm_egpd1(y, mu, sigma, nu),
    "egpd1_2" = .deriv_cache$dldd_egpd1(y, mu, sigma, nu),
    "egpd1_3" = .deriv_cache$dldv_egpd1(y, mu, sigma, nu),

    "egpd3_1" = .deriv_cache$dldm_egpd3(y, mu, sigma, nu),
    "egpd3_2" = .deriv_cache$dldd_egpd3(y, mu, sigma, nu),
    "egpd3_3" = .deriv_cache$dldv_egpd3(y, mu, sigma, nu),

    "egpd5_1" = .deriv_cache$dldm_egpd5(y, mu, sigma, nu),
    "egpd5_2" = .deriv_cache$dldd_egpd5(y, mu, sigma, nu),
    "egpd5_3" = .deriv_cache$dldv_egpd5(y, mu, sigma, nu),

    "egpd6_1" = .deriv_cache$dldm_egpd6(y, mu, sigma, nu),
    "egpd6_2" = .deriv_cache$dldd_egpd6(y, mu, sigma, nu),
    "egpd6_3" = .deriv_cache$dldv_egpd6_main(y, mu, sigma, nu) -
                .dlogC_dnu_egpd6(nu),

    ## --- Discrete EGPD ---
    "degpd1_1" = .deriv_cache$dldm_degpd1(y, mu, sigma, nu),
    "degpd1_2" = .deriv_cache$dldd_degpd1(y, mu, sigma, nu),
    "degpd1_3" = {
      dl <- .deriv_cache$dldv_degpd1(y, mu, sigma, nu)
      ## Fix y=0 boundary: u0=0 causes 0^nu*log(0)=NaN in Deriv output.
      ## True value: d/dnu log(u1^nu) = log(u1)
      if (any(y == 0)) {
        v1_0 <- 1 + sigma[y == 0] * 1 / mu[y == 0]
        u1_0 <- 1 - v1_0^(-1 / sigma[y == 0])
        dl[y == 0] <- log(u1_0)
      }
      dl
    },

    "degpd3_1" = .deriv_cache$dldm_degpd3(y, mu, sigma, nu),
    "degpd3_2" = .deriv_cache$dldd_degpd3(y, mu, sigma, nu),
    "degpd3_3" = .egpd_nd3(y, mu, sigma, nu, dDEGPD3, idx = 3L),

    "degpd5_1" = .deriv_cache$dldm_degpd5(y, mu, sigma, nu),
    "degpd5_2" = .deriv_cache$dldd_degpd5(y, mu, sigma, nu),
    "degpd5_3" = .deriv_cache$dldv_degpd5(y, mu, sigma, nu),

    "degpd6_1" = .deriv_cache$dldm_degpd6(y, mu, sigma, nu),
    "degpd6_2" = .deriv_cache$dldd_degpd6(y, mu, sigma, nu),
    "degpd6_3" = .egpd_nd3(y, mu, sigma, nu, dDEGPD6, idx = 3L),

    stop("Unknown family/idx: ", family, "_", idx)
  )
  dl[!is.finite(dl)] <- 0
  dl
}

## ---------------------------------------------------------------
## Dispatcher: 4-parameter analytical derivatives
## ---------------------------------------------------------------

#' @keywords internal
.egpd_ad4 <- function(y, mu, sigma, nu, tau, family, idx) {
  .ensure_deriv_cache()

  dl <- switch(
    paste0(family, "_", idx),

    ## --- Continuous EGPD4 ---
    "egpd4_1" = .deriv_cache$dldm_egpd4_d(y, mu, sigma, nu, tau) +
                .deriv_cache$dldm_egpd4_G4(y, mu, sigma, nu, tau),
    "egpd4_2" = .deriv_cache$dldd_egpd4_d(y, mu, sigma, nu, tau) +
                .deriv_cache$dldd_egpd4_G4(y, mu, sigma, nu, tau),
    "egpd4_3" = {
      deriv_part <- .deriv_cache$dldv_egpd4_d(y, mu, sigma, nu, tau)
      h <- pmax(abs(nu) * 1e-6, 1e-8)
      G4_corr <- (.logG4_term(y, mu, sigma, nu + h, tau) -
                   .logG4_term(y, mu, sigma, nu - h, tau)) / (2 * h)
      deriv_part + G4_corr
    },
    "egpd4_4" = .deriv_cache$dldt_egpd4_d(y, mu, sigma, nu, tau) +
                .ad_dldt_egpd4_G4_fn(y, mu, sigma, nu, tau),

    ## --- Discrete EGPD4 ---
    "degpd4_1" = .deriv_cache$dldm_degpd4(y, mu, sigma, nu, tau),
    "degpd4_2" = .deriv_cache$dldd_degpd4(y, mu, sigma, nu, tau),
    "degpd4_3" = .egpd_nd4(y, mu, sigma, nu, tau, dDEGPD4, idx = 3L),
    "degpd4_4" = .egpd_nd4(y, mu, sigma, nu, tau, dDEGPD4, idx = 4L),

    ## --- Zero-Inflated Continuous EGPD ---
    "ziegpd1_1" = .zi_ad(y, mu, sigma, nu, tau, "egpd1", 1L),
    "ziegpd1_2" = .zi_ad(y, mu, sigma, nu, tau, "egpd1", 2L),
    "ziegpd1_3" = .zi_ad(y, mu, sigma, nu, tau, "egpd1", 3L),
    "ziegpd1_4" = .zi_ad_tau(y, mu, sigma, nu, tau, dEGPD1),

    "ziegpd3_1" = .zi_ad(y, mu, sigma, nu, tau, "egpd3", 1L),
    "ziegpd3_2" = .zi_ad(y, mu, sigma, nu, tau, "egpd3", 2L),
    "ziegpd3_3" = .zi_ad(y, mu, sigma, nu, tau, "egpd3", 3L),
    "ziegpd3_4" = .zi_ad_tau(y, mu, sigma, nu, tau, dEGPD3),

    "ziegpd5_1" = .zi_ad(y, mu, sigma, nu, tau, "egpd5", 1L),
    "ziegpd5_2" = .zi_ad(y, mu, sigma, nu, tau, "egpd5", 2L),
    "ziegpd5_3" = .zi_ad(y, mu, sigma, nu, tau, "egpd5", 3L),
    "ziegpd5_4" = .zi_ad_tau(y, mu, sigma, nu, tau, dEGPD5),

    "ziegpd6_1" = .zi_ad(y, mu, sigma, nu, tau, "egpd6", 1L),
    "ziegpd6_2" = .zi_ad(y, mu, sigma, nu, tau, "egpd6", 2L),
    "ziegpd6_3" = .zi_ad(y, mu, sigma, nu, tau, "egpd6", 3L),
    "ziegpd6_4" = .zi_ad_tau(y, mu, sigma, nu, tau, dEGPD6),

    ## --- Zero-Inflated Discrete EGPD ---
    "zidegpd1_1" = .zid_ad(y, mu, sigma, nu, tau, "degpd1", 1L, dDEGPD1),
    "zidegpd1_2" = .zid_ad(y, mu, sigma, nu, tau, "degpd1", 2L, dDEGPD1),
    "zidegpd1_3" = .zid_ad(y, mu, sigma, nu, tau, "degpd1", 3L, dDEGPD1),
    "zidegpd1_4" = .zid_ad_tau(y, mu, sigma, nu, tau, dDEGPD1),

    "zidegpd3_1" = .zid_ad(y, mu, sigma, nu, tau, "degpd3", 1L, dDEGPD3),
    "zidegpd3_2" = .zid_ad(y, mu, sigma, nu, tau, "degpd3", 2L, dDEGPD3),
    "zidegpd3_3" = .zid_ad(y, mu, sigma, nu, tau, "degpd3", 3L, dDEGPD3),
    "zidegpd3_4" = .zid_ad_tau(y, mu, sigma, nu, tau, dDEGPD3),

    "zidegpd5_1" = .zid_ad(y, mu, sigma, nu, tau, "degpd5", 1L, dDEGPD5),
    "zidegpd5_2" = .zid_ad(y, mu, sigma, nu, tau, "degpd5", 2L, dDEGPD5),
    "zidegpd5_3" = .zid_ad(y, mu, sigma, nu, tau, "degpd5", 3L, dDEGPD5),
    "zidegpd5_4" = .zid_ad_tau(y, mu, sigma, nu, tau, dDEGPD5),

    "zidegpd6_1" = .zid_ad(y, mu, sigma, nu, tau, "degpd6", 1L, dDEGPD6),
    "zidegpd6_2" = .zid_ad(y, mu, sigma, nu, tau, "degpd6", 2L, dDEGPD6),
    "zidegpd6_3" = .zid_ad(y, mu, sigma, nu, tau, "degpd6", 3L, dDEGPD6),
    "zidegpd6_4" = .zid_ad_tau(y, mu, sigma, nu, tau, dDEGPD6),

    stop("Unknown family/idx: ", family, "_", idx)
  )
  dl[!is.finite(dl)] <- 0
  dl
}

## ---------------------------------------------------------------
## Zero-Inflated helpers
## ---------------------------------------------------------------

## ZI continuous: d/d(mu,sigma,nu)
## For y > 0: same as base family
## For y = 0: 0 (only pi contributes)
.zi_ad <- function(y, mu, sigma, nu, tau, base_family, idx) {
  dl_base <- .egpd_ad3(y, mu, sigma, nu, base_family, idx)
  ifelse(y > 0, dl_base, 0)
}

## ZI continuous: d/dtau (= d/dpi)
.zi_ad_tau <- function(y, mu, sigma, nu, tau, dfun_base) {
  ifelse(y > 0, -1 / (1 - tau), 1 / tau)
}

## ZI discrete: d/d(mu,sigma,nu)
## For y > 0: same as base family
## For y = 0: (1-pi)*P_base(0)*dl_base_0 / [pi + (1-pi)*P_base(0)]
.zid_ad <- function(y, mu, sigma, nu, tau, base_family, idx, dfn_base) {
  dl_base <- .egpd_ad3(y, mu, sigma, nu, base_family, idx)
  p0 <- dfn_base(0, mu = mu, sigma = sigma, nu = nu)
  denom <- tau + (1 - tau) * p0

  dl_base_0 <- .egpd_ad3(rep(0, length(y)), mu, sigma, nu, base_family, idx)
  dl_y0 <- (1 - tau) * p0 * dl_base_0 / denom

  dl <- ifelse(y > 0, dl_base, dl_y0)
  dl[!is.finite(dl)] <- 0
  dl
}

## ZI discrete: d/dtau
.zid_ad_tau <- function(y, mu, sigma, nu, tau, dfn_base) {
  p0 <- dfn_base(0, mu = mu, sigma = sigma, nu = nu)
  denom <- tau + (1 - tau) * p0
  ifelse(y > 0, -1 / (1 - tau), (1 - p0) / denom)
}

## ---------------------------------------------------------------
## Override helpers: plug analytical derivatives into a gamlss family
## ---------------------------------------------------------------

## 3-parameter family override
.override_derivs3 <- function(fam, family_key) {
  fam$dldm <- eval(bquote(function(y, mu, sigma, nu, ...) {
    egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L)
  }))
  fam$dldd <- eval(bquote(function(y, mu, sigma, nu, ...) {
    egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L)
  }))
  fam$dldv <- eval(bquote(function(y, mu, sigma, nu, ...) {
    egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
  }))

  fam$d2ldm2 <- eval(bquote(function(y, mu, sigma, nu, ...) {
    dl <- egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldd2 <- eval(bquote(function(y, mu, sigma, nu, ...) {
    dl <- egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldv2 <- eval(bquote(function(y, mu, sigma, nu, ...) {
    dl <- egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))

  fam$d2ldmdd <- eval(bquote(function(y, mu, sigma, nu, ...) {
    -egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L) *
      egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L)
  }))
  fam$d2ldmdv <- eval(bquote(function(y, mu, sigma, nu, ...) {
    -egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L) *
      egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
  }))
  fam$d2ldddv <- eval(bquote(function(y, mu, sigma, nu, ...) {
    -egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L) *
      egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
  }))

  fam
}

## 4-parameter family override
.override_derivs4 <- function(fam, family_key) {
  fam$dldm <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L)
  }))
  fam$dldd <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L)
  }))
  fam$dldv <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
  }))
  fam$dldt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))

  fam$d2ldm2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldd2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldv2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldt2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))

  fam$d2ldmdd <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L)
  }))
  fam$d2ldmdv <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
  }))
  fam$d2ldmdt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))
  fam$d2ldddv <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
  }))
  fam$d2ldddt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))
  fam$d2ldvdt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))

  fam
}
