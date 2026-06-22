# Analytic derivatives for DEGPD models 5 & 6 — investigation + implementation

**Status: IMPLEMENTED — model 5 and model 6 are now FULLY analytic in every family**
(DEGPD, ZIDEGPD, EGPD). The model-6 κ-derivatives use a stable quadrature for
`∂I(x;κ,κ)/∂κ` (below); no finite differences remain. Verified vs FD oracles:

| family | model 5 (grad / Hess) | model 6 (grad / Hess) |
|---|---|---|
| DEGPD   | 1.1e-9 / 7.9e-10 | 2.9e-9 / 3.2e-9 |
| ZIDEGPD | 5.0e-9 / 4.6e-10 | 7.6e-9 / 3.2e-9 |
| EGPD    | 4.4e-9 / 6.5e-9  | 3.3e-9 / 6.5e-9 |

All 626 tests pass; all fit end-to-end with finite logLik/AIC. (Earlier intermediate
versions used a σ/ξ-analytic + κ-FD *hybrid* for model 6, ~1e-5 Hessian; the quadrature
replaced that with the fully-analytic κ derivatives.)

**Original question.** `degpd5d12` (truncated-normal G) and `degpd6d12` (truncated-beta
G) computed their gradient and Hessian by **central finite differences** of the
per-observation NLL (`nll5`/`nll6`), unlike models 1–4 which use analytic `ee`/`eee`
blocks. Can they be made analytic? (This surfaced because the bounded-ξ FD oracle
showed m5/m6 Hessians match only to ~1e-4, vs ~1e-9 for m1–4.)

## Why it matters

The internally finite-differenced Hessian is only ~1e-4 accurate. That degrades exactly
the things the Hessian feeds: **REML smoothing-parameter estimation** (`.reml*` uses the
penalised Hessian) and **standard errors** (`Vp`/`Vc` invert it). Analytic derivatives
would give machine-precision Hessians, faster fits (no 4-point stencils per Hessian
entry), and better-conditioned REML — the same quality models 1–4 already enjoy.

## Structure (shared)

Every DEGPD model is `nll = -log(G(F_hi) - G(F_lo))`, where

- `F_lo, F_hi` are the **GPD CDF** at `y`, `y+1`. Its derivatives w.r.t. `lsigma`, `lxi`
  are **already derived analytically** in models 1/3/4 (the GPD part is shared). Closed
  form (verified in the prototypes):
  - `dF/dlsigma = -(v/xi) * t^(-1/xi - 1)`
  - `dF/dlxi    = -t^(-1/xi) * (log(t)/xi - v/(xi*t))`,  with `t = 1 + xi*m/sigma`, `v = xi*m/sigma`.
- `G(u; kappa)` is the **model-specific transform**. So the only new pieces are the
  G-function derivatives: `G'(u) = dG/du` (drives the σ, ξ directions via the chain rule
  `dG(F)/dθ = G'(F)·dF/dθ`) and `∂G/∂kappa` (the κ direction).

## Model 5 — truncated normal: **fully analytic, feasible** ✅

`G(u) = (Φ(√κ·(u−1)) − Φ(−√κ)) / (0.5 − Φ(−√κ))`.

- `G'(u) = √κ · φ(√κ·(u−1)) / D`,  `D = 0.5 − Φ(−√κ)`.
- `∂G/∂κ`: with `sk = √κ`, `dsk/dκ = 1/(2sk)`, `z = sk(u−1)`,
  `∂G/∂κ = [ (φ(z)(u−1)/(2sk) + φ(−sk)/(2sk))·D − (Φ(z)−Fmin)·(φ(−sk)/(2sk)) ] / D²`,
  and `∂G/∂lkappa = ∂G/∂κ · κ`.
- Everything reduces to the normal pdf/cdf (`φ`, `Φ`); second derivatives use
  `φ'(x) = −x·φ(x)` — also closed-form.

**Verified:** an R prototype of the analytic gradient matches numerical differences of
`nll5` to **1.4e-7** over 40 (y, θ) cases (`/tmp/m5analytic.R`, reproduced in the session).
→ A full analytic `degpd5d12` (gradient + 6-entry Hessian) is straightforward, just
tedious. Recommended.

## Model 6 — truncated beta: **σ, ξ analytic; κ is the obstacle** ⚠️

`G(u) = (I(u;κ,κ) − I(c2;κ,κ)) / (I(0.5;κ,κ) − I(c2;κ,κ))`, `I` = regularized incomplete
beta (`pbeta`), `u = c1·F + c2`, `c1 = 0.5 − 1/32`, `c2 = 1/32`.

- `G'(u) = c1 · dbeta(u;κ,κ) / D` — **closed form** (the beta density). So the σ, ξ
  derivatives are analytic. **Verified:** analytic σ/ξ gradient matches FD to **2.9e-7**
  (`/tmp/m6analytic.R`).
- `∂G/∂κ` needs `∂/∂κ I(x;κ,κ)` — the **parameter-derivative of the regularized
  incomplete beta**, which has no elementary closed form. **Now implemented** via the
  identities
  - `∂I/∂κ  = K₁ − I·β₁`
  - `∂²I/∂κ² = K₂ − 2K₁·β₁ + I·(β₁² − β₁′)`,

  where `Kₘ = ∫₀ˣ dbeta(t;κ,κ)·[ln(t(1−t))]ᵐ dt`, `β₁ = 2(ψ(κ)−ψ(2κ))`,
  `β₁′ = 2(ψ′(κ)−2ψ′(2κ))` (digamma/trigamma). The integrals `K₁,K₂` are computed by
  the substitution `t = x·z^(1/κ)` (which removes the `t^(κ-1)` endpoint singularity)
  followed by **adaptive Gauss-Kronrod 7-15** quadrature with relative-tolerance +
  minimum-width stopping. Both `K₁` and `K₂` share one adaptive subdivision. Validated:
  `∂I/∂κ ~3e-9`, `∂²I/∂κ²` essentially exact vs `pbeta` finite differences across
  `κ ∈ [0.25, 50]`, `x ∈ [1/32, 0.5]`; the full G derivatives `G_κ, G_κκ, G_Fκ` match FD
  to ~1e-9/1e-8.

All of this lives in `src/gfunc_derivs.h` (`incbeta_kk_derivs`, `tb_kappa`, `tb_Gtheta`)
and is shared by DEGPD, ZIDEGPD, and EGPD model 6.

## Scope note

The continuous **EGPD** and **ZIDEGPD** families have the identical truncated-normal /
truncated-beta G-functions (`egpd5/6`, `zidegpd5/6`), so the same verdict and formulas
carry over: m5 fully analytic, m6 σ/ξ analytic + κ special-cased.

## What was implemented (`src/degpd.cpp`)

- **`degpd5d12` — fully analytic.** Replaced the FD block with the closed-form gradient
  and Hessian: GPD-CDF derivatives via `w = t^(-1/xi)` (first and second order in
  lsigma/lxi), truncated-normal G derivatives `G_F, G_FF, G_κ, G_κκ, G_Fκ` via
  `φ`/`Φ`, assembled by the chain rule `nll = -log(G_hi - G_lo)`.
- **`degpd6d12` — fully analytic** (via `tb_Gtheta` in `src/gfunc_derivs.h`). The σ/ξ
  block uses the beta density (`G_F = c1·dbeta(u;κ,κ)/D`,
  `G_FF = c1²·(κ-1)·dbeta·(1-2u)/(u(1-u))/D`); the κ derivatives `G_κ, G_κκ, G_Fκ` use
  `incbeta_kk_derivs` (adaptive Gauss-Kronrod quadrature for `∂pbeta(x;κ,κ)/∂κ`). Same
  for `zidegpd6d12`; `egpd6d12` uses the same quadrature for the normalizer term.
- Both keep the bounded-ξ forward map (`bounded_lxi`), so the R chain rule
  (`.bounded_xi_chain`) still applies on top — `xi.max` works for m5/m6 unchanged.

## Follow-ups — ZIDEGPD done; EGPD already done upstream

- **ZIDEGPD models 5 & 6 — IMPLEMENTED this round.** `zidegpd5d12` is now **fully
  analytic** (grad ~5e-9, Hessian ~5e-10, both y>0 and y=0 branches) and `zidegpd6d12`
  is a **hybrid** (σ/ξ/π exact ~7e-9; κ via 1-D FD ~9e-6). They reuse the shared header
  `src/gfunc_derivs.h` (`gpd_Fderiv`, `tn_kappa`/`tn_Gtheta`, `tb_GthetaSX`). The
  zero-inflation enters in closed form: for y>0, `nll = -log(1-π) - log(P)` so the π
  gradient/Hessian are `π` and `π(1-π)` with zero σ/ξ/κ cross-terms; for y=0,
  `nll = -log(π + (1-π)G(F(1)))` assembled from the single-point G derivatives.
- **EGPD (continuous) models 5 & 6 were already analytic/hybrid upstream** — `egpd5d12`
  is fully analytic (verified ~6e-9) and `egpd6d12` is already a σ/ξ-analytic + κ-FD
  hybrid (verified ~2e-5). No change was needed; an earlier note here wrongly listed them
  as FD. (EGPD's GPD part uses the identity ξ link, handled inline there.)
- **Fully analytic model 6 κ-derivative — DONE for all three families** (DEGPD,
  ZIDEGPD, EGPD) via the `∂I(x;κ,κ)/∂κ` quadrature in `src/gfunc_derivs.h`
  (`incbeta_kk_derivs`). `degpd6d12`/`zidegpd6d12` build the full G derivatives from
  `tb_Gtheta`; `egpd6d12` (continuous) only needs the normalizer's κ-derivatives at the
  fixed bounds `{1/32, 0.5}`, so it calls `incbeta_kk_derivs` there. Hessians dropped
  from the ~1e-5 hybrid level to ~3e-9 (DEGPD/ZIDEGPD) and ~6e-9 (EGPD).
- *Performance note:* the adaptive quadrature runs per observation in `d12` (DEGPD/
  ZIDEGPD call it at the two/one variable G-points; EGPD only at the two fixed bounds).
  Cost is modest (the full test suite runs in ~22 s) but non-zero; if a fit is
  quadrature-bound, the σ/ξ-analytic + κ-FD hybrid in git history is a faster fallback.
