# Simulation studies

Long-running Monte Carlo studies of the `egpd` fitting functions. They are not
tests: they print bias, RMSE, Wald-interval coverage and standard-error
calibration tables for inspection, and assert nothing. They live outside
`tests/` because `R CMD check` runs every `tests/*.R` file and these scripts
take hours. The directory is excluded from the package build
(`.Rbuildignore`).

| Script | What it does | Default size |
|---|---|---|
| `test_simulation.R` | Parameter recovery for EGPD, DEGPD and ZIDEGPD models 1–6 (18 studies) on the link scale: bias, RMSE, 95% Wald coverage and SE/empirical-SD ratio | `nsim = 100` replicates of `n = 2000` per study |
| `sim_study_degpd1_identity.R` | Recovery of the DEGPD model 1 shape with the identity link (which allows ξ < 0) versus the default log link, for ξ from −0.3 to 0.2 | 15 replicates of `n = 3000` per ξ |

Run them from the package root against the installed package, for example

```sh
Rscript simulations/test_simulation.R 50 1000   # nsim = 50, n = 1000
Rscript simulations/sim_study_degpd1_identity.R
```

or against a source tree by setting `EGPD_SOURCE`, which loads the package
with `pkgload::load_all()`:

```sh
EGPD_SOURCE=. Rscript simulations/test_simulation.R 10 500
```

The Monte Carlo results shown in the *Parameter Estimate Coverage
Assessment* vignette are produced separately by
`vignettes/precompute/parameter-coverage.R`.
