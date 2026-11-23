# egpd: Extended Generalized Pareto Distribution GAMs

<!-- badges: start -->
[![R-CMD-check](https://github.com/sdwfrost/egpd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sdwfrost/egpd/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The **egpd** package fits Extended Generalized Pareto Distribution (EGPD),
Discrete EGPD (DEGPD), and Zero-Inflated Discrete EGPD (ZIDEGPD) models
within a GAM (Generalized Additive Model) framework. It also provides
zero-inflated continuous EGPD distribution helpers plus neural-Bayes
estimators for BEGPD, BDEGPD, BZIDEGPD, MDGPD, and ZIMDGPD models, with
additional support for fitting via
[gamlss](https://CRAN.R-project.org/package=gamlss) and
[bamlss](https://CRAN.R-project.org/package=bamlss). The fitting
infrastructure is adapted from the
[evgam](https://CRAN.R-project.org/package=evgam) package by Ben Youngman.
The package also provides univariate distribution fitting (`distfit()`),
compound Poisson EGPD models for aggregated data, blended generalized Pareto
and composite Pareto distributions, and overdispersed count families (Bell,
Poisson-inverse Gaussian, generalised Poisson-inverse Gaussian and others)
that can be fitted with the same GAM machinery.

Models 1--6 are supported for each distribution family. Models 1--4
follow the parameterizations of Naveau et al. (2016). Models 5--6 use the
truncated normal and truncated beta transformations of Ahmad & Arshad
(2026), which satisfy the conditions of Gamet & Jalbert (2022). The
discrete models follow Ahmad, Gaetan & Naveau (2025), and the zero-inflated
models follow Ahmad & Arshad (2026) (discrete) and Abbas, Ahmad & Ahmad
(2025) (continuous). Each model extends the standard GPD by composing it
with a transformation function *G*:

| Model | G-function | Parameters | Description |
|-------|------------|------------|-------------|
| 1 | Power: *G*(*u*) = *u*^*κ* | *σ*, *ξ*, *κ* | Simplest; reduces to standard GPD when *κ* = 1 |
| 2 | Mixture: *G*(*u*) = *p* *u*^*κ*₁ + (1−*p*) *u*^*κ*₂ | *σ*, *ξ*, *κ*₁, Δ*κ*, *p* | Mixture of two power transformations |
| 3 | Incomplete beta | *σ*, *ξ*, *δ* | Beta CDF transformation; more flexible bulk |
| 4 | Power-beta | *σ*, *ξ*, *δ*, *κ* | Combines power and beta; most flexible |
| 5 | Truncated normal | *σ*, *ξ*, *κ* | Truncated normal CDF transformation |
| 6 | Truncated beta | *σ*, *ξ*, *κ* | Truncated beta CDF transformation |

Here *σ* is the GPD scale, *ξ* is the GPD shape, and the remaining
parameters control the *G*-function.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("sdwfrost/egpd")
```

## Quick start

```r
library(egpd)

# Fit a DEGPD model to insurance complaint counts
data(ny_complaints)
fit <- egpd(
  list(upheld ~ s(year), ~ 1, ~ 1),
  data = ny_complaints,
  family = "degpd",
  degpd.args = list(m = 1)
)
summary(fit)
```

## Aggregated data

The compound Poisson EGPD (`dcpegpd()`, `pcpegpd()`, `qcpegpd()`,
`rcpegpd()`) describes the total of a Poisson number of EGPD events, each
recorded on a grid of width `h` as by a rain gauge (default `h = 0.2`;
Ailliot, Gaetan & Naveau, 2026). `rcpegpd()` simulates the same
gauge-discretised law as the d/p/q functions; `rcpegpd(h = 0)` gives the
continuous sum. `fit_aggregate_egpd()` fits the multiscale model across
aggregation durations by composite likelihood, and
`bootstrap_aggregate_egpd()` gives moving-block bootstrap standard errors
and percentile intervals that account for serial dependence.

## Bivariate EGPD (BEGPD)

The package also supports bivariate EGPD (BEGPD) fitting
via neural Bayes estimation, using trained neural networks. This
requires [Julia](https://julialang.org/) (>= 1.11) and additional
dependencies.

Neural models must be retrained for the corrected simulation and parameter
contracts. Legacy unversioned assets are rejected rather than silently reused;
keep the metadata sidecars generated alongside new models. Julia training and
calibration have not been verified by the R-only regression checks.

NPE posterior draws come from adaptive importance sampling on the learned
likelihood-to-evidence ratio; the effective sample size (ESS) and the number
of distinct draws are reported, with a warning when the ESS is small.

### Additional dependencies

```r
# R packages
install.packages("JuliaConnectoR")
remotes::install_github("msainsburydale/NeuralEstimators")
```

Then install the required Julia packages (training with the R package
NeuralEstimators 0.2.2 needs NeuralEstimators.jl 0.3 or later; see
`?train_begpd`):

```julia
using Pkg
Pkg.add(["NeuralEstimators", "Flux", "BSON", "Optimisers"])
```

### Quick start (bivariate BEGPD)

```r
library(egpd)

# Simulate bivariate BEGPD data (no Julia needed)
Y <- rbegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)

# Train compatible assets and fit (requires Julia)
paths <- train_begpd(savepath = "my_models", quick = TRUE)
fit <- distfit(Y, family = "begpd", method = "neuralbayes",
               model.path = paths$npe, estimator = "npe")
summary(fit)
plot(fit)
```

### Discrete extensions (Experimental)

The package also provides **experimental** bivariate discrete EGPD
(`family="bdegpd"`) and zero-inflated bivariate discrete EGPD
(`family="bzidegpd"`) distributions, constructed by applying `floor()` to
continuous BEGPD samples. These use the same neural Bayes estimation framework;
see the corresponding vignettes for training compatible assets before fitting.

```r
# Simulate bivariate discrete EGPD
Y <- rbdegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)

# Zero-inflated version
Y_zi <- rbzidegpd(1000, kappa = 2, sigma = 1, xi = 0.1,
                   thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
```

### MDGPD: Multivariate Discrete GPD (Experimental)

The package also includes the Aka-Kratz-Naveau (2025) Multivariate Discrete
Generalized Pareto Distribution (`family="mdgpd"`), using the paper's
threshold-excess construction. These are signed integer excesses, with at least
one positive component, not raw nonnegative counts. The standardized maximum
is positive geometric. The identifiable Poisson-difference parameter is
`rate = (1 - rho) * lambda`; the redundant legacy `lambda` and `rho` arguments
now error.

```r
# Simulate bivariate MDGPD
Y <- rmdgpd(1000, sigma = 2, xi = 0.2, rate = 0.5)
```

For neural fitting, retrain with `train_mdgpd()` and refit using the new assets;
old fits cannot be relabeled as the corrected model. Neural MDGPD estimation is
homogeneous only: the trained estimator assumes a common `sigma`, `xi` and
`rate` for all margins.

For details, see Aka, Kratz & Naveau (2025),
[arXiv:2506.19361](https://arxiv.org/abs/2506.19361).

For details on the BEGPD, see Alotaibi, Sainsbury-Dale, Naveau, Gaetan & Huser (2025),
[arXiv:2509.05982](https://arxiv.org/abs/2509.05982).

## Vignettes

Rendered versions are on the [package website](https://sdwfrost.github.io/egpd/articles/);
the sources are in [`vignettes/`](vignettes/). After installing with
`remotes::install_github("sdwfrost/egpd", build_vignettes = TRUE)`, run
`browseVignettes("egpd")`.

Univariate EGPD and DEGPD GAMs:

- [Discrete EGPD Models for Insurance Complaint Counts](https://sdwfrost.github.io/egpd/articles/insurance-complaints.html)
- [Zero-Inflated Discrete EGPD Models for Doctor Visit Counts](https://sdwfrost.github.io/egpd/articles/doctor-visits.html)
- [Threshold Exceedance Modeling with DEGPD](https://sdwfrost.github.io/egpd/articles/gaming-offenses.html)
- [Continuous EGPD Models for Temperature Extremes](https://sdwfrost.github.io/egpd/articles/temperature-extremes.html)
- [Simulation Examples](https://sdwfrost.github.io/egpd/articles/simulation.html)

Univariate distribution fitting:

- [Fitting Continuous Distributions with distfit](https://sdwfrost.github.io/egpd/articles/distfit-continuous.html)
- [Fitting Discrete Distributions with distfit](https://sdwfrost.github.io/egpd/articles/distfit-discrete.html)
- [Blended Generalized Pareto Distributions](https://sdwfrost.github.io/egpd/articles/blended-gp.html)
- [Composite Pareto Distributions and GAMs](https://sdwfrost.github.io/egpd/articles/composite-pareto.html)
- [Multiscale Aggregated EGPD Models](https://sdwfrost.github.io/egpd/articles/aggregated-egpd.html)

Overdispersed count families:

- [Bell Distribution Models for Overdispersed Counts](https://sdwfrost.github.io/egpd/articles/bell-counts.html)
- [Poisson-Inverse Gaussian Models for Overdispersed Counts](https://sdwfrost.github.io/egpd/articles/pig-counts.html)
- [Generalised Poisson-Inverse Gaussian Models for Heavy-Tailed Counts](https://sdwfrost.github.io/egpd/articles/gpig-counts.html)

Multivariate models (neural Bayes estimation):

- [Bivariate BEGPD via Neural Bayes Estimation](https://sdwfrost.github.io/egpd/articles/multivariate-egpd.html)
- [Bivariate Discrete EGPD via Neural Bayes Estimation](https://sdwfrost.github.io/egpd/articles/bivariate-discrete-egpd.html)
- [Multivariate Discrete Generalized Pareto Distributions](https://sdwfrost.github.io/egpd/articles/mdgpd.html)

Comparisons and calibration:

- [Comparing egpd and bamlss Fits](https://sdwfrost.github.io/egpd/articles/bamlss-comparison.html)
- [Comparing egpd and gamlss Fits](https://sdwfrost.github.io/egpd/articles/gamlss-comparison.html)
- [Predictive Coverage Assessment](https://sdwfrost.github.io/egpd/articles/predictive-coverage.html)
- [Parameter Estimate Coverage Assessment](https://sdwfrost.github.io/egpd/articles/parameter-coverage.html)

## References

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *Communications in Statistics - Simulation and
Computation*, 1--17. <https://doi.org/10.1080/03610918.2025.2585398>

Ahmad, T. and Arshad, I. A. (2026). New flexible versions of the extended
generalized Pareto distribution for count data. *Journal of Applied
Statistics*, 1--21. <https://doi.org/10.1080/02664763.2026.2634785>
(preprint: arXiv:2409.18719, <https://arxiv.org/abs/2409.18719>)

Ahmad, T., Gaetan, C. and Naveau, P. (2025). An extended generalized Pareto
regression model for count data. *Statistical Modelling*, 25(5), 416--431.
<https://doi.org/10.1177/1471082X241266729>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>

Ailliot, P., Gaetan, C. and Naveau, P. (2026). A parsimonious tail compliant
multiscale statistical model for aggregated rainfall. *Advances in Water
Resources*, 208, 105216. <https://doi.org/10.1016/j.advwatres.2026.105216>

Aka, S., Kratz, M. and Naveau, P. (2025). Multivariate discrete generalized
Pareto distributions: theory, simulation, and applications to dry spells.
*arXiv preprint* arXiv:2506.19361. <https://arxiv.org/abs/2506.19361>

Alotaibi, N., Sainsbury-Dale, M., Naveau, P., Gaetan, C. and Huser, R. (2025).
Joint modeling of low and high extremes using a multivariate extended
generalized Pareto distribution. *arXiv preprint* arXiv:2509.05982.
<https://arxiv.org/abs/2509.05982>

Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto
distribution for tail estimation. *Environmetrics*, 33(6), e2744.
<https://doi.org/10.1002/env.2744>

Le Carrer, N. and Gaetan, C. (2022). Distributional regression models for
extended generalized Pareto distributions. *arXiv preprint* arXiv:2209.04660.
<https://doi.org/10.48550/arXiv.2209.04660>

Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016). Modeling
jointly low, moderate, and heavy rainfall intensities without a
threshold selection. *Water Resources Research*, 52(4), 2753--2769.
<https://doi.org/10.1002/2015WR018552>
