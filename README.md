# egpd: Extended Generalized Pareto Distribution GAMs

The **egpd** package fits Extended Generalized Pareto Distribution (EGPD),
Discrete EGPD (DEGPD), and Zero-Inflated Discrete EGPD (ZIDEGPD) models
within a GAM (Generalized Additive Model) framework, with additional
support for fitting via
[gamlss](https://CRAN.R-project.org/package=gamlss) and
[bamlss](https://CRAN.R-project.org/package=bamlss). The fitting
infrastructure is adapted from the
[evgam](https://CRAN.R-project.org/package=evgam) package by Ben Youngman.

Models 1--6 are supported for each distribution family. Models 1--4
follow the parameterizations of Naveau et al. (2016), and Models 5--6
follow Gamet & Jalbert (2022). Each model extends the standard GPD by
composing it with a transformation function *G*:

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

## Bivariate EGPD (BEGPD)

The package also supports bivariate Bivariate EGPD (BEGPD) fitting
via neural Bayes estimation, using pre-trained neural networks. This
requires [Julia](https://julialang.org/) (>= 1.11) and additional
dependencies.

### Additional dependencies

```r
# R packages
install.packages("JuliaConnectoR")
remotes::install_github("msainsburydale/NeuralEstimators")
```

Then install the required Julia packages:

```julia
using Pkg
Pkg.add(["NeuralEstimators", "Flux"])
```

### Quick start (bivariate BEGPD)

```r
library(egpd)

# Simulate bivariate BEGPD data (no Julia needed)
Y <- rbegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)

# Fit using neural posterior estimation (requires Julia)
fit <- fitegpd(Y, family = "begpd", method = "neuralbayes", estimator = "npe")
summary(fit)
plot(fit)

# Train your own model (optional)
paths <- train_begpd(savepath = "my_models", quick = TRUE)
fit2 <- fitegpd(Y, family = "begpd", method = "neuralbayes",
                model.path = paths$npe, estimator = "npe")
```

For details, see Alotaibi, Sainsbury-Dale, Naveau, Gaetan & Huser (2025),
[arXiv:2509.05982](https://arxiv.org/abs/2509.05982).

## Vignettes

- [Discrete EGPD Models for Insurance Complaint Counts](articles/insurance-complaints.md)
- [Zero-Inflated Discrete EGPD Models for Doctor Visit Counts](articles/doctor-visits.md)
- [Threshold Exceedance Modeling with DEGPD](articles/gaming-offenses.md)
- [Continuous EGPD Models for Temperature Extremes](articles/temperature-extremes.md)
- [Fitting Continuous Distributions with fitegpd](articles/fitegpd-continuous.md)
- [Fitting Discrete Distributions with fitegpd](articles/fitegpd-discrete.md)
- [Bivariate BEGPD via Neural Bayes Estimation](articles/multivariate-egpd.md)
- [Simulation Examples](articles/simulation.md)
- [Comparing egpd and bamlss Fits](articles/bamlss-comparison.md)
- [Comparing egpd and gamlss Fits](articles/gamlss-comparison.md)
- [Predictive Coverage Assessment](articles/predictive-coverage.md)
- [Parameter Estimate Coverage Assessment](articles/parameter-coverage.md)

## References

Alotaibi, N., Sainsbury-Dale, M., Naveau, P., Gaetan, C., and Huser, R. (2025).
Joint modeling of low and high extremes using a multivariate extended generalized
Pareto distribution. *arXiv preprint* arXiv:2509.05982.
<https://arxiv.org/abs/2509.05982>

Ailliot, P., Gaetan, C., & Naveau, P. (2026). A parsimonious tail compliant multiscale statistical 
model for aggregated rainfall. *Advances in Water Resources*, 208, 105216.<https://doi.org/10.1016/j.advwatres.2026.105216>

Abbas, A., Ahmad, T. and Ahmad, I. (2025). Modeling zero-inflated
precipitation extremes. *Communications in Statistics-Simulation and Computation,* 1-17.
<https://doi.org/10.1080/03610918.2025.2585398>

Ahmad, T. and Arshad, I. A. (2026). New flexible versions of the extended
generalized Pareto model for count data. *Journal of Applied Statistics*. <https://arxiv.org/abs/2409.18719>

Ahmad, T. and Hussain, A. (2025). Flexible model for varying levels of
zeros and outliers in count data. *arXiv preprint* arXiv:2510.27365.
<https://arxiv.org/abs/2510.27365>

Ahmad, T., Gaetan, C., & Naveau, P. (2025). An extended generalized Pareto
regression model for count data. *Statistical Modelling*, 25(5), 416-431.
<https://doi.org/10.1177/1471082X241266729>

Carrer, N. L., & Gaetan, C. (2022). Distributional regression models for Extended Generalized 
Pareto distributions. *arXiv preprint* arXiv:2209.04660.<https://doi.org/10.48550/arXiv.2209.04660>

Gamet, P. and Jalbert, J. (2022). A flexible extended generalized Pareto distribution for tail estimation.
*Environmetrics*, 33(4), e2714.<https://doi.org/10.1002/env.2744>

Naveau, P., Huser, R., Ribereau, P., and Hannart, A. (2016). Modeling
jointly low, moderate, and heavy rainfall intensities without a
threshold selection. *Water Resources Research*, 52(4), 2897--2911.
