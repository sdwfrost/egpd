Multivariate Discrete GPD (MDGPD) via Neural Bayes Estimation
================

The `egpd` package provides **experimental** support for the
Multivariate Discrete Generalized Pareto Distribution (MDGPD) of Aka,
Kratz & Naveau (2025). Unlike the BDEGPD (which discretizes continuous
BEGPD samples via `floor()`), the MDGPD is constructed directly in the
discrete domain using a bivariate Poisson generator and geometric
maximum, yielding theoretically rigorous discrete GPD marginals with
threshold stability.

This vignette covers:

1.  **The MDGPD construction** – from Poisson generator to discrete GPD
2.  **Simulating multivariate data** – `rmdgpd()` and `rzimdgpd()` for
    $d \geq 2$ dimensions
3.  **Exploring parameters and dependence** – the role of $\rho$,
    $\lambda$, $\sigma$, and $\xi$
4.  **Higher dimensions** – trivariate and beyond
5.  **Neural Bayes estimation** – fitting via `fitegpd()`
6.  **Training custom models** – `train_mdgpd()`

## 1. The MDGPD construction

### Standard MDGPD (geometric marginals)

The Aka-Kratz-Naveau construction builds a $d$-dimensional discrete GPD
from three components:

**Step 1: Equicorrelated Poisson generator.** Generate $d$ correlated
Poisson random variables via a common random effect:

$$T_j = X_j + Z, \quad j = 1, \ldots, d$$

where $Z \sim \mathrm{Poisson}(\rho\lambda)$ is shared across all
components and $X_j \sim \mathrm{Poisson}((1-\rho)\lambda)$ are
independent. The correlation between any pair $(T_i, T_j)$ is $\rho$.

**Step 2: Spectral differences.** For each component $i$, compute:

$$\Delta_i = T_i - \max_{j \neq i} T_j$$

This captures how much component $i$ exceeds the maximum of all others.
Note that $\Delta_i \leq 0$ when component $i$ is not the overall
maximum.

**Step 3: Geometric maximum.** Generate
$G \sim \mathrm{Geometric}(1 - e^{-1})$, independently of the Poisson
generator.

**Step 4: Standard MDGPD.** Combine:

$$N_i = G + \min(\Delta_i, 0), \quad i = 1, \ldots, d$$

The resulting vector $(N_1, \ldots, N_d)$ has geometric marginals
$\mathrm{Geom}(1 - e^{-1})$ and the property that
$\max(N_1, \ldots, N_d) = G$.

### Non-standard MDGPD (discrete GPD marginals)

To obtain discrete GPD marginals with parameters $\sigma > 0$ and
$\xi \geq 0$, apply the quantile transform:

$$M_i = \left\lfloor \frac{\sigma}{\xi}\left(e^{\xi \max(N_i, 0)} - 1\right) \right\rfloor$$

When $\xi = 0$, this reduces to $M_i = \lfloor \sigma N_i \rfloor$
(scaled geometric marginals).

### Parameter summary

| Parameter | Symbol    | Range         | Role                                  |
|-----------|-----------|---------------|---------------------------------------|
| `sigma`   | $\sigma$  | $(0, \infty)$ | GPD scale (marginal spread)           |
| `xi`      | $\xi$     | $[0, \infty)$ | GPD shape (tail heaviness)            |
| `lambda`  | $\lambda$ | $(0, \infty)$ | Poisson rate (spectral spread)        |
| `rho`     | $\rho$    | $[0, 1)$      | Equicorrelation (dependence strength) |
| `pi0`     | $\pi_0$   | $(0, 1)$      | Joint zero-inflation (ZIMDGPD only)   |

## 2. Simulating data

### Bivariate MDGPD (default)

The `rmdgpd()` function generates samples from the $d$-dimensional
MDGPD. It is pure R and does not require Julia.

``` r
library(egpd)
set.seed(42)

Y <- rmdgpd(2000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
head(Y)
#>      Y1 Y2
#> [1,]  0  0
#> [2,]  0  0
#> [3,]  2  0
#> [4,]  0  0
#> [5,]  2  0
#> [6,]  0  0
cat("Dimensions:", nrow(Y), "x", ncol(Y), "\n")
#> Dimensions: 2000 x 2
cat("Storage mode:", storage.mode(Y), "\n")
#> Storage mode: integer
cat("Range Y1:", range(Y[, 1]), "  Range Y2:", range(Y[, 2]), "\n")
#> Range Y1: 0 39   Range Y2: 0 30
```

``` r
plot(jitter(Y[, 1]), jitter(Y[, 2]), pch = 20, cex = 0.3,
     xlab = expression(Y[1]), ylab = expression(Y[2]),
     main = "Simulated bivariate MDGPD (jittered)",
     col = adjustcolor("steelblue", 0.4))
```

![](mdgpd_files/figure-gfm/scatter-basic-1.png)<!-- -->

### Marginal distributions

The marginals of the non-standard MDGPD are discrete GPD. When
$\sigma = 1$ and $\xi = 0$, they reduce to $\mathrm{Geom}(1 - e^{-1})$:

``` r
set.seed(42)
Y_geom <- rmdgpd(5000, sigma = 1, xi = 0, lambda = 1, rho = 0.5)

op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
barplot(table(Y_geom[, 1]) / nrow(Y_geom),
        main = expression("Marginal " * Y[1] * " (sigma=1, xi=0)"),
        xlab = expression(Y[1]), ylab = "Proportion",
        col = "lightblue", border = "grey")

## Overlay theoretical Geom(1-e^{-1})
p_geom <- 1 - exp(-1)
k_vals <- 0:max(Y_geom[, 1])
points(seq_along(k_vals), dgeom(k_vals, prob = p_geom),
       pch = 16, col = "red", cex = 0.8)
legend("topright", "Geom(1-e^{-1})", pch = 16, col = "red",
       cex = 0.8, bg = "white")

## With sigma=2, xi=0.2: heavier tails
barplot(table(Y[, 1]) / nrow(Y),
        main = expression("Marginal " * Y[1] * " (sigma=2, xi=0.2)"),
        xlab = expression(Y[1]), ylab = "Proportion",
        col = "lightblue", border = "grey")
```

![](mdgpd_files/figure-gfm/marginals-1.png)<!-- -->

``` r
par(op)
```

### Maximum component is geometric

A key theoretical property: the component-wise maximum
$\max(N_1, \ldots, N_d) = G \sim \mathrm{Geom}(1-e^{-1})$. We verify
this for the standard MDGPD ($\sigma=1$, $\xi=0$):

``` r
maxY <- pmax(Y_geom[, 1], Y_geom[, 2])

barplot(table(maxY) / length(maxY),
        main = expression("max(" * Y[1] * ", " * Y[2] * ") vs Geom(1-" * e^{-1} * ")"),
        xlab = "max(Y1, Y2)", ylab = "Proportion",
        col = "lightblue", border = "grey")
k_max <- 0:max(maxY)
points(seq_along(k_max), dgeom(k_max, prob = p_geom),
       pch = 16, col = "red", cex = 0.8)
legend("topright", "Geom(1-e^{-1})", pch = 16, col = "red",
       cex = 0.8, bg = "white")
```

![](mdgpd_files/figure-gfm/max-geom-1.png)<!-- -->

``` r

cat("Empirical mean of max:", round(mean(maxY), 3), "\n")
#> Empirical mean of max: 0.578
cat("Theoretical mean:", round(exp(-1) / (1 - exp(-1)), 3), "\n")
#> Theoretical mean: 0.582
```

### Zero-inflated MDGPD

``` r
set.seed(42)
Y_zi <- rzimdgpd(2000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5,
                  pi0 = 0.3)

joint_zeros <- mean(rowSums(Y_zi) == 0)
cat("Proportion of all-zero rows:", round(joint_zeros, 3), "\n")
#> Proportion of all-zero rows: 0.754
cat("(Expected >= pi0 = 0.3 due to natural zeros)\n")
#> (Expected >= pi0 = 0.3 due to natural zeros)
```

## 3. The role of parameters

### Dependence strength ($\rho$)

The parameter $\rho$ directly controls the equicorrelation of the
Poisson generator. Higher $\rho$ means the common component $Z$
dominates, making all $T_j$ nearly identical and thus producing
near-identical $(Y_1, \ldots, Y_d)$:

``` r
set.seed(42)
n_demo <- 3000

rho_vals <- c(0, 0.3, 0.7, 0.95)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (rho_v in rho_vals) {
  Y_rho <- rmdgpd(n_demo, sigma = 2, xi = 0.2, lambda = 1, rho = rho_v)
  cor_val <- cor(Y_rho[, 1], Y_rho[, 2])
  plot(jitter(Y_rho[, 1]), jitter(Y_rho[, 2]), pch = 20, cex = 0.3,
       main = bquote(rho == .(rho_v) ~ " (cor = " * .(round(cor_val, 2)) * ")"),
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("steelblue", 0.4))
}
```

![](mdgpd_files/figure-gfm/rho-effect-1.png)<!-- -->

``` r
par(op)
```

### Spectral spread ($\lambda$)

The parameter $\lambda$ controls the Poisson rate. Larger $\lambda$
produces larger spectral differences $\Delta_i$, which in turn makes the
$N_i$ more variable relative to $G$:

``` r
set.seed(42)

lambda_vals <- c(0.1, 1, 5)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (lam_v in lambda_vals) {
  Y_lam <- rmdgpd(n_demo, sigma = 2, xi = 0.2, lambda = lam_v, rho = 0.5)
  cor_val <- cor(Y_lam[, 1], Y_lam[, 2])
  plot(jitter(Y_lam[, 1]), jitter(Y_lam[, 2]), pch = 20, cex = 0.3,
       main = bquote(lambda == .(lam_v) ~ " (cor = " * .(round(cor_val, 2)) * ")"),
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("steelblue", 0.4))
}
```

![](mdgpd_files/figure-gfm/lambda-effect-1.png)<!-- -->

``` r
par(op)
```

### Tail heaviness ($\xi$)

The shape parameter $\xi$ controls the tail behaviour. Larger $\xi$
produces heavier tails (occasional very large values):

``` r
set.seed(42)

xi_vals <- c(0, 0.1, 0.5)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (xi_v in xi_vals) {
  Y_xi <- rmdgpd(n_demo, sigma = 2, xi = xi_v, lambda = 1, rho = 0.5)
  plot(jitter(Y_xi[, 1]), jitter(Y_xi[, 2]), pch = 20, cex = 0.3,
       main = bquote(xi == .(xi_v)),
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("steelblue", 0.4))
}
```

![](mdgpd_files/figure-gfm/xi-effect-1.png)<!-- -->

``` r
par(op)
```

## 4. Higher dimensions

The MDGPD generalises naturally to $d \geq 2$ dimensions. Simply pass
the `d` parameter to `rmdgpd()`:

### Trivariate MDGPD ($d = 3$)

``` r
set.seed(42)

Y3 <- rmdgpd(2000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)
cat("Dimensions:", nrow(Y3), "x", ncol(Y3), "\n")
#> Dimensions: 2000 x 3
cat("Column names:", colnames(Y3), "\n")
#> Column names: Y1 Y2 Y3

## Pairwise scatter plots
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
pairs_list <- list(c(1, 2), c(1, 3), c(2, 3))
for (pr in pairs_list) {
  plot(jitter(Y3[, pr[1]]), jitter(Y3[, pr[2]]), pch = 20, cex = 0.3,
       main = paste0("Y", pr[1], " vs Y", pr[2]),
       xlab = paste0("Y", pr[1]), ylab = paste0("Y", pr[2]),
       col = adjustcolor("steelblue", 0.4))
}
```

![](mdgpd_files/figure-gfm/trivariate-1.png)<!-- -->

``` r
par(op)
```

### Pairwise correlation in $d = 3$

``` r
cat("Pairwise correlations:\n")
#> Pairwise correlations:
print(round(cor(Y3), 3))
#>       Y1    Y2    Y3
#> Y1 1.000 0.875 0.871
#> Y2 0.875 1.000 0.853
#> Y3 0.871 0.853 1.000
```

### Dependence at high $\rho$ in $d = 3$

``` r
set.seed(42)
Y3_high <- rmdgpd(2000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.95, d = 3)

op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (pr in pairs_list) {
  plot(jitter(Y3_high[, pr[1]]), jitter(Y3_high[, pr[2]]), pch = 20, cex = 0.3,
       main = paste0("Y", pr[1], " vs Y", pr[2], " (rho=0.95)"),
       xlab = paste0("Y", pr[1]), ylab = paste0("Y", pr[2]),
       col = adjustcolor("steelblue", 0.4))
}
```

![](mdgpd_files/figure-gfm/trivariate-high-rho-1.png)<!-- -->

``` r
par(op)

cat("Pairwise correlations (rho=0.95):\n")
#> Pairwise correlations (rho=0.95):
print(round(cor(Y3_high), 3))
#>       Y1    Y2    Y3
#> Y1 1.000 0.988 0.986
#> Y2 0.988 1.000 0.982
#> Y3 0.986 0.982 1.000
```

### Even higher dimensions

The construction works for any $d \geq 2$:

``` r
set.seed(42)

Y5 <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 5)
cat("5-dimensional MDGPD:", nrow(Y5), "x", ncol(Y5), "\n")
#> 5-dimensional MDGPD: 1000 x 5
cat("Columns:", colnames(Y5), "\n")
#> Columns: Y1 Y2 Y3 Y4 Y5
cat("\nPairwise correlations:\n")
#> 
#> Pairwise correlations:
print(round(cor(Y5), 3))
#>       Y1    Y2    Y3    Y4    Y5
#> Y1 1.000 0.823 0.831 0.836 0.821
#> Y2 0.823 1.000 0.827 0.821 0.831
#> Y3 0.831 0.827 1.000 0.788 0.813
#> Y4 0.836 0.821 0.788 1.000 0.841
#> Y5 0.821 0.831 0.813 0.841 1.000
```

## 5. Zero-inflated MDGPD in higher dimensions

The `rzimdgpd()` function also supports arbitrary $d$:

``` r
set.seed(42)
Y3_zi <- rzimdgpd(2000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5,
                    pi0 = 0.3, d = 3)

joint_zeros_3d <- mean(rowSums(Y3_zi) == 0)
cat("Proportion of all-zero rows (d=3):", round(joint_zeros_3d, 3), "\n")
#> Proportion of all-zero rows (d=3): 0.744
```

## 6. Neural Bayes estimation

### Fitting bivariate MDGPD

The MDGPD has no closed-form likelihood (the construction involves a
discrete maximum and floor operations). Neural Bayes estimation is used
instead, via pre-trained neural networks.

``` r
set.seed(42)

Y <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
fit <- fitegpd(Y, family = "mdgpd", method = "neuralbayes",
               estimator = "npe", nsamples = 2000)
#> Starting Julia ...
summary(fit)
#> Fitting of bivariate MDGPD (Aka-Kratz-Naveau) [Experimental]
#> Method: neuralbayes (npe)  [2000 posterior samples]
#> 
#> Posterior summary:
#>        Median Post.SD   2.5%  97.5%
#> sigma  2.0516  0.2531 1.5609 2.5299
#> xi     0.2234  0.1362 0.0252 0.5066
#> lambda 1.5100  1.2584 0.4087 4.7220
#> rho    0.6767  0.2632 0.0651 0.9167
#> 
#> Note: log-likelihood, AIC, and BIC are not available for neural estimation
#> Number of observations:  1000
```

``` r
plot(fit)
```

![](mdgpd_files/figure-gfm/fit-mdgpd-diagnostics-1.png)<!-- -->

### Fitting ZIMDGPD

``` r
set.seed(42)

Y_zi <- rzimdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5,
                  pi0 = 0.3)
fit_zi <- fitegpd(Y_zi, family = "zimdgpd", method = "neuralbayes",
                   estimator = "npe")
summary(fit_zi)
#> Fitting of zero-inflated bivariate MDGPD (ZIMDGPD) [Experimental]
#> Method: neuralbayes (npe)  [1000 posterior samples]
#> 
#> Posterior summary:
#>        Median Post.SD   2.5%  97.5%
#> sigma  1.7303  0.2993 1.2101 2.3859
#> xi     0.2180  0.1206 0.0135 0.4263
#> lambda 1.6216  1.1322 0.3820 4.6136
#> rho    0.6098  0.2326 0.1104 0.9047
#> pi0    0.2449  0.1049 0.0443 0.4537
#> 
#> Note: log-likelihood, AIC, and BIC are not available for neural estimation
#> Number of observations:  1000
```

### S3 methods

``` r
coef(fit)                # Posterior median estimates
#>     sigma        xi    lambda       rho 
#> 2.0515832 0.2233847 1.5099955 0.6766685
confint(fit)             # 95% credible intervals
#>             2.5 %    97.5 %
#> sigma  1.56091809 2.5299328
#> xi     0.02523016 0.5065719
#> lambda 0.40872400 4.7220143
#> rho    0.06507751 0.9166604
vcov(fit)                # Posterior covariance matrix
#>               sigma            xi      lambda          rho
#> sigma   0.064070258 -0.0238900084 0.002214538 0.0041159956
#> xi     -0.023890008  0.0185562674 0.005013470 0.0008536071
#> lambda  0.002214538  0.0050134705 1.583644524 0.2763197354
#> rho     0.004115996  0.0008536071 0.276319735 0.0692714675
nobs(fit)                # Number of observations
#> [1] 1000
```

### Higher-dimensional estimation

For $d \geq 3$, a separate neural network must be trained for each data
dimension, since the network architecture depends on the input
dimension. Use `train_mdgpd()` with the `data_dim` parameter:

``` r
# Fit 3D data using a pre-trained model
Y3_fit <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)

model_3d_path <- system.file("models", "MDGPD_3D_NPE.bson", package = "egpd")
fit_3d <- fitegpd(Y3_fit, family = "mdgpd", method = "neuralbayes",
                   model.path = model_3d_path, estimator = "npe")
summary(fit_3d)
#> Fitting of 3-variate MDGPD (Aka-Kratz-Naveau) [Experimental]
#> Method: neuralbayes (npe)  [1000 posterior samples]
#> 
#> Posterior summary:
#>        Median Post.SD   2.5%  97.5%
#> sigma  1.9399  0.2863 1.3830 2.4875
#> xi     0.1635  0.1285 0.0145 0.4449
#> lambda 1.3413  1.2467 0.3777 4.8796
#> rho    0.6712  0.2404 0.0995 0.9153
#> 
#> Note: log-likelihood, AIC, and BIC are not available for neural estimation
#> Number of observations:  1000
```

## 7. Comparison: MDGPD vs BDEGPD

The package provides two approaches to bivariate discrete extreme value
modelling. They are suited to different settings:

| Aspect | BDEGPD | MDGPD |
|----|----|----|
| Construction | `floor(continuous BEGPD)` | Poisson generator + geometric max |
| Theoretical basis | Ad hoc discretization | Threshold-stable discrete GPD |
| Parameters | 6 ($\kappa, \sigma, \xi, \theta_L, \theta_U, \theta_\omega$) | 4 ($\sigma, \xi, \lambda, \rho$) |
| Dependence | Asymmetric (separate lower/upper tail) | Exchangeable (equicorrelated) |
| Dimensions | 2 only | $d \geq 2$ |
| Marginals | Discrete EGPD | Discrete GPD |
| Threshold stability | Not guaranteed | Yes (by construction) |

**When to use BDEGPD:** When the data exhibit asymmetric dependence
between lower and upper tails, or when you want to directly compare
continuous and discrete BEGPD fits with the same parameter set.

**When to use MDGPD:** When you want a parsimonious model with
theoretical threshold stability, exchangeable dependence, and the
ability to extend to $d > 2$ dimensions.

``` r
set.seed(42)
n_cmp <- 3000

Y_bdegpd <- rbdegpd(n_cmp, kappa = 2, sigma = 1.5, xi = 0.1,
                     thL = 3, thU = 3, thw = 0.2)
Y_mdgpd  <- rmdgpd(n_cmp, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)

op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(jitter(Y_bdegpd[, 1]), jitter(Y_bdegpd[, 2]), pch = 20, cex = 0.3,
     main = "BDEGPD (6 params)",
     xlab = expression(Y[1]), ylab = expression(Y[2]),
     col = adjustcolor("steelblue", 0.4))
plot(jitter(Y_mdgpd[, 1]), jitter(Y_mdgpd[, 2]), pch = 20, cex = 0.3,
     main = "MDGPD (4 params)",
     xlab = expression(Y[1]), ylab = expression(Y[2]),
     col = adjustcolor("firebrick", 0.4))
```

![](mdgpd_files/figure-gfm/compare-1.png)<!-- -->

``` r
par(op)
```

## 8. Training custom models

### MDGPD training

``` r
# Quick training for 2D MDGPD
paths <- train_mdgpd(
  savepath = tempdir(),
  family = "mdgpd",
  data_dim = 2L,
  estimator = "both",
  quick = TRUE,
  verbose = TRUE
)

# Quick training for 3D MDGPD
paths_3d <- train_mdgpd(
  savepath = tempdir(),
  family = "mdgpd",
  data_dim = 3L,
  estimator = "npe",
  quick = TRUE
)

# ZIMDGPD training
paths_zi <- train_mdgpd(
  savepath = tempdir(),
  family = "zimdgpd",
  data_dim = 2L,
  estimator = "npe",
  quick = TRUE
)
```

### Prior distributions

The training procedure samples parameters from the following priors:

| Parameter         | Prior           | Transform        |
|-------------------|-----------------|------------------|
| $\sigma$          | $U(0.1, 10)$    | $\log$           |
| $\xi$             | $U(0.01, 0.5)$  | $\log$           |
| $\lambda$         | $U(0.01, 5)$    | $\log$           |
| $\rho$            | $U(0.01, 0.99)$ | $\mathrm{logit}$ |
| $\pi_0$ (ZIMDGPD) | $U(0.01, 0.9)$  | $\mathrm{logit}$ |

Model files are named with a dimension-specific prefix: `MDGPD_NPE.bson`
(2D), `MDGPD3D_NPE.bson` (3D), etc.

## References

Aka, S., Kratz, M., and Naveau, P. (2025). Multivariate discrete
generalized Pareto distributions: theory, simulation, and applications
to dry spells. *arXiv preprint* arXiv:2506.19361.
<https://arxiv.org/abs/2506.19361>

Sainsbury-Dale, M., Zammit-Mangion, A., and Huser, R. (2024).
Likelihood-free parameter estimation with neural Bayes estimators. *The
American Statistician*, 78(1), 1–14.
