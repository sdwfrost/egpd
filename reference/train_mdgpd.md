# \\Experimental\\ Train a neural estimator for MDGPD

Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
(NBE) networks for MDGPD inference. Based on the Aka, Kratz & Naveau
(2025) framework. Requires Julia (\>= 1.11) with NeuralEstimators.jl and
Flux.jl installed, plus the R packages JuliaConnectoR and
NeuralEstimators.

## Usage

``` r
train_mdgpd(
  savepath = "inst/models",
  family = c("mdgpd", "zimdgpd"),
  data_dim = 2L,
  estimator = c("both", "npe", "nbe"),
  K = NULL,
  m = NULL,
  epochs = NULL,
  stopping_epochs = 10,
  quick = FALSE,
  mc.cores = parallel::detectCores() - 1L,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- savepath:

  character; directory to save trained model files (.bson). Created if
  it does not exist.

- family:

  character; `"mdgpd"` (4 parameters) or `"zimdgpd"` (5 parameters, adds
  pi0).

- data_dim:

  integer \>= 2; dimension of the multivariate data (default `2L`). A
  separate neural network is trained for each data dimension.

- estimator:

  character; which estimator(s) to train: `"both"` (default), `"npe"`,
  or `"nbe"`.

- K:

  integer; number of training parameter sets to sample from the prior.

- m:

  integer vector; range of sample sizes for simulated datasets.

- epochs:

  integer; maximum number of training epochs.

- stopping_epochs:

  integer; early stopping patience.

- quick:

  logical; if `TRUE`, uses reduced training settings.

- mc.cores:

  integer; number of parallel cores for data simulation.

- seed:

  integer or `NULL`; random seed for reproducibility.

- verbose:

  logical; if `TRUE`, print progress messages.

## Value

A named list of file paths to saved .bson model files (invisible).

## Details

This function is **experimental** and its interface may change.

The training procedure:

1.  Samples parameters from uniform priors: `sigma ~ U(0.1, 10)`,
    `xi ~ U(0.01, 0.5)`, `lambda ~ U(0.01, 5)`, `rho ~ U(0.01, 0.99)`,
    and for ZIMDGPD: `pi0 ~ U(0.01, 0.9)`.

2.  Simulates `data_dim`-dimensional MDGPD data using
    [`rmdgpd`](https://sdwfrost.github.io/egpd/reference/rmdgpd.md) or
    [`rzimdgpd`](https://sdwfrost.github.io/egpd/reference/rzimdgpd.md).

3.  Applies variance-stabilizing (signed log) and parameter transforms
    (`log` for sigma/xi/lambda, `logit` for rho and pi0).

4.  Trains the neural network(s).

5.  Saves trained states as .bson files.

Model files are named with a dimension-specific prefix (e.g.,
`MDGPD_NPE.bson` for `data_dim = 2`, `MDGPD_3D_NPE.bson` for
`data_dim = 3`).

## References

Aka, S., Kratz, M., and Naveau, P. (2025). Multivariate discrete
generalized Pareto distributions: theory, simulation, and applications
to dry spells. *arXiv preprint* arXiv:2506.19361.

## See also

[`rmdgpd`](https://sdwfrost.github.io/egpd/reference/rmdgpd.md),
[`rzimdgpd`](https://sdwfrost.github.io/egpd/reference/rzimdgpd.md),
[`train_begpd`](https://sdwfrost.github.io/egpd/reference/train_begpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 2D (default)
paths <- train_mdgpd(savepath = tempdir(), family = "mdgpd", quick = TRUE)
# 3D
paths3 <- train_mdgpd(savepath = tempdir(), family = "mdgpd",
                       data_dim = 3L, quick = TRUE)
} # }
```
