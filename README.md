# dcsbm

An R package for fitting Bayesian weighted degree-corrected stochastic block
models (DC-SBMs) via Gibbs sampling. Developed as part of my PhD dissertation
research on community detection in weighted networks, and generalized here
into a standalone tool.

## Features

- Gaussian DC-SBM likelihood for weighted (not just binary) networks, with a
  per-node degree-correction term.
- A fixed number of blocks `K`, chosen by the user (see `select_K()` below
  for comparing a grid of K values).
- Two optional edge-level covariate effects, each with its own block-pair-
  specific coefficient: a continuous covariate and/or a binary indicator.
- An optional "reallocate" Metropolis-Hastings move that jointly
  redistributes node memberships between two existing blocks in a single
  step, to improve within-chain mixing beyond what single-node Gibbs updates
  achieve alone.

## Installation

```r
# install.packages("remotes")
remotes::install_github("<github-username>/dcsbm")
```

## Usage

```r
library(dcsbm)

# toy 3-block network
set.seed(1)
n <- 30
z_true <- rep(1:3, each = 10)
G_true <- matrix(c(2, -1, -1, -1, 2, -1, -1, -1, 2), 3, 3)
W <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
  W[i, j] <- rnorm(1, G_true[z_true[i], z_true[j]], 0.3)
}
W <- (W + t(W)) / 2

fit <- fit_dcsbm(W, K = 3, n_iter = 2000, burn = 1000, thin = 5, seed = 1)
print(fit)
summary(fit)
plot(fit, type = "heatmap")

# compare a grid of K values
grid <- select_K(W, K_grid = 2:6, n_iter = 2000, burn = 1000, thin = 5)
plot(grid)
```

Add an edge covariate and the reallocate move:

```r
fit <- fit_dcsbm(W, K = 3, X = my_edge_covariate, reallocate = TRUE)
```

## Background

The reallocate move addresses a specific mixing problem: single-node Gibbs
updates only move a node if that node's individual reassignment looks
favorable in isolation, so a chain can get stuck when several nodes would
only want to swap blocks together. In one real-data comparison, adding this
move increased the effective sample size of the log-likelihood trace by
roughly 24x on an otherwise identical chain -- see `?fit_dcsbm` for details
on how it's implemented (a Dahl (2003) sequentially-allocated
restricted-Gibbs launch, scored as one Metropolis-Hastings step).

## License

MIT
