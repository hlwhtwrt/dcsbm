\name{select_K}
\alias{select_K}
\alias{plot.dcsbm_Kgrid}
\title{Compare a grid of K values}
\description{
Repeatedly calls \code{\link{fit_dcsbm}} across a grid of \code{K} (and,
optionally, multiple seeds per K) and collects a criterion table to help
choose K -- e.g. by an elbow in mean log-likelihood, or by checking how
often blocks come out empty/tiny at a given K.
}
\usage{
select_K(W, K_grid, seeds = 1L, X = NULL, X_bin = NULL,
         reallocate = FALSE, n_iter = 2000, burn = 1000, thin = 5,
         verbose = TRUE, ...)

\method{plot}{dcsbm_Kgrid}(x, ...)
}
\arguments{
  \item{W}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{K_grid}{Integer vector of K values to try.}
  \item{seeds}{Integer vector of seeds to try per K value; default \code{1L}.}
  \item{X}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{X_bin}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{reallocate}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{n_iter}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{burn}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{thin}{Passed through to each \code{\link{fit_dcsbm}} call.}
  \item{verbose}{Logical; print a one-line progress message before each fit (each individual fit itself runs with \code{verbose = FALSE}).}
  \item{x}{A \code{dcsbm_Kgrid} object from \code{select_K}, for the \code{plot} method.}
  \item{...}{Additional arguments passed to \code{\link{fit_dcsbm}} (for \code{select_K}) or to the underlying \code{plot()} call (for the \code{plot} method).}
}
\value{
\code{select_K} returns a \code{data.frame} (class \code{dcsbm_Kgrid}) with
one row per (K, seed): \code{K}, \code{seed}, \code{mean_loglik},
\code{max_loglik}, \code{n_occupied_blocks}, \code{min_block_size}, and
\code{realloc_accept_rate} (\code{NA} unless \code{reallocate = TRUE}). If a
(K, seed) fit fails outright -- e.g. no post-burn draw ever had all K blocks
occupied -- its row is filled with \code{NA} and a warning is raised, but the
sweep continues over the rest of the grid. The individual fitted
\code{dcsbm_fit} objects are kept in \code{attr(., "fits")} (\code{NULL} for
any failed fit).

The \code{plot} method invisibly returns \code{x}; called for its plotting
side effect.
}
\examples{
\donttest{
set.seed(1)
n <- 30
z_true <- rep(1:3, each = 10)
G_true <- matrix(c(2, -1, -1, -1, 2, -1, -1, -1, 2), 3, 3)
W <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
  W[i, j] <- rnorm(1, G_true[z_true[i], z_true[j]], 0.3)
}
W <- (W + t(W)) / 2

grid <- select_K(W, K_grid = 2:4, n_iter = 100, burn = 50, thin = 2, verbose = FALSE)
plot(grid)
}
}
