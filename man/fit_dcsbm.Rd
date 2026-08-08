\name{fit_dcsbm}
\alias{fit_dcsbm}
\title{Fit a weighted degree-corrected stochastic block model}
\description{
Fits a Gaussian degree-corrected stochastic block model to a symmetric
weighted network via Gibbs sampling, with K (the number of blocks) fixed
and specified by the caller. Two optional edge-level covariate effects can
be added, each with its own block-pair-specific coefficient: a continuous
covariate \code{X} (its effect is estimated on values centered within each
block pair, among positive entries) and/or a binary covariate \code{X_bin}.
An optional "reallocate" Metropolis-Hastings move can be added after the
per-node update to improve within-chain mixing; it jointly redistributes
node memberships between two existing blocks without changing K.
}
\usage{
fit_dcsbm(W, K, X = NULL, X_bin = NULL,
          reallocate = FALSE,
          s2_eta = 1.0, s2_gamma = 1.0, s2_beta = 1.0, s2_delta = 1.0,
          a_tau = 2, b_tau = 1, alpha_pi = NULL,
          n_iter = 2000, burn = 1000, thin = 5,
          seed = NULL, init_z = NULL, verbose = TRUE)
}
\arguments{
  \item{W}{Symmetric numeric weighted adjacency matrix (n x n). The diagonal is ignored.}
  \item{K}{Integer number of blocks, fixed for this fit. See \code{\link{select_K}} to compare a grid of K values.}
  \item{X}{Optional symmetric numeric edge-covariate matrix, same dimensions as \code{W}. Values <= 0 are treated as "not applicable" for that edge.}
  \item{X_bin}{Optional symmetric binary (0/1) edge-covariate matrix, same dimensions as \code{W}, independent of \code{X}.}
  \item{reallocate}{Logical; add the reallocate MH move after each sweep's per-node update. Default \code{FALSE}.}
  \item{s2_eta}{Prior variance for the degree effect.}
  \item{s2_gamma}{Prior variance for the block-pair mean.}
  \item{s2_beta}{Prior variance for the continuous-covariate slope.}
  \item{s2_delta}{Prior variance for the binary-covariate slope.}
  \item{a_tau}{Inverse-gamma prior shape for the residual variance \code{tau2}.}
  \item{b_tau}{Inverse-gamma prior rate for the residual variance \code{tau2}.}
  \item{alpha_pi}{Dirichlet prior concentration for block proportions; defaults to \code{rep(1, K)}.}
  \item{n_iter}{Total Gibbs iterations.}
  \item{burn}{Burn-in iterations discarded before saving draws.}
  \item{thin}{Thinning interval for saved draws.}
  \item{seed}{Optional integer seed, set via \code{set.seed()} before sampling.}
  \item{init_z}{Optional integer vector of initial block assignments; defaults to a random assignment covering all K blocks.}
  \item{verbose}{Logical; print progress every 50 iterations.}
}
\value{
An object of class \code{dcsbm_fit}: a list with the model call, dimensions,
full posterior draws (\code{$posterior}), and posterior summaries
(\code{$summary}), including the MAP block assignment \code{$summary$z_hat}.
}
\seealso{
\code{\link{select_K}}
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

fit <- fit_dcsbm(W, K = 3, n_iter = 200, burn = 100, thin = 2, seed = 1, verbose = FALSE)
print(fit)
}
}
