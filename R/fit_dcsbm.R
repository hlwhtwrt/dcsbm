#' Fit a weighted degree-corrected stochastic block model
#'
#' Fits a Gaussian degree-corrected stochastic block model to a symmetric
#' weighted network via Gibbs sampling, with K (the number of blocks) fixed
#' and specified by the caller. Two optional edge-level covariate effects can
#' be added, each with its own block-pair-specific coefficient: a continuous
#' covariate `X` (its effect is estimated on values centered within each
#' block pair, among positive entries) and/or a binary covariate `X_bin`. An
#' optional "reallocate" Metropolis-Hastings move can be added after the
#' per-node update to improve within-chain mixing; it jointly redistributes
#' node memberships between two existing blocks without changing K.
#'
#' @param W Symmetric numeric weighted adjacency matrix (n x n). The diagonal
#'   is ignored.
#' @param K Integer number of blocks, fixed for this fit. See
#'   \code{\link{select_K}} to compare a grid of K values.
#' @param X Optional symmetric numeric edge-covariate matrix, same dimensions
#'   as `W`. Values <= 0 are treated as "not applicable" for that edge.
#' @param X_bin Optional symmetric binary (0/1) edge-covariate matrix, same
#'   dimensions as `W`, independent of `X`.
#' @param reallocate Logical; add the reallocate MH move after each sweep's
#'   per-node update. Default `FALSE`.
#' @param s2_eta,s2_gamma,s2_beta,s2_delta Prior variances for the degree
#'   effect, block-pair mean, continuous-covariate slope, and binary-covariate
#'   slope respectively.
#' @param a_tau,b_tau Inverse-gamma prior shape/rate for the residual
#'   variance `tau2`.
#' @param alpha_pi Dirichlet prior concentration for block proportions;
#'   defaults to `rep(1, K)`.
#' @param n_iter,burn,thin Total Gibbs iterations, burn-in, and thinning
#'   interval for saved draws.
#' @param seed Optional integer seed, set via `set.seed()` before sampling.
#' @param init_z Optional integer vector of initial block assignments;
#'   defaults to a random assignment covering all K blocks.
#' @param verbose Logical; print progress every 50 iterations.
#'
#' @return An object of class `dcsbm_fit`: a list with the model call,
#'   dimensions, full posterior draws (`$posterior`), and posterior summaries
#'   (`$summary`), including the MAP block assignment `$summary$z_hat`.
#'
#' @seealso \code{\link{select_K}}, \code{\link{print.dcsbm_fit}},
#'   \code{\link{summary.dcsbm_fit}}, \code{\link{plot.dcsbm_fit}}
#'
#' @export
fit_dcsbm <- function(W, K, X = NULL, X_bin = NULL,
                       reallocate = FALSE,
                       s2_eta = 1.0, s2_gamma = 1.0, s2_beta = 1.0, s2_delta = 1.0,
                       a_tau = 2, b_tau = 1, alpha_pi = NULL,
                       n_iter = 2000, burn = 1000, thin = 5,
                       seed = NULL, init_z = NULL, verbose = TRUE) {

  W <- as.matrix(W)
  if (!is.numeric(W)) stop("W must be numeric.")
  if (nrow(W) != ncol(W)) stop("W must be square.")
  if (any(abs(W - t(W)) > 1e-8, na.rm = TRUE)) stop("W must be symmetric.")
  n <- nrow(W)

  node_names <- rownames(W)
  if (is.null(node_names)) node_names <- paste0("node_", seq_len(n))

  K <- as.integer(K)
  if (is.na(K) || K < 2L) stop("K must be an integer >= 2.")
  if (K > n) stop("K cannot exceed the number of nodes.")

  if (!is.null(X)) {
    X <- as.matrix(X)
    if (!all(dim(X) == dim(W))) stop("X must have the same dimensions as W.")
    if (any(abs(X - t(X)) > 1e-8, na.rm = TRUE)) stop("X must be symmetric.")
  }
  if (!is.null(X_bin)) {
    X_bin <- as.matrix(X_bin)
    if (!all(dim(X_bin) == dim(W))) stop("X_bin must have the same dimensions as W.")
    if (any(abs(X_bin - t(X_bin)) > 1e-8, na.rm = TRUE)) stop("X_bin must be symmetric.")
    if (any(!X_bin[upper.tri(X_bin)] %in% c(0, 1)))  stop("X_bin must contain only 0/1 entries.")
  }
  if (!is.null(init_z) && length(init_z) != n) stop("init_z must have length n.")

  if (!is.null(seed)) set.seed(seed)

  raw <- .dcsbm_gibbs(
    W = W, X = X, X_bin = X_bin, K = K,
    s2_eta = s2_eta, s2_gamma = s2_gamma, s2_beta = s2_beta, s2_delta = s2_delta,
    a_tau = a_tau, b_tau = b_tau, alpha_pi = alpha_pi,
    niter = n_iter, burn = burn, thin = thin,
    init_z = init_z, reallocate = reallocate, verbose = verbose
  )

  structure(
    list(
      call = match.call(),
      K = K, n = n, node_names = node_names,
      W = W, X = X, X_bin = X_bin,
      n_iter = n_iter, burn = burn, thin = thin, seed = seed,
      has_x = raw$has_x, has_xbin = raw$has_xbin, reallocate = raw$reallocate,
      posterior = raw$posterior,
      summary = raw$summary
    ),
    class = "dcsbm_fit"
  )
}
