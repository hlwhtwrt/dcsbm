# Base-R smoke tests (no testthat dependency) -- run automatically by
# R CMD check. Uses stopifnot() throughout; a silent exit means all passed.

library(dcsbm)

set.seed(42)
n <- 24
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)
G_true <- matrix(-1.5, K_true, K_true); diag(G_true) <- 2.5

W <- matrix(0, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
  W[i, j] <- rnorm(1, G_true[z_true[i], z_true[j]], 0.25)
}
W <- (W + t(W)) / 2

X <- matrix(0, n, n)
X[upper.tri(X)] <- runif(sum(upper.tri(X)), 0, 1)
X <- X + t(X)
diag(X) <- 0

X_bin <- matrix(rbinom(n * n, 1, 0.3), n, n)
X_bin[lower.tri(X_bin)] <- t(X_bin)[lower.tri(X_bin)]
diag(X_bin) <- 0

## burn/n_iter are generous relative to this tiny n=24 example: with a short
## burn-in, some covariate/reallocate combinations can (deterministically,
## for a fixed seed) hit a run of empty-block iterations and never save a
## draw -- see the K=4 case in the select_K() block below, which is left
## short on purpose to exercise that failure path.
fit_args <- list(n_iter = 400, burn = 200, thin = 2, seed = 1, verbose = FALSE)

## ---- base model (no covariates, no reallocate) ----
fit0 <- do.call(fit_dcsbm, c(list(W = W, K = K_true), fit_args))
stopifnot(
  inherits(fit0, "dcsbm_fit"),
  !fit0$has_x, !fit0$has_xbin, !fit0$reallocate,
  nrow(fit0$posterior$z_save) <= floor((fit_args$n_iter - fit_args$burn) / fit_args$thin),
  nrow(fit0$posterior$z_save) > 0L,
  ncol(fit0$posterior$z_save) == n,
  length(fit0$summary$z_hat) == n,
  all(fit0$summary$block_sizes_hat >= 0),
  sum(fit0$summary$block_sizes_hat) == n,
  is.finite(fit0$summary$mean_loglik)
)

## recovers the planted block structure up to label permutation: each true
## block should map to a single dominant recovered block (contingency table
## has one large entry per row)
tab <- table(z_true, fit0$summary$z_hat)
row_purity <- apply(tab, 1, function(r) max(r) / sum(r))
stopifnot(all(row_purity > 0.7))

## ---- base model + reallocate ----
fit_ra <- do.call(fit_dcsbm, c(list(W = W, K = K_true, reallocate = TRUE), fit_args))
stopifnot(
  fit_ra$reallocate,
  is.numeric(fit_ra$summary$realloc_accept_rate),
  fit_ra$summary$realloc_accept_rate >= 0, fit_ra$summary$realloc_accept_rate <= 1
)

## ---- continuous covariate only ----
fit_x <- do.call(fit_dcsbm, c(list(W = W, K = K_true, X = X), fit_args))
stopifnot(
  fit_x$has_x, !fit_x$has_xbin,
  is.matrix(fit_x$summary$beta_post_mean),
  all(dim(fit_x$summary$beta_post_mean) == K_true),
  is.null(fit_x$summary$delta_post_mean)
)

## ---- binary covariate only ----
fit_xb <- do.call(fit_dcsbm, c(list(W = W, K = K_true, X_bin = X_bin), fit_args))
stopifnot(
  !fit_xb$has_x, fit_xb$has_xbin,
  is.matrix(fit_xb$summary$delta_post_mean),
  is.null(fit_xb$summary$beta_post_mean)
)

## ---- both covariates + reallocate together ----
fit_full <- do.call(fit_dcsbm, c(list(W = W, K = K_true, X = X, X_bin = X_bin, reallocate = TRUE), fit_args))
stopifnot(fit_full$has_x, fit_full$has_xbin, fit_full$reallocate)

## ---- print / summary / plot don't error ----
utils::capture.output(print(fit0))
utils::capture.output(print(summary(fit0)))
grDevices::pdf(NULL)
plot(fit0, type = "loglik")
plot(fit0, type = "block_sizes")
plot(fit0, type = "heatmap")
grDevices::dev.off()

## ---- input validation ----
bad_W <- W; bad_W[1, 2] <- bad_W[1, 2] + 1  # break symmetry
err <- tryCatch({ fit_dcsbm(bad_W, K = K_true, n_iter = 10, burn = 5, verbose = FALSE); "no error" },
                 error = function(e) "error")
stopifnot(err == "error")

err2 <- tryCatch({ fit_dcsbm(W, K = n + 1, n_iter = 10, burn = 5, verbose = FALSE); "no error" },
                  error = function(e) "error")
stopifnot(err2 == "error")

## ---- select_K ----
## K=4 with this short a burn-in deterministically never saves a draw for
## this seed/data (exercises the "one K in the grid fails" graceful-degradation
## path); K=2 and K=3 should still succeed and the sweep should not abort.
warned <- FALSE
grid <- withCallingHandlers(
  select_K(W, K_grid = 2:4, n_iter = 80, burn = 30, thin = 2, verbose = FALSE),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }
)
stopifnot(
  inherits(grid, "dcsbm_Kgrid"),
  nrow(grid) == 3L,
  all(c("K", "seed", "mean_loglik", "max_loglik", "n_occupied_blocks", "min_block_size") %in% names(grid)),
  length(attr(grid, "fits")) == 3L,
  warned,
  is.na(grid$mean_loglik[grid$K == 4]),
  is.null(attr(grid, "fits")[[which(grid$K == 4)]]),
  !anyNA(grid$mean_loglik[grid$K != 4])
)
grDevices::pdf(NULL)
plot(grid)
grDevices::dev.off()

cat("All dcsbm smoke tests passed.\n")
