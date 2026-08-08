#' Compare a grid of K values
#'
#' Repeatedly calls \code{\link{fit_dcsbm}} across a grid of `K` (and,
#' optionally, multiple seeds per K) and collects a criterion table to help
#' choose K -- e.g. by an elbow in mean log-likelihood, or by checking how
#' often blocks come out empty/tiny at a given K.
#'
#' @param W,X,X_bin,reallocate,n_iter,burn,thin,... Passed through to each
#'   \code{\link{fit_dcsbm}} call.
#' @param K_grid Integer vector of K values to try.
#' @param seeds Integer vector of seeds to try per K value; default `1L`.
#' @param verbose Logical; print a one-line progress message before each fit
#'   (each individual fit itself runs with `verbose = FALSE`).
#'
#' @return A `data.frame` (class `dcsbm_Kgrid`) with one row per (K, seed):
#'   `K`, `seed`, `mean_loglik`, `max_loglik`, `n_occupied_blocks`,
#'   `min_block_size`, and `realloc_accept_rate` (`NA` unless
#'   `reallocate = TRUE`). If a (K, seed) fit fails outright -- e.g. no
#'   post-burn draw ever had all K blocks occupied -- its row is filled with
#'   `NA` and a warning is raised, but the sweep continues over the rest of
#'   the grid. The individual fitted `dcsbm_fit` objects are kept in
#'   `attr(., "fits")` (`NULL` for any failed fit). Has a
#'   \code{\link{plot.dcsbm_Kgrid}} method.
#'
#' @export
select_K <- function(W, K_grid, seeds = 1L, X = NULL, X_bin = NULL,
                      reallocate = FALSE, n_iter = 2000, burn = 1000, thin = 5,
                      verbose = TRUE, ...) {
  K_grid <- sort(unique(as.integer(K_grid)))
  seeds  <- as.integer(seeds)
  grid   <- expand.grid(K = K_grid, seed = seeds, KEEP.OUT.ATTRS = FALSE)

  fits <- vector("list", nrow(grid))
  rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    K <- grid$K[i]; seed <- grid$seed[i]
    if (verbose) message(sprintf("select_K: fitting K=%d seed=%d (%d/%d)", K, seed, i, nrow(grid)))

    fit <- tryCatch(
      fit_dcsbm(W = W, K = K, X = X, X_bin = X_bin, reallocate = reallocate,
                n_iter = n_iter, burn = burn, thin = thin, seed = seed,
                verbose = FALSE, ...),
      error = function(e) { warning("K=", K, " seed=", seed, " failed: ", conditionMessage(e)); NULL }
    )
    fits[i] <- list(fit)  # NB: fits[[i]] <- fit would delete the slot entirely when fit is NULL
    rows[[i]] <- if (is.null(fit)) {
      data.frame(K = K, seed = seed, mean_loglik = NA_real_, max_loglik = NA_real_,
                 n_occupied_blocks = NA_integer_, min_block_size = NA_integer_,
                 realloc_accept_rate = NA_real_)
    } else {
      data.frame(
        K = K, seed = seed,
        mean_loglik = fit$summary$mean_loglik,
        max_loglik  = fit$summary$max_loglik,
        n_occupied_blocks = sum(fit$summary$block_sizes_hat > 0),
        min_block_size    = min(fit$summary$block_sizes_hat),
        realloc_accept_rate = if (reallocate) fit$summary$realloc_accept_rate else NA_real_
      )
    }
  }

  out <- do.call(rbind, rows)
  attr(out, "fits") <- fits
  class(out) <- c("dcsbm_Kgrid", class(out))
  out
}

#' Plot a K-selection criterion curve
#'
#' @param x A `dcsbm_Kgrid` object from \code{\link{select_K}}.
#' @param ... Additional arguments passed to the underlying `plot()` call.
#' @return Invisibly returns `x`; called for its plotting side effect.
#' @export
plot.dcsbm_Kgrid <- function(x, ...) {
  agg <- stats::aggregate(mean_loglik ~ K, data = x, FUN = mean)
  graphics::plot(agg$K, agg$mean_loglik, type = "b",
                  xlab = "K", ylab = "mean log-likelihood (averaged over seeds)",
                  main = "K selection", ...)
  if (length(unique(x$seed)) > 1L) {
    graphics::points(x$K, x$mean_loglik, pch = 1, col = "grey50")
  }
  invisible(x)
}
