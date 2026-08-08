#' @export
print.dcsbm_fit <- function(x, ...) {
  cat("Weighted degree-corrected SBM fit\n")
  cat("  n =", x$n, "  K =", x$K, "\n")
  cat("  covariates: ",
      if (x$has_x) "continuous " else "",
      if (x$has_xbin) "binary " else "",
      if (!x$has_x && !x$has_xbin) "none",
      "\n", sep = "")
  cat("  reallocate move:", x$reallocate, "\n")
  cat("  n_iter =", x$n_iter, " burn =", x$burn, " thin =", x$thin,
      " (", nrow(x$posterior$z_save), " saved draws)\n", sep = "")
  cat("  mean log-likelihood:", round(x$summary$mean_loglik, 2), "\n")
  cat("  block sizes (MAP):", paste(x$summary$block_sizes_hat, collapse = ", "), "\n")
  if (x$reallocate) {
    cat("  reallocate acceptance rate:", round(x$summary$realloc_accept_rate, 3), "\n")
  }
  invisible(x)
}

#' Summarize a dcsbm_fit object
#'
#' @param object A `dcsbm_fit` object from \code{\link{fit_dcsbm}}.
#' @param ... Unused.
#' @return An object of class `summary.dcsbm_fit`, printed for its side effect.
#' @export
summary.dcsbm_fit <- function(object, ...) {
  structure(
    list(
      call = object$call,
      n = object$n, K = object$K,
      has_x = object$has_x, has_xbin = object$has_xbin, reallocate = object$reallocate,
      n_saved = nrow(object$posterior$z_save),
      block_sizes_hat = object$summary$block_sizes_hat,
      tau2_post_mean = object$summary$tau2_post_mean,
      mean_loglik = object$summary$mean_loglik,
      max_loglik = object$summary$max_loglik,
      realloc_accept_rate = object$summary$realloc_accept_rate,
      gamma_post_mean = object$summary$gamma_post_mean,
      beta_post_mean = object$summary$beta_post_mean,
      delta_post_mean = object$summary$delta_post_mean
    ),
    class = "summary.dcsbm_fit"
  )
}

#' @export
print.summary.dcsbm_fit <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nn =", x$n, " K =", x$K, "\n")
  cat("Block sizes (MAP):", paste(x$block_sizes_hat, collapse = ", "), "\n")
  cat("tau2 posterior mean:", round(x$tau2_post_mean, 4), "\n")
  cat("Mean / max log-likelihood:", round(x$mean_loglik, 2), "/", round(x$max_loglik, 2), "\n")
  if (x$reallocate) cat("Reallocate acceptance rate:", round(x$realloc_accept_rate, 3), "\n")
  cat("\nPosterior mean block-pair intercept (gamma):\n")
  print(round(x$gamma_post_mean, 3))
  if (x$has_x) { cat("\nPosterior mean continuous-covariate slope (beta):\n"); print(round(x$beta_post_mean, 3)) }
  if (x$has_xbin) { cat("\nPosterior mean binary-covariate slope (delta):\n"); print(round(x$delta_post_mean, 3)) }
  invisible(x)
}

#' Plot diagnostics for a dcsbm_fit object
#'
#' @param x A `dcsbm_fit` object.
#' @param type One of `"loglik"` (trace plot of saved log-likelihood draws,
#'   the default), `"block_sizes"` (bar plot of MAP block sizes), or
#'   `"heatmap"` (adjacency matrix reordered by MAP block assignment).
#' @param ... Additional arguments passed to the underlying base plot function.
#' @return Invisibly returns `x`; called for its plotting side effect.
#' @export
plot.dcsbm_fit <- function(x, type = c("loglik", "block_sizes", "heatmap"), ...) {
  type <- match.arg(type)
  if (type == "loglik") {
    graphics::plot(x$posterior$loglik_save, type = "l",
                    xlab = "saved draw", ylab = "log-likelihood",
                    main = "Log-likelihood trace", ...)
  } else if (type == "block_sizes") {
    graphics::barplot(x$summary$block_sizes_hat,
                       names.arg = seq_along(x$summary$block_sizes_hat),
                       xlab = "block", ylab = "size",
                       main = "MAP block sizes", ...)
  } else if (type == "heatmap") {
    ord <- order(x$summary$z_hat)
    graphics::image(x$W[ord, ord], axes = FALSE, main = "Adjacency reordered by MAP blocks", ...)
  }
  invisible(x)
}
