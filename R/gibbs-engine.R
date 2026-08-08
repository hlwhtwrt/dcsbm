# Internal Gibbs sampler shared by every fit_dcsbm() call.
#
# Handles four combinations via has_x/has_xbin flags rather than four
# duplicated loops: no covariates, continuous edge covariate only (X),
# binary edge covariate only (X_bin), or both. The `reallocate` flag adds an
# extra Metropolis-Hastings move after the per-node z update that jointly
# redistributes members between two randomly chosen existing blocks (K is
# never changed) via a Dahl (2003) sequentially-allocated restricted-Gibbs
# launch. This mirrors validated exploratory samplers where single-node
# Gibbs updates alone mixed slowly within a mode (~14.5/400 effective sample
# size on log-likelihood without it, ~346/400 with it, on one real network).

.make_sym_mat <- function(vec, K) {
  M <- matrix(0, K, K)
  idx <- 1L
  for (k in seq_len(K)) for (l in k:K) { M[k, l] <- M[l, k] <- vec[idx]; idx <- idx + 1L }
  M
}

.sym_mat_to_vec <- function(M, K) {
  out <- numeric(K * (K + 1L) / 2L)
  idx <- 1L
  for (k in seq_len(K)) for (l in k:K) { out[idx] <- M[k, l]; idx <- idx + 1L }
  out
}

.dcsbm_gibbs <- function(W, X = NULL, X_bin = NULL, K,
                          s2_eta = 1.0, s2_gamma = 1.0, s2_beta = 1.0, s2_delta = 1.0,
                          a_tau = 2, b_tau = 1, alpha_pi = NULL,
                          niter = 2000, burn = 1000, thin = 5,
                          init_z = NULL, reallocate = FALSE, verbose = TRUE) {

  n <- nrow(W)
  diag(W) <- 0
  has_x    <- !is.null(X)
  has_xbin <- !is.null(X_bin)
  if (has_x)    diag(X) <- 0
  if (has_xbin) diag(X_bin) <- 0

  if (is.null(alpha_pi)) alpha_pi <- rep(1, K)

  ut <- which(upper.tri(W), arr.ind = TRUE)
  ii <- ut[, 1L]; jj <- ut[, 2L]
  m_edges <- length(ii)
  w_vec  <- W[cbind(ii, jj)]
  x_vec  <- if (has_x)    X[cbind(ii, jj)]     else NULL
  xb_vec <- if (has_xbin) X_bin[cbind(ii, jj)] else NULL

  P <- K * (K + 1L) / 2L

  z <- if (!is.null(init_z)) as.integer(init_z) else sample(seq_len(K), n, replace = TRUE)
  while (length(unique(z)) < K) z <- sample(seq_len(K), n, replace = TRUE)

  eta   <- rep(0, n)
  G     <- matrix(0, K, K)
  B_mat <- if (has_x)    matrix(0, K, K) else NULL
  D_mat <- if (has_xbin) matrix(0, K, K) else NULL
  tau2  <- stats::var(w_vec)
  pi_k  <- rep(1 / K, K)

  n_save <- floor((niter - burn) / thin)
  if (n_save <= 0L) stop("niter must exceed burn by at least `thin`.")

  z_save      <- matrix(NA_integer_, n_save, n)
  eta_save    <- matrix(NA_real_,    n_save, n)
  gamma_save  <- matrix(NA_real_,    n_save, P)
  beta_save   <- if (has_x)    matrix(NA_real_, n_save, P) else NULL
  delta_save  <- if (has_xbin) matrix(NA_real_, n_save, P) else NULL
  tau2_save   <- rep(NA_real_, n_save)
  loglik_save <- rep(NA_real_, n_save)
  realloc_accept <- if (reallocate) logical(niter) else NULL
  save_idx <- 0L

  # ---- covariate-centering helper (only meaningful when has_x) ----
  compute_xbar_cis <- function(z) {
    xbar <- matrix(0, K, K)
    for (k in seq_len(K)) {
      for (l in k:K) {
        ik <- which(z == k); il <- which(z == l)
        if (k == l) {
          if (length(ik) < 2L) next
          x_kl <- X[ik, ik][upper.tri(X[ik, ik, drop = FALSE])]
        } else {
          if (length(ik) == 0L || length(il) == 0L) next
          x_kl <- as.vector(X[ik, il, drop = FALSE])
        }
        pos <- x_kl[x_kl > 0]
        v <- if (length(pos) > 0L) mean(pos) else 0
        xbar[k, l] <- xbar[l, k] <- v
      }
    }
    xbar
  }

  # ---- mean function for one block pair, given current parameters ----
  # e: eta_i + eta_j vector; xt: centered continuous covariate (or NULL); xb: binary covariate (or NULL)
  .mu_pair <- function(k, l, e, xt, xb) {
    mu <- G[k, l] + e
    if (has_x)    mu <- mu + B_mat[k, l] * xt
    if (has_xbin) mu <- mu + D_mat[k, l] * xb
    mu
  }

  # NB: `xt_vec` is the FROZEN (start-of-iteration, pre-z-update) centered
  # covariate vector, deliberately not recomputed against the post-update z
  # here -- this reproduces the validated original sampler's behavior exactly
  # (which combines frozen covariate-centering with the freshly updated z for
  # the G/beta/delta block-pair lookups).
  compute_loglik <- function(xt_vec) {
    e <- eta[ii] + eta[jj]
    mu <- .mu_pair_vec(z[ii], z[jj], e, xt_vec, xb_vec)
    sum(stats::dnorm(w_vec, mu, sqrt(tau2), log = TRUE))
  }

  # vectorized version of .mu_pair for the full edge set (used by compute_loglik)
  .mu_pair_vec <- function(zi, zj, e, xt, xb) {
    mu <- G[cbind(zi, zj)] + e
    if (has_x)    mu <- mu + B_mat[cbind(zi, zj)] * xt
    if (has_xbin) mu <- mu + D_mat[cbind(zi, zj)] * xb
    mu
  }

  # ==== reallocate-move helpers (only built/used when reallocate == TRUE) ====

  block_pair_ll_pointval <- function(k, l, z, xbar_cis) {
    ik <- which(z == k); il <- which(z == l)
    if (k == l) {
      if (length(ik) < 2L) return(0)
      subW <- W[ik, ik, drop = FALSE]; subE <- outer(eta[ik], eta[ik], "+")
      utk <- upper.tri(subW)
      w <- subW[utk]; e <- subE[utk]
      xt <- if (has_x)    { xr <- X[ik, ik, drop = FALSE][utk]; ifelse(xr > 0, xr - xbar_cis[k, l], 0) } else NULL
      xb <- if (has_xbin) X_bin[ik, ik, drop = FALSE][utk] else NULL
    } else {
      if (length(ik) == 0L || length(il) == 0L) return(0)
      w <- as.vector(W[ik, il, drop = FALSE]); e <- as.vector(outer(eta[ik], eta[il], "+"))
      xt <- if (has_x)    { xr <- as.vector(X[ik, il, drop = FALSE]); ifelse(xr > 0, xr - xbar_cis[k, l], 0) } else NULL
      xb <- if (has_xbin) as.vector(X_bin[ik, il, drop = FALSE]) else NULL
    }
    mu <- .mu_pair(k, l, e, xt, xb)
    sum(stats::dnorm(w, mu, sqrt(tau2), log = TRUE))
  }

  affected_pairs_ll <- function(focal_labels, z, xbar_cis) {
    total <- 0
    for (fi in seq_along(focal_labels)) {
      k <- focal_labels[fi]
      total <- total + block_pair_ll_pointval(k, k, z, xbar_cis)
      for (l in seq_len(K)) {
        if (l == k) next
        if (fi > 1L && l %in% focal_labels[seq_len(fi - 1L)]) next
        total <- total + block_pair_ll_pointval(k, l, z, xbar_cis)
      }
    }
    total
  }

  restricted_gibbs_launch <- function(anchor_i, anchor_j, free_members, c1, c2, xbar_cis, force_target = NULL) {
    LAB1 <- -1L; LAB2 <- -2L
    z_temp <- rep(0L, n)
    z_temp[anchor_i] <- LAB1
    z_temp[anchor_j] <- LAB2

    s_pair_ll <- function(zt) {
      ik1 <- which(zt == LAB1); ik2 <- which(zt == LAB2)
      ll <- 0
      if (length(ik1) >= 2L) {
        subW <- W[ik1, ik1, drop = FALSE]; subE <- outer(eta[ik1], eta[ik1], "+")
        utk <- upper.tri(subW); w <- subW[utk]; e <- subE[utk]
        xt <- if (has_x)    { xr <- X[ik1, ik1, drop = FALSE][utk]; ifelse(xr > 0, xr - xbar_cis[c1, c1], 0) } else NULL
        xb <- if (has_xbin) X_bin[ik1, ik1, drop = FALSE][utk] else NULL
        ll <- ll + sum(stats::dnorm(w, .mu_pair(c1, c1, e, xt, xb), sqrt(tau2), log = TRUE))
      }
      if (length(ik2) >= 2L) {
        subW <- W[ik2, ik2, drop = FALSE]; subE <- outer(eta[ik2], eta[ik2], "+")
        utk <- upper.tri(subW); w <- subW[utk]; e <- subE[utk]
        xt <- if (has_x)    { xr <- X[ik2, ik2, drop = FALSE][utk]; ifelse(xr > 0, xr - xbar_cis[c2, c2], 0) } else NULL
        xb <- if (has_xbin) X_bin[ik2, ik2, drop = FALSE][utk] else NULL
        ll <- ll + sum(stats::dnorm(w, .mu_pair(c2, c2, e, xt, xb), sqrt(tau2), log = TRUE))
      }
      if (length(ik1) >= 1L && length(ik2) >= 1L) {
        w <- as.vector(W[ik1, ik2, drop = FALSE]); e <- as.vector(outer(eta[ik1], eta[ik2], "+"))
        xt <- if (has_x)    { xr <- as.vector(X[ik1, ik2, drop = FALSE]); ifelse(xr > 0, xr - xbar_cis[c1, c2], 0) } else NULL
        xb <- if (has_xbin) as.vector(X_bin[ik1, ik2, drop = FALSE]) else NULL
        ll <- ll + sum(stats::dnorm(w, .mu_pair(c1, c2, e, xt, xb), sqrt(tau2), log = TRUE))
      }
      ll
    }

    ord <- sample(free_members)
    log_q <- 0
    for (m in ord) {
      z_try1 <- z_temp; z_try1[m] <- LAB1
      z_try2 <- z_temp; z_try2[m] <- LAB2
      logp1 <- log(pi_k[c1] + 1e-300) + s_pair_ll(z_try1)
      logp2 <- log(pi_k[c2] + 1e-300) + s_pair_ll(z_try2)
      mx <- max(logp1, logp2)
      p1 <- exp(logp1 - mx); p2 <- exp(logp2 - mx); p1 <- p1 / (p1 + p2)

      if (!is.null(force_target)) {
        chosen <- force_target[[as.character(m)]]
        pchoice <- if (chosen == LAB1) p1 else (1 - p1)
        log_q <- log_q + log(pchoice)
      } else {
        chosen <- if (stats::runif(1) < p1) LAB1 else LAB2
        log_q <- log_q + log(if (chosen == LAB1) p1 else (1 - p1))
      }
      z_temp[m] <- chosen
    }
    list(z_temp = z_temp, log_q = log_q)
  }

  propose_reallocate <- function(z, xbar_cis) {
    ij <- sample.int(n, 2L)
    i <- ij[1]; j <- ij[2]
    c1 <- z[i]; c2 <- z[j]
    if (c1 == c2) return(list(z = z, accepted = FALSE))

    members1 <- which(z == c1); members2 <- which(z == c2)
    free <- setdiff(c(members1, members2), c(i, j))

    launch_fwd <- restricted_gibbs_launch(i, j, free, c1, c2, xbar_cis)
    z_temp <- launch_fwd$z_temp
    n1_new <- sum(z_temp == -1L); n2_new <- sum(z_temp == -2L)
    if (n1_new == 0L || n2_new == 0L) return(list(z = z, accepted = FALSE))

    z_new <- z
    if (length(free) > 0L) {
      z_new[free[z_temp[free] == -1L]] <- c1
      z_new[free[z_temp[free] == -2L]] <- c2
    }

    force_target <- as.list(ifelse(z[free] == c1, -1L, -2L))
    names(force_target) <- as.character(free)
    launch_rev <- restricted_gibbs_launch(i, j, free, c1, c2, xbar_cis, force_target = force_target)

    ll_before <- affected_pairs_ll(c(c1, c2), z,     xbar_cis)
    ll_after  <- affected_pairs_ll(c(c1, c2), z_new, xbar_cis)

    reassigned <- free[z_new[free] != z[free]]
    log_prior_ratio <- if (length(reassigned) > 0L) {
      sum(log(pi_k[z_new[reassigned]] + 1e-300) - log(pi_k[z[reassigned]] + 1e-300))
    } else 0

    log_alpha <- log_prior_ratio + (ll_after - ll_before) + launch_rev$log_q - launch_fwd$log_q

    if (log(stats::runif(1)) < log_alpha) list(z = z_new, accepted = TRUE) else list(z = z, accepted = FALSE)
  }

  # ==== main Gibbs loop ====
  for (iter in seq_len(niter)) {

    xbar_cis <- if (has_x) compute_xbar_cis(z) else NULL
    xt_vec   <- if (has_x) ifelse(x_vec > 0, x_vec - xbar_cis[cbind(z[ii], z[jj])], 0) else NULL

    ## 1. block-pair parameters: gamma always, beta if has_x, delta if has_xbin
    for (k in seq_len(K)) {
      for (l in k:K) {
        ik <- which(z == k); il <- which(z == l)

        if (k == l) {
          if (length(ik) < 2L) next
          subW <- W[ik, ik, drop = FALSE]; subE <- outer(eta[ik], eta[ik], "+")
          utk <- upper.tri(subW)
          w_kl <- subW[utk]; e_kl <- subE[utk]
          x_kl  <- if (has_x)    X[ik, ik, drop = FALSE][utk]     else NULL
          xb_kl <- if (has_xbin) X_bin[ik, ik, drop = FALSE][utk] else NULL
        } else {
          if (length(ik) == 0L || length(il) == 0L) next
          w_kl <- as.vector(W[ik, il, drop = FALSE]); e_kl <- as.vector(outer(eta[ik], eta[il], "+"))
          x_kl  <- if (has_x)    as.vector(X[ik, il, drop = FALSE])     else NULL
          xb_kl <- if (has_xbin) as.vector(X_bin[ik, il, drop = FALSE]) else NULL
        }
        n_kl <- length(w_kl)

        xt_kl <- if (has_x) { xbar_kl <- if (any(x_kl > 0)) mean(x_kl[x_kl > 0]) else 0; ifelse(x_kl > 0, x_kl - xbar_kl, 0) } else NULL

        old_bd <- 0
        if (has_x)    old_bd <- old_bd + B_mat[k, l] * xt_kl
        if (has_xbin) old_bd <- old_bd + D_mat[k, l] * xb_kl
        r_gam  <- w_kl - old_bd - e_kl
        pp_gam <- n_kl / tau2 + 1 / s2_gamma
        g_new  <- stats::rnorm(1L, sum(r_gam) / tau2 / pp_gam, sqrt(1 / pp_gam))
        G[k, l] <- G[l, k] <- g_new

        if (has_x) {
          d_term <- if (has_xbin) D_mat[k, l] * xb_kl else 0
          r_bet  <- w_kl - g_new - d_term - e_kl
          pp_bet <- sum(xt_kl^2) / tau2 + 1 / s2_beta
          b_new  <- stats::rnorm(1L, sum(r_bet * xt_kl) / tau2 / pp_bet, sqrt(1 / pp_bet))
          B_mat[k, l] <- B_mat[l, k] <- b_new
        }
        if (has_xbin) {
          b_term <- if (has_x) B_mat[k, l] * xt_kl else 0
          r_del  <- w_kl - g_new - b_term - e_kl
          pp_del <- sum(xb_kl^2) / tau2 + 1 / s2_delta
          d_new  <- stats::rnorm(1L, sum(r_del * xb_kl) / tau2 / pp_del, sqrt(1 / pp_del))
          D_mat[k, l] <- D_mat[l, k] <- d_new
        }
      }
    }

    ## 2. eta_i
    for (i in seq_len(n)) {
      zi_i <- z[i]; j_idx <- seq_len(n)[-i]; zj <- z[j_idx]
      e_no_eta_i <- eta[j_idx]
      xt_ij <- if (has_x) { xr <- X[i, j_idx]; ifelse(xr > 0, xr - xbar_cis[cbind(zi_i, zj)], 0) } else NULL
      xb_ij <- if (has_xbin) X_bin[i, j_idx] else NULL
      mu_no_eta_i <- .mu_pair_vec(rep(zi_i, length(zj)), zj, e_no_eta_i, xt_ij, xb_ij)
      resid_i <- W[i, j_idx] - mu_no_eta_i
      pp_eta  <- (n - 1L) / tau2 + 1 / s2_eta
      eta[i]  <- stats::rnorm(1L, sum(resid_i) / tau2 / pp_eta, sqrt(1 / pp_eta))
    }
    eta <- eta - mean(eta)

    ## 3. tau2
    e_all <- eta[ii] + eta[jj]
    mu_all <- .mu_pair_vec(z[ii], z[jj], e_all, xt_vec, xb_vec)
    rss <- sum((w_vec - mu_all)^2)
    tau2 <- 1 / stats::rgamma(1L, shape = a_tau + m_edges / 2, rate = b_tau + rss / 2)

    ## 4. pi
    pi_k <- as.numeric(gtools::rdirichlet(1L, alpha_pi + tabulate(z, nbins = K)))

    ## 5. z_i (frozen xbar_cis from top of iteration)
    for (i in seq_len(n)) {
      j_idx <- seq_len(n)[-i]; zj <- z[j_idx]
      wij <- W[i, j_idx]
      xij  <- if (has_x)    X[i, j_idx]     else NULL
      xbij <- if (has_xbin) X_bin[i, j_idx] else NULL
      eta_j <- eta[j_idx]

      logp <- vapply(seq_len(K), function(k) {
        xt_ij_k <- if (has_x) ifelse(xij > 0, xij - xbar_cis[k, zj], 0) else NULL
        mu_ij <- .mu_pair_vec(rep(k, length(zj)), zj, eta[i] + eta_j, xt_ij_k, xbij)
        log(pi_k[k] + 1e-300) - sum((wij - mu_ij)^2) / (2 * tau2)
      }, numeric(1L))

      logp <- logp - max(logp)
      p <- exp(logp); p <- p / sum(p)
      z[i] <- sample.int(K, 1L, prob = p)
    }
    if (any(tabulate(z, nbins = K) < 1L)) next

    ## 5b. optional reallocate move
    if (reallocate) {
      ra <- propose_reallocate(z, xbar_cis)
      realloc_accept[iter] <- ra$accepted
      if (ra$accepted) z <- ra$z
    }

    ## 6. save
    if (iter > burn && ((iter - burn) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      z_save[save_idx, ]   <- z
      eta_save[save_idx, ] <- eta
      gamma_save[save_idx, ] <- .sym_mat_to_vec(G, K)
      if (has_x)    beta_save[save_idx, ]  <- .sym_mat_to_vec(B_mat, K)
      if (has_xbin) delta_save[save_idx, ] <- .sym_mat_to_vec(D_mat, K)
      tau2_save[save_idx] <- tau2
      loglik_save[save_idx] <- compute_loglik(xt_vec)
    }

    if (verbose && iter %% 50L == 0L) {
      msg <- sprintf("iter %d  tau2=%.4f  sizes=%s", iter, tau2, paste(tabulate(z, nbins = K), collapse = ","))
      if (reallocate) {
        recent <- realloc_accept[max(1L, iter - 49L):iter]
        msg <- paste0(msg, sprintf("  realloc_accepts(last50)=%d/%d", sum(recent), length(recent)))
      }
      message(msg)
    }
  }

  if (save_idx == 0L) {
    stop("No saved draw had all ", K, " blocks occupied (every post-burn iteration hit an ",
         "empty block and was skipped). Try a smaller K, a longer burn-in, or more iterations.")
  }
  if (save_idx < n_save) {
    warning(save_idx, " of ", n_save, " scheduled draws were saved; the rest were skipped ",
            "because a block was empty at save time. Summaries below are based on the ",
            save_idx, " saved draws only.")
    z_save      <- z_save[seq_len(save_idx), , drop = FALSE]
    eta_save    <- eta_save[seq_len(save_idx), , drop = FALSE]
    gamma_save  <- gamma_save[seq_len(save_idx), , drop = FALSE]
    if (has_x)    beta_save  <- beta_save[seq_len(save_idx), , drop = FALSE]
    if (has_xbin) delta_save <- delta_save[seq_len(save_idx), , drop = FALSE]
    tau2_save   <- tau2_save[seq_len(save_idx)]
    loglik_save <- loglik_save[seq_len(save_idx)]
  }

  z_hat <- apply(z_save, 2L, function(x) { tb <- table(x); as.integer(names(tb)[which.max(tb)]) })

  list(
    has_x = has_x, has_xbin = has_xbin, reallocate = reallocate,
    posterior = list(
      z_save = z_save, eta_save = eta_save, gamma_save = gamma_save,
      beta_save = beta_save, delta_save = delta_save,
      tau2_save = tau2_save, loglik_save = loglik_save,
      realloc_accept = realloc_accept
    ),
    summary = list(
      z_hat = z_hat,
      block_sizes_hat = tabulate(z_hat, nbins = K),
      eta_post_mean = colMeans(eta_save),
      gamma_post_mean = .make_sym_mat(colMeans(gamma_save), K),
      beta_post_mean  = if (has_x)    .make_sym_mat(colMeans(beta_save), K)  else NULL,
      delta_post_mean = if (has_xbin) .make_sym_mat(colMeans(delta_save), K) else NULL,
      tau2_post_mean = mean(tau2_save),
      mean_loglik = mean(loglik_save),
      max_loglik = max(loglik_save),
      realloc_accept_rate = if (reallocate) mean(realloc_accept) else NULL
    )
  )
}
