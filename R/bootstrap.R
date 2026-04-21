utils::globalVariables(c(".ps_wts", ".pg_wts"))

#' @keywords internal
wdsm_bootstrap_ate <- function(point_est, boots = 200, alpha = 0.05,
                               model.ps = NULL, model.pg = NULL,
                               sampling = "retrospective") {
  Y <- point_est$Y
  Z <- point_est$Z
  w <- point_est$weights
  n <- length(Y)
  X_df <- point_est$X_internal

  K0_w <- point_est$K0_w
  K1_w <- point_est$K1_w
  matches_0 <- point_est$matches_0
  matches_1 <- point_est$matches_1

  boot_estimates <- numeric(boots)

  for (b in seq_len(boots)) {
    m_star <- as.numeric(stats::rmultinom(1, n, rep(1 / n, n)))
    w_boot <- m_star * w

    if (sum(w_boot[Z == 1]) == 0 || sum(w_boot[Z == 0]) == 0) next

    df_b <- data.frame(Z = Z, Y = Y, X_df)
    idx0 <- which(Z == 0); idx1 <- which(Z == 1)

    ps_fit <- eval(substitute(
      stats::glm(model.ps, data = df_b, family = stats::quasibinomial(link = "logit"),
                 weights = WTS), list(WTS = w_boot)))
    ps_star <- ps_fit$fitted.values

    pg0_fit <- eval(substitute(
      stats::lm(model.pg, data = df_b[idx0, ], weights = WTS),
      list(WTS = m_star[idx0])))
    mu0_star <- as.numeric(stats::predict(pg0_fit, newdata = df_b))

    pg1_fit <- eval(substitute(
      stats::lm(model.pg, data = df_b[idx1, ], weights = WTS),
      list(WTS = m_star[idx1])))
    mu1_star <- as.numeric(stats::predict(pg1_fit, newdata = df_b))

    ps_logit_star <- safe_logit(ps_star)
    sieve_star <- build_sieve_basis(mu0_star, mu1_star, ps_logit_star)

    bc0_fit <- eval(substitute(
      stats::lm(Y ~ ., data = data.frame(Y = Y[idx0], sieve_star[idx0, , drop = FALSE]),
                weights = WTS), list(WTS = w_boot[idx0])))
    mu0_bc_star <- as.numeric(stats::predict(bc0_fit, newdata = data.frame(sieve_star)))

    bc1_fit <- eval(substitute(
      stats::lm(Y ~ ., data = data.frame(Y = Y[idx1], sieve_star[idx1, , drop = FALSE]),
                weights = WTS), list(WTS = w_boot[idx1])))
    mu1_bc_star <- as.numeric(stats::predict(bc1_fit, newdata = data.frame(sieve_star)))

    tau_star <- mean(w_boot * mu1_bc_star) / mean(w_boot) -
                mean(w_boot * mu0_bc_star) / mean(w_boot)

    for (i in seq_len(n)) {
      if (Z[i] == 0) {
        nbrs <- matches_0[[i]]
        K_i <- K0_w[i]
        if (K_i > 0 && length(nbrs) > 0) {
          numer <- sum(w_boot[nbrs] * (Y[nbrs] - mu1_bc_star[nbrs]))
          denom <- sum(w_boot[nbrs])
          if (denom > 0) {
            residual <- numer / denom
            tau_star <- tau_star + w_boot[i] * (1 + K_i) * residual / sum(w_boot)
          }
        }
      } else {
        nbrs <- matches_1[[i]]
        K_i <- K1_w[i]
        if (K_i > 0 && length(nbrs) > 0) {
          numer <- sum(w_boot[nbrs] * (Y[nbrs] - mu0_bc_star[nbrs]))
          denom <- sum(w_boot[nbrs])
          if (denom > 0) {
            residual <- numer / denom
            tau_star <- tau_star - w_boot[i] * (1 + K_i) * residual / sum(w_boot)
          }
        }
      }
    }

    boot_estimates[b] <- tau_star
  }

  boot_estimates <- boot_estimates[!is.na(boot_estimates)]
  se <- stats::sd(boot_estimates)
  ci <- stats::quantile(boot_estimates, c(alpha / 2, 1 - alpha / 2))

  list(se = se, ci = as.numeric(ci), boot_estimates = boot_estimates)
}

#' @keywords internal
wdsm_bootstrap_att <- function(point_est, boots = 200, alpha = 0.05,
                               model.ps = NULL, model.pg = NULL,
                               sampling = "retrospective") {
  Y <- point_est$Y
  Z <- point_est$Z
  w <- point_est$weights
  n <- length(Y)
  X_df <- point_est$X_internal

  K0_patt <- point_est$K0_w
  matches_0 <- point_est$matches_0

  boot_estimates <- numeric(boots)

  for (b in seq_len(boots)) {
    m_star <- as.numeric(stats::rmultinom(1, n, rep(1 / n, n)))
    w_boot <- m_star * w

    if (sum(w_boot[Z == 1]) == 0 || sum(w_boot[Z == 0]) == 0) next

    df_b <- data.frame(Z = Z, Y = Y, X_df)
    idx0 <- which(Z == 0); idx1 <- which(Z == 1)

    ps_fit <- eval(substitute(
      stats::glm(model.ps, data = df_b, family = stats::quasibinomial(link = "logit"),
                 weights = WTS), list(WTS = w_boot)))
    ps_star <- ps_fit$fitted.values

    pg0_fit <- eval(substitute(
      stats::lm(model.pg, data = df_b[idx0, ], weights = WTS),
      list(WTS = m_star[idx0])))
    psi0_star <- as.numeric(stats::predict(pg0_fit, newdata = df_b))

    pg1_fit <- eval(substitute(
      stats::lm(model.pg, data = df_b[idx1, ], weights = WTS),
      list(WTS = m_star[idx1])))
    psi1_star <- as.numeric(stats::predict(pg1_fit, newdata = df_b))

    ps_logit_star <- safe_logit(ps_star)
    sieve_star <- build_sieve_basis(psi0_star, psi1_star, ps_logit_star)

    bc0_fit <- eval(substitute(
      stats::lm(Y ~ ., data = data.frame(Y = Y[idx0], sieve_star[idx0, , drop = FALSE]),
                weights = WTS), list(WTS = w_boot[idx0])))
    mu0_bc_star <- as.numeric(stats::predict(bc0_fit, newdata = data.frame(sieve_star)))

    bc1_fit <- eval(substitute(
      stats::lm(Y ~ ., data = data.frame(Y = Y[idx1], sieve_star[idx1, , drop = FALSE]),
                weights = WTS), list(WTS = w_boot[idx1])))
    mu1_bc_star <- as.numeric(stats::predict(bc1_fit, newdata = data.frame(sieve_star)))

    treated_sum <- sum(w_boot[idx1] * (Y[idx1] - mu0_bc_star[idx1]))
    control_sum <- sum(m_star[idx0] * K0_patt[idx0] * (Y[idx0] - mu0_bc_star[idx0]))
    norm_treated <- sum(w_boot[idx1])

    tau_star <- (treated_sum - control_sum) / norm_treated
    boot_estimates[b] <- tau_star
  }

  boot_estimates <- boot_estimates[!is.na(boot_estimates)]
  se <- stats::sd(boot_estimates)
  ci <- stats::quantile(boot_estimates, c(alpha / 2, 1 - alpha / 2))

  list(se = se, ci = as.numeric(ci), boot_estimates = boot_estimates)
}
