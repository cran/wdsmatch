#' @keywords internal
expit <- function(x) {
  1 / (1 + exp(-x))
}

#' @keywords internal
safe_logit <- function(p, eps = 0.001) {
  p <- pmax(pmin(p, 1 - eps), eps)
  log(p / (1 - p))
}

#' @keywords internal
standardize_scores <- function(S, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(S))
  w_sum <- sum(weights)
  S_mean <- colSums(weights * S) / w_sum
  S_centered <- sweep(S, 2, S_mean, "-")
  S_var <- colSums(weights * S_centered^2) / w_sum
  S_sd <- sqrt(S_var)
  S_sd[S_sd < 1e-10] <- 1
  sweep(S_centered, 2, S_sd, "/")
}

#' @keywords internal
find_M_nearest <- function(distances, candidate_idx, M) {
  ord <- order(distances)
  k <- min(M, length(ord))
  candidate_idx[ord[seq_len(k)]]
}

#' @keywords internal
build_sieve_basis <- function(psi0, psi1 = NULL, ps_logit) {
  if (!is.null(psi1)) {
    cbind(psi0, psi1, psi0^2, psi1^2, psi0 * psi1,
          ps_logit, ps_logit^2, ps_logit * psi0,
          ps_logit * psi1, ps_logit * psi0 * psi1)
  } else {
    cbind(psi0, psi0^2, ps_logit, ps_logit^2, ps_logit * psi0)
  }
}

#' @keywords internal
sieve_bias_correct <- function(Y, Z, sieve_x, weights = NULL) {
  colnames(sieve_x) <- paste0("s", seq_len(ncol(sieve_x)))
  w0 <- if (!is.null(weights)) weights[Z == 0] else NULL
  w1 <- if (!is.null(weights)) weights[Z == 1] else NULL

  bc0_fit <- stats::lm(Y ~ ., data = data.frame(Y = Y[Z == 0], sieve_x[Z == 0, , drop = FALSE]),
                        weights = w0)
  mu0_bc <- stats::predict(bc0_fit, newdata = data.frame(sieve_x))

  bc1_fit <- stats::lm(Y ~ ., data = data.frame(Y = Y[Z == 1], sieve_x[Z == 1, , drop = FALSE]),
                        weights = w1)
  mu1_bc <- stats::predict(bc1_fit, newdata = data.frame(sieve_x))

  list(mu0_bc = as.numeric(mu0_bc), mu1_bc = as.numeric(mu1_bc))
}
