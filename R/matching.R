#' @keywords internal
wdsm_match_ate <- function(Y, Z, weights, ps_logit, psi0, psi1, M,
                           use.bias.correction = TRUE) {
  n <- length(Y)
  treated <- which(Z == 1)
  control <- which(Z == 0)

  D0 <- cbind(ps_logit, psi0)
  D1 <- cbind(ps_logit, psi1)

  D0_std <- standardize_scores(D0, weights)
  D1_std <- standardize_scores(D1, weights)

  matches_0 <- vector("list", n)
  matches_1 <- vector("list", n)
  K0_w <- K1_w <- numeric(n)

  for (i in treated) {
    dists <- sqrt(rowSums((sweep(D0_std[control, , drop = FALSE], 2,
                                  D0_std[i, ]))^2))
    nbrs <- find_M_nearest(dists, control, M)
    matches_0[[i]] <- nbrs
    nbr_w <- weights[nbrs]
    for (j in nbrs) K0_w[j] <- K0_w[j] + weights[i] / sum(nbr_w)
  }

  for (i in control) {
    dists <- sqrt(rowSums((sweep(D1_std[treated, , drop = FALSE], 2,
                                  D1_std[i, ]))^2))
    nbrs <- find_M_nearest(dists, treated, M)
    matches_1[[i]] <- nbrs
    nbr_w <- weights[nbrs]
    for (j in nbrs) K1_w[j] <- K1_w[j] + weights[i] / sum(nbr_w)
  }

  sieve_x <- build_sieve_basis(psi0, psi1, ps_logit)
  if (use.bias.correction) {
    bc <- sieve_bias_correct(Y, Z, sieve_x, weights)
    mu0_bc <- bc$mu0_bc
    mu1_bc <- bc$mu1_bc
  } else {
    mu0_bc <- psi0
    mu1_bc <- psi1
  }

  Yhat0 <- Y
  Yhat1 <- Y
  for (i in treated) {
    nbrs <- matches_0[[i]]
    nbr_w <- weights[nbrs]
    Yhat0[i] <- sum(nbr_w * Y[nbrs]) / sum(nbr_w)
    Yhat1[i] <- Y[i]
  }
  for (i in control) {
    nbrs <- matches_1[[i]]
    nbr_w <- weights[nbrs]
    Yhat1[i] <- sum(nbr_w * Y[nbrs]) / sum(nbr_w)
    Yhat0[i] <- Y[i]
  }

  if (use.bias.correction) {
    for (i in treated) {
      nbrs <- matches_0[[i]]
      nbr_w <- weights[nbrs]
      Yhat0[i] <- sum(nbr_w * Y[nbrs]) / sum(nbr_w) +
                  mu0_bc[i] - sum(nbr_w * mu0_bc[nbrs]) / sum(nbr_w)
    }
    for (i in control) {
      nbrs <- matches_1[[i]]
      nbr_w <- weights[nbrs]
      Yhat1[i] <- sum(nbr_w * Y[nbrs]) / sum(nbr_w) +
                  mu1_bc[i] - sum(nbr_w * mu1_bc[nbrs]) / sum(nbr_w)
    }
  }

  estimate <- sum(weights * (Yhat1 - Yhat0)) / sum(weights)

  list(estimate = estimate,
       Y = Y, Z = Z, X_internal = NULL, weights = weights,
       e_score = expit(ps_logit), psi0 = psi0, psi1 = psi1,
       ps_logit = ps_logit,
       mu0_bc = mu0_bc, mu1_bc = mu1_bc,
       matches_0 = matches_0, matches_1 = matches_1,
       K0_w = K0_w, K1_w = K1_w,
       sieve_x = sieve_x)
}

#' @keywords internal
wdsm_match_att <- function(Y, Z, weights, ps_logit, psi0, psi1, M,
                           use.bias.correction = TRUE) {
  n <- length(Y)
  treated <- which(Z == 1)
  control <- which(Z == 0)

  D0 <- cbind(ps_logit, psi0)
  D0_std <- standardize_scores(D0, weights)

  matches_0 <- vector("list", n)
  K0_w <- numeric(n)

  for (i in treated) {
    dists <- sqrt(rowSums((sweep(D0_std[control, , drop = FALSE], 2,
                                  D0_std[i, ]))^2))
    nbrs <- find_M_nearest(dists, control, M)
    matches_0[[i]] <- nbrs
    nbr_w <- weights[nbrs]
    for (j in nbrs) K0_w[j] <- K0_w[j] + weights[i] * weights[j] / sum(nbr_w)
  }

  sieve_x <- build_sieve_basis(psi0, psi1, ps_logit)
  if (use.bias.correction) {
    bc <- sieve_bias_correct(Y, Z, sieve_x, weights)
    mu0_bc <- bc$mu0_bc
    mu1_bc <- bc$mu1_bc
  } else {
    mu0_bc <- psi0
    mu1_bc <- psi1
  }

  Yhat0 <- numeric(n)
  for (i in treated) {
    nbrs <- matches_0[[i]]
    nbr_w <- weights[nbrs]
    if (use.bias.correction) {
      Yhat0[i] <- sum(nbr_w * Y[nbrs]) / sum(nbr_w) +
                  mu0_bc[i] - sum(nbr_w * mu0_bc[nbrs]) / sum(nbr_w)
    } else {
      Yhat0[i] <- sum(nbr_w * Y[nbrs]) / sum(nbr_w)
    }
  }

  W <- weights
  estimate <- sum(W[treated] * (Y[treated] - Yhat0[treated])) / sum(W[treated])

  list(estimate = estimate,
       Y = Y, Z = Z, X_internal = NULL, weights = weights,
       e_score = expit(ps_logit), psi0 = psi0, psi1 = psi1,
       ps_logit = ps_logit,
       mu0_bc = mu0_bc, mu1_bc = mu1_bc,
       matches_0 = matches_0, K0_w = K0_w,
       sieve_x = sieve_x)
}
