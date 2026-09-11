wdsm_point <- function(Y, Z, weights, scores, M, estimand) {
  n <- length(Y)
  treated <- which(Z == 1)
  control <- which(Z == 0)
  matches_0 <- matches_1 <- vector("list", n)
  for (i in treated) {
    distance <- (scores$D0[control, 1L] - scores$D0[i, 1L])^2 +
                (scores$D0[control, 2L] - scores$D0[i, 2L])^2
    matches_0[[i]] <- control[order(distance, method = "radix")[seq_len(M)]]
  }
  if (estimand == "PATE") for (i in control) {
    distance <- (scores$D1[treated, 1L] - scores$D1[i, 1L])^2 +
                (scores$D1[treated, 2L] - scores$D1[i, 2L])^2
    matches_1[[i]] <- treated[order(distance, method = "radix")[seq_len(M)]]
  }
  K <- wdsm_reuse(Z, weights, matches_0, matches_1, estimand)
  contrast <- rep(NA_real_, n)
  for (i in treated) {
    donors <- matches_0[[i]]
    imputed <- scores$q0[i] +
      sum(weights[donors] * (Y[donors] - scores$q0[donors])) / sum(weights[donors])
    contrast[i] <- Y[i] - imputed
  }
  if (estimand == "PATE") for (i in control) {
    donors <- matches_1[[i]]
    imputed <- scores$q1[i] +
      sum(weights[donors] * (Y[donors] - scores$q1[donors])) / sum(weights[donors])
    contrast[i] <- imputed - Y[i]
  }
  targets <- if (estimand == "PATE") seq_len(n) else treated
  estimate <- sum(weights[targets] * contrast[targets]) / sum(weights[targets])
  check <- wdsm_replicate(Y, Z, weights, K, rep(1, n), scores$q0, scores$q1, estimand)
  identity_error <- abs(estimate - check)
  if (!is.finite(identity_error) || identity_error > 1e-9 * (1 + abs(estimate)))
    stop("Direct imputation and original-reuse representations disagree", call. = FALSE)
  list(estimate = estimate, Y = Y, Z = Z, weights = weights,
       K = K, matches_0 = matches_0, matches_1 = matches_1,
       scores = scores, estimand = estimand, identity_error = identity_error)
}
