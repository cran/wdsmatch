#' @keywords internal
estimate_scores <- function(Y, X, Z, sw = NULL, ps = NULL, pg = NULL,
                            model.ps = NULL, model.pg = NULL,
                            sampling = "retrospective") {
  n <- length(Y)
  df <- data.frame(Z = Z, Y = Y, X)

  if (is.null(ps)) {
    if (is.null(model.ps)) {
      xnames <- colnames(X)
      model.ps <- stats::as.formula(paste("Z ~", paste(xnames, collapse = " + ")))
    }
    use_sw <- (sampling == "retrospective" && !is.null(sw))
    if (use_sw) {
      ps_fit <- eval(substitute(
        stats::glm(model.ps, data = df, family = stats::quasibinomial(link = "logit"),
                   weights = WTS),
        list(WTS = sw)
      ))
    } else {
      ps_fit <- stats::glm(model.ps, data = df, family = stats::quasibinomial(link = "logit"))
    }
    ps <- ps_fit$fitted.values
  }

  if (is.null(pg)) {
    if (is.null(model.pg)) {
      xnames <- colnames(X)
      model.pg <- stats::as.formula(paste("Y ~", paste(xnames, collapse = " + ")))
    }
    psi0 <- rep(NA_real_, n)
    psi1 <- rep(NA_real_, n)
    pg0_fit <- stats::lm(model.pg, data = df[Z == 0, ])
    psi0 <- stats::predict(pg0_fit, newdata = df)
    pg1_fit <- stats::lm(model.pg, data = df[Z == 1, ])
    psi1 <- stats::predict(pg1_fit, newdata = df)
    pg <- cbind(psi0 = as.numeric(psi0), psi1 = as.numeric(psi1))
  } else {
    if (is.vector(pg)) stop("pg must be a 2-column matrix (psi0, psi1)")
    if (ncol(pg) == 1) pg <- cbind(pg, pg)
  }

  list(ps = as.numeric(ps), psi0 = pg[, 1], psi1 = pg[, 2],
       ps_logit = safe_logit(ps))
}
