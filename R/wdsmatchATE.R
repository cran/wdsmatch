#' Weighted Double Score Matching Estimator for Population Average Treatment Effect
#'
#' Estimates the population average treatment effect (PATE) using weighted
#' double score matching (WDSM) with survey design weights. The method
#' matches treated and control units on arm-specific double scores
#' D_z(X) = (e(X), Psi_z(X)) for z in \{0,1\}, imputes missing potential
#' outcomes via survey-weighted averaging within match sets, and aggregates
#' using Hajek normalization. Polynomial sieve bias correction removes the
#' finite-sample matching discrepancy. Variance estimation uses a
#' linearization-based multinomial bootstrap that re-estimates score
#' parameters while preserving the original matching structure and
#' survey-weighted reuse frequencies.
#'
#' @details
#' The estimator achieves double robustness: it is consistent when either the
#' propensity score model or the prognostic score model is correctly specified.
#'
#' Under retrospective sampling (\code{sampling = "retrospective"}), the
#' propensity score is estimated with survey weights to recover the
#' population-level treatment assignment mechanism. Under prospective sampling
#' (\code{sampling = "prospective"}), the propensity score is estimated
#' without survey weights. Prognostic scores are always estimated without
#' survey weights, as the conditional outcome mean is invariant to the
#' sampling design.
#'
#' The sieve basis uses log-odds of the propensity score to match the
#' coordinate system used in matching distance computation.
#'
#' @param Y Numeric vector of observed outcomes.
#' @param X Numeric matrix or data frame of covariates.
#' @param Z Binary treatment assignment indicator (1 = treated, 0 = control).
#' @param weights Numeric vector of survey design weights. Required.
#' @param M Number of nearest neighbors for matching (default 5).
#' @param ps Numeric vector of pre-estimated propensity scores. If
#'   \code{NULL} (default), estimated internally using \code{model.ps}.
#' @param pg Numeric matrix of pre-estimated prognostic scores with columns
#'   \code{psi0} (control) and \code{psi1} (treated). If \code{NULL}
#'   (default), estimated internally using \code{model.pg}.
#' @param model.ps Formula for propensity score model (e.g.,
#'   \code{Z ~ X1 + X2}). If \code{NULL}, uses all columns of \code{X}.
#' @param model.pg Formula for prognostic score model (e.g.,
#'   \code{Y ~ X1 + X2}). If \code{NULL}, uses all columns of \code{X}.
#' @param sampling Character: \code{"retrospective"} (default) for
#'   treatment-dependent sampling (survey-weighted PS estimation), or
#'   \code{"prospective"} for treatment-independent sampling
#'   (unweighted PS estimation).
#' @param use.bias.correction Logical: apply polynomial sieve bias correction
#'   (default \code{TRUE}).
#' @param varest Logical: compute bootstrap variance estimate and confidence
#'   interval (default \code{TRUE}).
#' @param boots Number of multinomial bootstrap replicates (default 200).
#' @param alpha Significance level for confidence intervals (default 0.05).
#'
#' @return A list with components:
#'   \item{estimate}{Point estimate of PATE.}
#'   \item{se}{Bootstrap standard error (if \code{varest = TRUE}).}
#'   \item{ci}{Confidence interval as \code{c(lower, upper)}
#'     (if \code{varest = TRUE}).}
#'   \item{boot.estimates}{Vector of bootstrap replicate estimates
#'     (if \code{varest = TRUE}).}
#'   \item{M}{Number of matches used.}
#'   \item{n}{Sample size.}
#'   \item{n.treated}{Number of treated units.}
#'   \item{n.control}{Number of control units.}
#'   \item{call}{The matched call.}
#'
#' @examples
#' data(survey_obs)
#' fit <- wdsmatchATE(
#'   Y = survey_obs$Y,
#'   X = survey_obs[, c("X1","X2","X3","X4","X5","X6")],
#'   Z = survey_obs$Z,
#'   weights = survey_obs$survey_weight,
#'   M = 3,
#'   model.ps = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
#'   model.pg = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
#'   sampling = "retrospective",
#'   varest = FALSE
#' )
#' fit
#'
#' @export
wdsmatchATE <- function(Y, X, Z, weights, M = 5,
                        ps = NULL, pg = NULL,
                        model.ps = NULL, model.pg = NULL,
                        sampling = c("retrospective", "prospective"),
                        use.bias.correction = TRUE,
                        varest = TRUE, boots = 200, alpha = 0.05) {
  cl <- match.call()
  sampling <- match.arg(sampling)

  if (missing(weights) || is.null(weights)) stop("'weights' are required for WDSM.")
  if (!is.numeric(Y)) stop("'Y' must be numeric.")
  if (!all(Z %in% c(0, 1))) stop("'Z' must be binary (0/1).")
  if (length(Y) != length(Z)) stop("'Y' and 'Z' must have the same length.")

  X <- as.data.frame(X)
  n <- length(Y)

  scores <- estimate_scores(Y, X, Z, sw = weights, ps = ps, pg = pg,
                            model.ps = model.ps, model.pg = model.pg,
                            sampling = sampling)

  pt <- wdsm_match_ate(Y, Z, weights, scores$ps_logit,
                       scores$psi0, scores$psi1, M, use.bias.correction)
  pt$X_internal <- X

  result <- list(
    estimate = pt$estimate,
    se = NA_real_, ci = c(NA_real_, NA_real_),
    boot.estimates = NULL,
    M = M, n = n,
    n.treated = sum(Z == 1), n.control = sum(Z == 0),
    call = cl
  )

  if (varest) {
    if (is.null(model.ps)) {
      xnames <- colnames(X)
      model.ps <- stats::as.formula(paste("Z ~", paste(xnames, collapse = " + ")))
    }
    if (is.null(model.pg)) {
      xnames <- colnames(X)
      model.pg <- stats::as.formula(paste("Y ~", paste(xnames, collapse = " + ")))
    }
    boot <- wdsm_bootstrap_ate(pt, boots = boots, alpha = alpha,
                               model.ps = model.ps, model.pg = model.pg,
                               sampling = sampling)
    result$se <- boot$se
    result$ci <- boot$ci
    result$boot.estimates <- boot$boot_estimates
  }

  class(result) <- "wdsmatch"
  result
}
