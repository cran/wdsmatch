#' Weighted Double Score Matching Estimator for Population Average Treatment Effect on the Treated
#'
#' Estimates the population average treatment effect on the treated (PATT)
#' using weighted double score matching (WDSM) with survey design weights.
#' Performs one-sided matching from treated to control units on the
#' control-side double score D_0(X) = (e(X), Psi_0(X)). Only the
#' counterfactual control outcome Y(0) needs imputation; treated outcomes
#' are directly observed. Aggregation uses Hajek normalization over
#' treated-side survey weights. Polynomial sieve bias correction and
#' linearization-based multinomial bootstrap with survey-weighted reuse
#' frequencies are applied. PATT requires only one-sided
#' unconfoundedness: Y(0) independent of Z given X.
#'
#' @inheritParams wdsmatchATE
#'
#' @return A list with components:
#'   \item{estimate}{Point estimate of PATT.}
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
#' fit <- wdsmatchATT(
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
wdsmatchATT <- function(Y, X, Z, weights, M = 5,
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

  X <- as.data.frame(X)
  n <- length(Y)

  scores <- estimate_scores(Y, X, Z, sw = weights, ps = ps, pg = pg,
                            model.ps = model.ps, model.pg = model.pg,
                            sampling = sampling)

  pt <- wdsm_match_att(Y, Z, weights, scores$ps_logit,
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
    boot <- wdsm_bootstrap_att(pt, boots = boots, alpha = alpha,
                               model.ps = model.ps, model.pg = model.pg,
                               sampling = sampling)
    result$se <- boot$se
    result$ci <- boot$ci
    result$boot.estimates <- boot$boot_estimates
  }

  class(result) <- "wdsmatch"
  result
}
