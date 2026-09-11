# wdsmatch 0.2.0

This release corrects the numerical implementation. Point estimates, standard
errors, and confidence intervals can differ from 0.1.1. Analyses made with
0.1.1 should be rerun before numerical results are reused.

* Match on standardized propensity probabilities and arm-specific prognostic
  scores. Replace the old logit and joint-score polynomial with a complete
  quadratic basis in each arm's own two-dimensional double score.
* Correct PATE's replication residual contributions using the original weighted
  matching reuse coefficients. Use the corresponding fixed-reuse PATT formula.
* Apply survey-weighted propensity fitting for retrospective sampling and
  unweighted propensity fitting for prospective sampling, including replication
  refits. Use finite-checked, zero-initialized logistic fits.
* Use centered normal Wald intervals based on replication variance with divisor
  B, replacing the previous percentile intervals and B-1 variance convention.
* Respect disabled bias correction in both the point estimate and replication.
  This option retains matching discrepancies and does not provide the same
  asymptotic bias-correction guarantee.
* Hold supplied propensity/prognostic scores fixed in replication; corresponding
  standard errors are conditional on those supplied scores. PATT only requires
  the control-side prognostic model.
* Reject invalid inputs and unsuccessful numerical fits explicitly. Do not
  insert zero replicates or silently discard failed replicates.
* Retain the existing public arguments and result components, and add variance,
  confidence-level and procedure metadata. Print the requested confidence level
  rather than always labeling an interval as 95%.
* Add regression checks for the corrected calculations, improve the method and
  inference documentation, and add checks on Linux, Windows, and macOS.

# wdsmatch 0.1.1

* Add missing return-value documentation for the print and summary methods for
  CRAN resubmission. No numerical estimator changes from 0.1.0.

# wdsmatch 0.1.0

* Initial package version.
* `wdsmatchATE()`: Population average treatment effect estimation.
* `wdsmatchATT()`: Population average treatment effect on the treated.
* Supports retrospective and prospective sampling designs.
* Polynomial sieve bias correction with the original 10-term logit basis.
* Linearization-based multinomial bootstrap variance estimation.
* Bundled example dataset `survey_obs`.
