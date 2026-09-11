test_that("multinomial score scaling agrees with a literal expanded sample", {
  scores <- cbind(c(-4, 0, 3, 7), c(2, -1, 8, 5))
  counts <- c(2, 0, 3, 1)
  index <- rep(seq_len(nrow(scores)), counts)
  expanded <- scores[index, , drop = FALSE]
  # scale() uses the sample SD. Convert it to the population-moment
  # convention using the literal expanded sample size.
  expected <- scale(expanded) * sqrt(nrow(expanded) / (nrow(expanded) - 1))
  actual <- wdsmatch:::wdsm_standardize(scores, counts)
  expect_equal(as.vector(actual[index, ]), as.vector(expected), tolerance = 1e-12)
  expect_equal(colMeans(actual[index, ]), c(0, 0), tolerance = 1e-12)
  expect_equal(colMeans(actual[index, ]^2), c(1, 1), tolerance = 1e-12)
  omitted_changed <- scores
  omitted_changed[2, ] <- c(1000, -1000)
  changed <- wdsmatch:::wdsm_standardize(omitted_changed, counts)
  expect_equal(changed[index, ], actual[index, ], tolerance = 1e-12)
  expect_error(wdsmatch:::wdsm_standardize(scores, c(0, 0, 0, 0)),
               "Invalid scores or multiplicities")
  expect_error(wdsmatch:::wdsm_standardize(cbind(rep(1, 4), 1:4), counts),
               "Degenerate matching score")
})

test_that("exact distance ties select donors in original row order", {
  Y <- c(2, 10, 6, 12, 8, 4)
  Z <- c(0, 1, 0, 1, 0, 1)
  coordinate <- c(-1, 0, 1, 0, 1, -1)
  scores <- list(D0 = cbind(coordinate, coordinate),
                 D1 = cbind(coordinate, coordinate),
                 q0 = rep(0, 6), q1 = rep(0, 6))
  fit <- wdsmatch:::wdsm_point(Y, Z, rep(1, 6), scores, 1L, "PATE")
  # These donor expectations are specified by original row number.
  # Units 2 and 4 are equidistant from all three controls; units 3 and 5
  # are equidistant from treated units 2 and 4.
  expect_identical(fit$matches_0[c(2, 4, 6)], list(1L, 1L, 1L))
  expect_identical(fit$matches_1[c(1, 3, 5)], list(6L, 2L, 2L))
  expect_equal(fit$estimate, 14 / 3, tolerance = 1e-12)
  fit_two <- wdsmatch:::wdsm_point(Y, Z, rep(1, 6), scores, 2L, "PATE")
  expect_identical(fit_two$matches_0[[2]], c(1L, 3L))
  expect_identical(fit_two$matches_1[[3]], c(2L, 4L))
})

test_that("propensity fitting rejects unidentified and saturated solutions", {
  singular <- data.frame(Z = rep(c(0, 1), 10), x = rep(-4:5, each = 2))
  singular$duplicate <- singular$x
  expect_error(wdsmatch:::wdsm_fit_ps(Z ~ x + duplicate, singular, rep(1, 20)),
               "unidentified or nonfinite coefficients")
  separated <- data.frame(Z = rep(c(0, 1), each = 10),
                          x = c(-10:-1, 1:10))
  expect_error(suppressWarnings(wdsmatch:::wdsm_fit_ps(
    Z ~ x, separated, rep(1, 20))),
    "numerically saturated probabilities|did not converge")
  expect_error(wdsmatch:::wdsm_fit_ps(
    Z ~ x, singular, as.numeric(singular$Z == 0)),
    "both positive-weight treatment arms")
})

test_that("regression failures cannot silently replace nuisance predictions", {
  d <- data.frame(Y = c(1, 3, 2, 5, 4, 7), x = 1:6, duplicate = 2 * (1:6))
  expect_error(wdsmatch:::wdsm_fit_regression(
    Y ~ x + duplicate, d, d, rep(1, 6), "Control prognostic"),
    "Control prognostic fit has unidentified or nonfinite coefficients")
  expect_error(wdsmatch:::wdsm_fit_regression(
    Y ~ x, d, d, rep(0, 6), "Control prognostic"), "positive total")
  missing <- d
  missing$Y[1] <- NA_real_
  expect_error(wdsmatch:::wdsm_fit_regression(
    Y ~ x, missing, d, rep(1, 6), "Control prognostic"), "missing values")
})
