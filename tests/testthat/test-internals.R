test_that("safe_logit handles boundary values", {
  expect_true(is.finite(wdsmatch:::safe_logit(0)))
  expect_true(is.finite(wdsmatch:::safe_logit(1)))
  expect_true(is.finite(wdsmatch:::safe_logit(0.5)))
  expect_equal(wdsmatch:::safe_logit(0.5), 0, tolerance = 1e-6)
})

test_that("expit is inverse of logit", {
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  expect_equal(wdsmatch:::expit(wdsmatch:::safe_logit(p)), p, tolerance = 1e-4)
})

test_that("standardize_scores produces unit variance", {
  set.seed(1)
  S <- matrix(rnorm(200), ncol = 2)
  S_std <- wdsmatch:::standardize_scores(S)
  expect_equal(apply(S_std, 2, mean), c(0, 0), tolerance = 1e-10)
  expect_equal(apply(S_std, 2, sd), c(1, 1), tolerance = 0.1)
})

test_that("standardize_scores respects weights", {
  set.seed(2)
  S <- matrix(rnorm(200), ncol = 2)
  w <- runif(100, 1, 5)
  S_std <- wdsmatch:::standardize_scores(S, w)
  wmean <- colSums(w * S_std) / sum(w)
  expect_equal(wmean, c(0, 0), tolerance = 1e-10)
})

test_that("build_sieve_basis ATE has 10 columns", {
  psi0 <- rnorm(50); psi1 <- rnorm(50); ps <- rnorm(50)
  basis <- wdsmatch:::build_sieve_basis(psi0, psi1, ps)
  expect_equal(ncol(basis), 10)
  expect_equal(nrow(basis), 50)
})

test_that("build_sieve_basis ATT has 5 columns", {
  psi0 <- rnorm(50); ps <- rnorm(50)
  basis <- wdsmatch:::build_sieve_basis(psi0, NULL, ps)
  expect_equal(ncol(basis), 5)
})

test_that("find_M_nearest returns correct count", {
  dists <- c(3, 1, 4, 1, 5)
  idx <- 1:5
  expect_equal(length(wdsmatch:::find_M_nearest(dists, idx, 3)), 3)
  expect_equal(length(wdsmatch:::find_M_nearest(dists, idx, 1)), 1)
  expect_equal(wdsmatch:::find_M_nearest(dists, idx, 1), 2)
})

test_that("sieve_bias_correct returns predictions for all units", {
  set.seed(3)
  n <- 100
  Y <- rnorm(n)
  Z <- rep(c(0, 1), each = 50)
  sx <- matrix(rnorm(n * 5), ncol = 5)
  bc <- wdsmatch:::sieve_bias_correct(Y, Z, sx)
  expect_equal(length(bc$mu0_bc), n)
  expect_equal(length(bc$mu1_bc), n)
})
