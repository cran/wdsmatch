test_that("wdsmatchATE replicates simulation-style estimate (CorCor, retrospective)", {
  skip_on_cran()

  pop_file <- system.file("extdata", "sim_pop_small.RData", package = "wdsmatch")
  if (!nzchar(pop_file)) skip("Simulation population data not bundled")

  load(pop_file)
  set.seed(2026)

  n_sample <- 1000
  N <- nrow(population)
  expit <- function(x) 1 / (1 + exp(-x))

  c0 <- log(n_sample/N / (1 - n_sample/N))
  formula_S <- c0 + log(0.9) * population$Z +
    log(1.05) * population$X1 + log(1.10) * population$X2 +
    log(1.15) * population$X3 + log(1.10) * population$X4 +
    log(1.05) * population$X5 + log(1.10) * population$X6
  prob_S <- expit(formula_S)
  S_ind <- rbinom(N, 1, prob_S)
  samp <- population[S_ind == 1, ]
  samp$survey_weight <- 1 / prob_S[S_ind == 1]

  X <- samp[, c("X1","X2","X3","X4","X5","X6")]
  fit_ate <- wdsmatchATE(
    Y = samp$Y_obs, X = X, Z = samp$Z,
    weights = samp$survey_weight, M = 5,
    model.ps = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
    model.pg = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
    sampling = "retrospective",
    use.bias.correction = TRUE, varest = FALSE
  )

  true_pate <- mean(population$Y_1) - mean(population$Y_0)
  expect_true(is.finite(fit_ate$estimate))
  expect_true(abs(fit_ate$estimate - true_pate) / abs(true_pate) < 0.20)

  fit_att <- wdsmatchATT(
    Y = samp$Y_obs, X = X, Z = samp$Z,
    weights = samp$survey_weight, M = 5,
    model.ps = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
    model.pg = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
    sampling = "retrospective",
    use.bias.correction = TRUE, varest = FALSE
  )

  true_patt <- mean(population$Y_1[population$Z == 1]) -
               mean(population$Y_0[population$Z == 1])
  expect_true(is.finite(fit_att$estimate))
  expect_true(abs(fit_att$estimate - true_patt) / abs(true_patt) < 0.20)
})
