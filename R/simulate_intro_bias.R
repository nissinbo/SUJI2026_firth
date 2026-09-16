# Independent introduction experiment. Run explicitly from the project root.
stopifnot(requireNamespace("logistf"), requireNamespace("lpSolve"))
set.seed(20261015)
intro_n <- 60L
intro_B <- 2000L
intro_truth <- c(0, 1, 1, 1) # intercept, x1, x2, x3
intro_estimates <- array(NA_real_, c(intro_B, 3L, 2L),
  dimnames = list(NULL, c("x1", "x2", "x3"), c("MLE", "Firth")))
intro_separated <- logical(intro_B)
for (b in seq_len(intro_B)) {
  dat <- as.data.frame(matrix(rnorm(intro_n * 3L), intro_n, 3L))
  names(dat) <- c("x1", "x2", "x3")
  X <- model.matrix(~ x1 + x2 + x3, dat)
  dat$y <- rbinom(intro_n, 1, plogis(drop(X %*% intro_truth)))
  # Nonnegative separating margins with a positive sum detect complete/quasi separation.
  A <- qr.Q(qr(X)) * (2 * dat$y - 1)
  A <- cbind(A, -A)
  sep <- lpSolve::lp("max", colSums(A), rbind(A, 1),
    c(rep(">=", intro_n), "<="), c(rep(0, intro_n), 1))
  stopifnot(sep$status == 0L)
  intro_separated[b] <- sep$objval > 1e-8
  if (intro_separated[b]) next # Paired comparison on the same finite-MLE samples.
  ml <- glm(y ~ x1 + x2 + x3, data = dat, family = binomial(),
    control = glm.control(epsilon = 1e-10, maxit = 100L))
  ff <- logistf::logistf(y ~ x1 + x2 + x3, data = dat, pl = FALSE,
    control = logistf::logistf.control(maxit = 100L,
      lconv = 1e-8, gconv = 1e-8, xconv = 1e-8))
  stopifnot(ml$converged, !ml$boundary, all(is.finite(coef(ml))),
    all(is.finite(coef(ff))), max(abs(ff$conv)) < 1e-6)
  intro_estimates[b, , "MLE"] <- coef(ml)[-1]
  intro_estimates[b, , "Firth"] <- coef(ff)[-1]
}

dir.create("results/raw", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(n = intro_n, repetitions = intro_B, truth = intro_truth,
  estimates = intro_estimates, separated = intro_separated, seed = 20261015L,
  versions = c(R = R.version.string, logistf = as.character(packageVersion("logistf")),
    lpSolve = as.character(packageVersion("lpSolve")))), "results/raw/intro_bias.rds")
cat("Saved", intro_B, "replications; separated:", sum(intro_separated), "\n")
