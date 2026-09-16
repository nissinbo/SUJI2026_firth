# Model settings used in the presentation experiments.
finite_all <- function(x) length(x) > 0L && all(is.finite(x))

lf_variant <- function(id, fit = "NR", maxit = 25L, fit_tight = FALSE, pl_tight = FALSE) {
  ctl <- logistf::logistf.control(fit = fit, maxit = maxit)
  if (fit_tight) ctl[c("lconv", "gconv", "xconv")] <- list(1e-12, 1e-12, 1e-12)
  plctl <- logistf::logistpl.control()
  if (pl_tight) plctl[c("maxit", "lconv", "xconv")] <- list(200L, 1e-12, 1e-12)
  list(id = id, implementation = "logistf", algorithm = fit,
       initial_mode = "all_coefficients_zero", control = ctl, plcontrol = plctl)
}

bg_variant <- function(id, type, maxit = 100L, epsilon = 1e-6) {
  list(id = id, implementation = "brglm2", algorithm = type,
       initial_mode = "package_default_start_NULL", plcontrol = NULL,
       control = brglm2::brglmControl(type = type, a = 0.5, maxit = maxit, epsilon = epsilon))
}

variants <- list(
  lf_variant("lf_nr_default"),
  lf_variant("lf_irls_default", fit = "IRLS"),
  lf_variant("lf_nr_maxit5000", maxit = 5000L),
  lf_variant("lf_nr_fit_tight", maxit = 5000L, fit_tight = TRUE),
  lf_variant("lf_nr_pl_tight", pl_tight = TRUE),
  lf_variant("lf_nr_fit_pl_tight", maxit = 5000L, fit_tight = TRUE, pl_tight = TRUE),
  lf_variant("lf_irls_fit_pl_tight", fit = "IRLS", maxit = 5000L, fit_tight = TRUE, pl_tight = TRUE),
  bg_variant("bg_asmean_default", "AS_mean"),
  bg_variant("bg_asmean_maxit5000", "AS_mean", maxit = 5000L),
  bg_variant("bg_asmean_epsilon12", "AS_mean", maxit = 5000L, epsilon = 1e-12),
  bg_variant("bg_mpl_default", "MPL_Jeffreys"),
  bg_variant("bg_mpl_maxit5000", "MPL_Jeffreys", maxit = 5000L),
  bg_variant("bg_mpl_epsilon12", "MPL_Jeffreys", maxit = 5000L, epsilon = 1e-12)
)

run_logistf <- function(dat, meta, v, plconf = 2L) {
  form <- model_formula(meta)
  logistf::logistf(form, data = dat, firth = TRUE, pl = TRUE, alpha = .05,
    init = rep(0, ncol(model.matrix(form, dat))), control = v$control,
    plcontrol = v$plcontrol, plconf = plconf)
}

fit_one <- function(meta, v, dat) {
  timing <- proc.time()
  cap <- with_conditions(if (v$implementation == "logistf") run_logistf(dat, meta, v)
    else stats::glm(model_formula(meta), data = dat, family = binomial("logit"),
                   method = brglm2::brglmFit, control = v$control))
  elapsed <- proc.time() - timing
  audit <- data.frame(
    config_id = meta$config_id, scenario = meta$scenario, n = meta$n,
    n_exposed = meta$n_exposed, n_unexposed = meta$n_unexposed,
    a = meta$a, b = meta$b, c = meta$c, d = meta$d,
    separation = meta$separation, sparsity_class = meta$sparsity_class,
    variant_id = v$id, implementation = v$implementation, algorithm = v$algorithm,
    initial_mode = v$initial_mode, fit_id = paste(meta$config_id, v$id, sep = "::"),
    beta_reference = meta$beta_reference, beta = NA_real_, standard_error = NA_real_,
    beta_error = NA_real_, fit_criterion_met = NA, software_converged = FALSE,
    full_iteration_limit = NA, iterations = NA_real_, global_null_iterations = NA_real_,
    ll_change = NA_real_, max_abs_score = NA_real_, beta_change = NA_real_,
    ci_lower_iterations = NA_real_, ci_upper_iterations = NA_real_, plrt_iterations = NA_real_,
    ci_lower_iteration_limit = NA, ci_upper_iteration_limit = NA, plrt_iteration_limit = NA,
    ci_lower_ll_change = NA_real_, ci_lower_beta_change = NA_real_,
    ci_upper_ll_change = NA_real_, ci_upper_beta_change = NA_real_,
    ci_criterion_met = NA, ci_converged = NA,
    pl_low = NA_real_, pl_high = NA_real_, pl_p = NA_real_,
    wald_low = NA_real_, wald_high = NA_real_, wald_p = NA_real_,
    fit_status = if (cap$error) "error" else "pending",
    effective_control_matches = NA, warning_count = length(cap$warnings),
    warning_or_error = collapse_messages(c(cap$warnings, if (cap$error) cap$value$message)),
    elapsed_seconds = unname(elapsed[["elapsed"]]), cpu_user_seconds = unname(elapsed[["user.self"]]),
    cpu_system_seconds = unname(elapsed[["sys.self"]]),
    software_version = as.character(packageVersion(v$implementation)), stringsAsFactors = FALSE
  )
  if (!cap$error) {
    fit <- cap$value
    k <- match("exposure", names(coef(fit)))
    audit$beta <- unname(coef(fit)[k])
    if (v$implementation == "logistf") {
      keys <- setdiff(names(v$control), "call")
      audit$effective_control_matches <- isTRUE(all.equal(fit$control[keys], v$control[keys], check.attributes = FALSE))
      audit$standard_error <- sqrt(fit$var[k, k])
      cv <- as.numeric(fit$conv)
      audit[c("ll_change", "max_abs_score", "beta_change")] <- as.list(cv)
      audit$fit_criterion_met <- finite_all(cv) && all(abs(cv) <= unlist(v$control[c("lconv", "gconv", "xconv")]))
      audit$iterations <- unname(fit$iter["full"])
      audit$global_null_iterations <- unname(fit$iter["null"])
      audit$full_iteration_limit <- audit$iterations >= v$control$maxit
      audit$software_converged <- audit$fit_criterion_met && !audit$full_iteration_limit
      # pl.conv rows follow plconf, while pl.iter rows follow the full model.
      pr <- match(as.character(k), rownames(fit$pl.conv))
      pc <- as.numeric(fit$pl.conv[pr, ])
      audit[c("ci_lower_ll_change", "ci_lower_beta_change", "ci_upper_ll_change", "ci_upper_beta_change")] <- as.list(pc)
      audit$ci_criterion_met <- finite_all(pc) && all(abs(pc) <= unlist(v$plcontrol[c("lconv", "xconv", "lconv", "xconv")]))
      pi <- as.numeric(fit$pl.iter[k, ])
      audit[c("ci_lower_iterations", "ci_upper_iterations", "plrt_iterations")] <- as.list(pi)
      audit$ci_lower_iteration_limit <- pi[1] >= v$plcontrol$maxit
      audit$ci_upper_iteration_limit <- pi[2] >= v$plcontrol$maxit
      audit$plrt_iteration_limit <- pi[3] >= v$control$maxit
      audit$ci_converged <- audit$software_converged && audit$ci_criterion_met &&
        !audit$ci_lower_iteration_limit && !audit$ci_upper_iteration_limit
      audit$pl_low <- unname(fit$ci.lower[k])
      audit$pl_high <- unname(fit$ci.upper[k])
      audit$pl_p <- unname(fit$prob[k])
    } else {
      keys <- setdiff(names(v$control), "call")
      audit$effective_control_matches <- identical(fit$type, v$control$type) &&
        isTRUE(all.equal(fit$control[keys], v$control[keys], check.attributes = FALSE))
      audit$standard_error <- sqrt(vcov(fit)[k, k])
      audit$software_converged <- isTRUE(fit$converged)
      audit$iterations <- fit$iter
      audit$full_iteration_limit <- fit$iter >= v$control$maxit
      audit$ci_converged <- audit$software_converged
    }
    stopifnot(audit$effective_control_matches)
    audit$wald_low <- audit$beta - qnorm(.975) * audit$standard_error
    audit$wald_high <- audit$beta + qnorm(.975) * audit$standard_error
    audit$wald_p <- 2 * pnorm(-abs(audit$beta / audit$standard_error))
    audit$beta_error <- abs(audit$beta - audit$beta_reference)
    audit$fit_status <- if (!is.finite(audit$beta)) "estimate_nonfinite"
      else if (isTRUE(audit$full_iteration_limit)) "iteration_limit"
      else if (!isTRUE(audit$software_converged)) "criterion_not_met"
      else "converged"
  }
  outputs <- if (v$implementation == "logistf") c("PL", "Wald") else "Wald"
  result <- do.call(rbind, lapply(outputs, function(inference) {
    z <- audit[c("config_id", "scenario", "n", "n_exposed", "n_unexposed", "a", "b", "c", "d",
      "separation", "sparsity_class", "variant_id", "implementation", "algorithm", "initial_mode", "fit_id",
      "beta", "standard_error", "beta_reference", "beta_error", "software_converged", "full_iteration_limit",
      "iterations", "fit_status", "warning_count", "warning_or_error", "elapsed_seconds", "software_version")]
    z$inference <- inference
    z$ci_low <- if (inference == "PL") audit$pl_low else audit$wald_low
    z$ci_high <- if (inference == "PL") audit$pl_high else audit$wald_high
    z$p_value <- if (inference == "PL") audit$pl_p else audit$wald_p
    z$ci_converged <- if (inference == "PL") audit$ci_converged else audit$software_converged
    z$plrt_iteration_limit <- if (inference == "PL") audit$plrt_iteration_limit else NA
    z$ci_finite <- is.finite(z$ci_low) && is.finite(z$ci_high)
    z$ci_excludes_zero <- if (z$ci_finite) z$ci_low > 0 || z$ci_high < 0 else NA
    z$p_below_005 <- if (is.finite(z$p_value)) z$p_value < .05 else NA
    z
  }))
  list(audit = audit, result = result, fit = if (cap$error) NULL else cap$value)
}

