# Coverage of the exposure coefficient in the saved covariate simulations.
# No input generation or model fitting. Run from the project root.
source("R/project_helpers.R", local = TRUE)
set_project_locale()
source("R/analysis_helpers.R", local = TRUE)
ensure_directories(c("results/tables", "results/figures", "results/checks"))
design <- read.csv("data/covariate_scenarios.csv", stringsAsFactors = FALSE)
meta <- readRDS("data/covariate_metadata.rds")
fits <- readRDS("results/model_fits.rds")
fits <- fits[fits$config_id %in% meta$config_id, ]
stopifnot(nrow(design) == 8L, nrow(meta) == 20000L,
          !anyDuplicated(meta$config_id), all(design$OR == 2),
          all(design$repetitions == 2500L))
methods <- sort(unique(fits$key))
eligible <- meta$config_id[meta$data_state == "ok"]
stopifnot(length(methods) == 36L, nrow(fits) == length(eligible) * 36L,
          !anyDuplicated(paste(fits$config_id, fits$key)))
for (method in methods)
  stopifnot(setequal(fits$config_id[fits$key == method], eligible))
fits$scenario_id <- meta$design_id[match(fits$config_id, meta$config_id)]
fits$truth <- log(design$OR[match(fits$scenario_id, design$scenario_id)])
fits$available <- ci_ok(fits)
fits$covered <- fits$available & fits$low <= fits$truth & fits$truth <= fits$high
fits$point_ok <- confirmed_convergence(fits$converged)
fits$ci_ok <- confirmed_convergence(fits$ci_converged)

summarize_group <- function(z) {
  s <- design[match(z$scenario_id[1], design$scenario_id), ]
  generated <- meta[meta$design_id == s$scenario_id, ]
  stopifnot(nrow(generated) == s$repetitions,
            nrow(z) == sum(generated$data_state == "ok"))
  cohorts <- list(available = z$available,
                  point_converged = z$available & z$point_ok,
                  point_and_ci_converged = z$available & z$point_ok & z$ci_ok)
  do.call(rbind, lapply(names(cohorts), function(cohort) {
    keep <- cohorts[[cohort]]
    denominator <- sum(keep)
    covered <- sum(z$covered & keep)
    # The two noncoverage tails partition the valid intervals with coverage.
    below <- sum(keep & z$high < z$truth)
    above <- sum(keep & z$low > z$truth)
    stopifnot(covered + below + above == denominator)
    estimate <- if (denominator) covered / denominator else NA_real_
    mcse <- if (denominator) sqrt(estimate * (1 - estimate) / denominator) else NA_real_
    data.frame(scenario_id = s$scenario_id, n = s$n, n_exposed = s$n_exposed,
      rho = s$rho, true_beta = log(s$OR), method = z$key[1],
      engine = z$engine[1], variant = z$variant[1], inference = z$inference[1],
      cohort = cohort, generated_n = nrow(generated),
      excluded_one_outcome_n = sum(generated$data_state != "ok"), fitted_n = nrow(z),
      ci_available_n = sum(z$available), ci_unavailable_fitted_n = sum(!z$available),
      ci_available_fraction = sum(z$available) / nrow(generated),
      point_nonconverged_available_n = sum(z$available & !is.na(z$converged) & !z$converged),
      point_unknown_available_n = sum(z$available & is.na(z$converged)),
      ci_nonconverged_available_n = sum(z$available & !is.na(z$ci_converged) & !z$ci_converged),
      ci_unknown_available_n = sum(z$available & is.na(z$ci_converged)),
      denominator_n = denominator, covered_n = covered,
      interval_below_truth_n = below, interval_above_truth_n = above,
      coverage = estimate, mcse = mcse,
      mc_low = max(0, estimate - 1.96 * mcse),
      mc_high = min(1, estimate + 1.96 * mcse))
  }))
}
groups <- split(fits, interaction(fits$scenario_id, fits$key, drop = TRUE))
coverage <- do.call(rbind, lapply(groups, summarize_group))
rownames(coverage) <- NULL
coverage <- coverage[order(coverage$n, coverage$rho, coverage$method, coverage$cohort), ]
stopifnot(nrow(coverage) == 8L * 36L * 3L)

# Independently aggregate OR-scale inclusion across valid individual outputs.
valid <- fits[fits$available, ]
valid$covered_or <- as.integer(exp(valid$low) <= exp(valid$truth) &
                                exp(valid$truth) <= exp(valid$high))
valid$one <- 1L
check <- aggregate(cbind(covered_or, one) ~ scenario_id + key, valid, sum)
observed <- coverage[coverage$cohort == "available", ]
idx <- match(paste(observed$scenario_id, observed$method), paste(check$scenario_id, check$key))
stopifnot(!anyNA(idx), all(observed$covered_n == check$covered_or[idx]),
          all(observed$denominator_n == check$one[idx]),
          all(coverage$excluded_one_outcome_n + coverage$fitted_n == coverage$generated_n),
          all(coverage$ci_available_n + coverage$ci_unavailable_fitted_n == coverage$fitted_n),
          all(coverage$denominator_n <= coverage$ci_available_n),
          all(is.na(coverage$coverage[coverage$denominator_n == 0])))
# Confirm denominators against the existing independently generated status table.
status <- read.csv("results/tables/method_status.csv", stringsAsFactors = FALSE)
si <- match(paste(observed$scenario_id, observed$method), paste(status$scenario_id, status$method))
stopifnot(!anyNA(si),
  all(abs(observed$ci_available_fraction -
    (1 - status$excluded_probability[si] - status$ci_missing_probability[si])) < 1e-12),
  all(coverage$coverage[coverage$denominator_n > 0] >= 0),
  all(coverage$coverage[coverage$denominator_n > 0] <= 1))
write.csv(coverage, "results/tables/covariate_coverage.csv", row.names = FALSE, na = "")

# Default fits, including the alternative inference output from the same fit.
selected <- c("sas_proc|fisher_default|Wald", "sas_proc|fisher_default|PL",
              "sas_fl|default|PL", "logistf|lf_nr_default|Wald",
              "logistf|lf_nr_default|PL", "brglm2|bg_asmean_default|Wald")
defaults <- observed[observed$method %in% selected, ]
defaults <- defaults[order(defaults$n, defaults$rho, match(defaults$method, selected)), ]
stopifnot(nrow(defaults) == 48L)
write.csv(defaults, "results/tables/covariate_coverage_default.csv", row.names = FALSE, na = "")

library(ggplot2)
plot_data <- defaults
plot_data$implementation <- factor(plot_data$engine,
  c("sas_proc", "sas_fl", "logistf", "brglm2"),
  c("PROC LOGISTIC", "%FL", "logistf", "brglm2"))
plot_data$background <- factor(ifelse(plot_data$n == 40, "40例・20:20", "100例・20:80"),
                               c("40例・20:20", "100例・20:80"))
plot_data$inference <- factor(plot_data$inference, c("Wald", "PL"), c("Wald区間", "PL区間"))
plot_data$x <- match(plot_data$rho, c(0, .5, .9, .99))
dodge <- position_dodge(width = .18)
p <- ggplot(plot_data, aes(x, 100 * coverage, colour = implementation, group = implementation)) +
  geom_hline(yintercept = 95, linetype = "dashed", colour = "#59676A") +
  geom_line(linewidth = .8, position = dodge) +
  geom_errorbar(aes(ymin = 100 * mc_low, ymax = 100 * mc_high,
    linetype = "エラーバー：Monte Carlo誤差に基づく95%信頼区間"),
    width = .07, position = dodge, show.legend = c(colour = FALSE, linetype = TRUE)) +
  geom_point(aes(shape = implementation), size = 2.6, position = dodge) +
  facet_grid(inference ~ background) +
  scale_x_continuous(breaks = 1:4, labels = c("0", "0.5", "0.9", "0.99")) +
  scale_colour_manual(values = c("#23373B", "#8A5177", "#E87722", "#007C91")) +
  scale_linetype_manual(values = 1) +
  guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, override.aes = list(colour = "#59676A"))) +
  labs(title = "共変量モデル：95%信頼区間のカバレッジ",
    subtitle = "既定の当てはめ設定／真の曝露効果 log(2)／破線：95%",
    x = "連続共変量間の相関", y = "区間取得例におけるカバレッジ（%）",
    colour = NULL, shape = NULL, linetype = NULL,
    caption = "有限で順序が正しい区間を分母とし、未収束でも出力があれば含む。\nアウトカム一水準・区間欠測は除外。区間取得率・収束状態は集計表に併記。") +
  theme_minimal(base_size = 16, base_family = "Yu Gothic") +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = grid::unit(2, "pt"),
        plot.caption = element_text(hjust = 0, size = 12),
        plot.margin = margin(12, 16, 12, 12))
ggsave("results/figures/covariate_coverage.png", p, width = 11, height = 8, dpi = 180, bg = "white")
writeLines(c("estimand=conditional exposure log odds ratio log(2), not the null value 0",
  "scope=8 covariate scenarios x 36 output settings x 3 availability/convergence cohorts",
  "coverage=P(interval includes log(2) | interval available and cohort selected)",
  "availability=finite endpoints in correct order; endpoints inclusive",
  "default_figure=6 outputs; default fits with Wald and PL output shown separately",
  "unknown_convergence=not treated as converged; empty cohorts yield NA, never zero coverage",
  "MCSE=sqrt(coverage*(1-coverage)/denominator); error bars=+/-1.96*MCSE, clipped to [0,1]",
  "checks=complete unique eligible IDs; count partitions; independent OR-scale inclusion; existing method-status denominators",
  "no_model_refits=TRUE"), "results/checks/coverage.txt")
print(defaults[c("scenario_id", "method", "ci_available_n", "coverage", "mcse")], row.names = FALSE)
