# Summarize and plot saved introduction fits; no simulation is rerun.
source("R/project_helpers.R", local = TRUE)
set_project_locale()
intro <- readRDS("results/raw/intro_bias.rds")
intro_n <- intro$n
intro_B <- intro$repetitions
intro_truth <- intro$truth
intro_estimates <- intro$estimates
intro_separated <- intro$separated
intro_used <- sum(!intro_separated)
intro_summary <- do.call(rbind, lapply(c("MLE", "Firth"), function(method) {
  estimates <- intro_estimates[!intro_separated, , method, drop = FALSE]
  dim(estimates) <- c(intro_used, 3L)
  data.frame(method = method, term = c("x1", "x2", "x3"),
    truth = intro_truth[-1], mean = colMeans(estimates),
    bias = colMeans(estimates) - intro_truth[-1],
    mcse = apply(estimates, 2, sd) / sqrt(intro_used),
    n = intro_n, repetitions = intro_B, used = intro_used,
    separated = sum(intro_separated), seed = intro$seed)
}))
stopifnot(intro_used > 0L, all(is.finite(intro_summary$mcse)))

write.csv(intro_summary, "results/tables/intro_bias.csv", row.names = FALSE)
library(ggplot2)
# Display the first coefficient; retain all three in the saved summary.
intro_plot <- subset(intro_summary, term == "x1")
intro_plot$method <- factor(intro_plot$method, c("MLE", "Firth"),
  c("通常のロジスティック回帰", "Firth法"))
intro_figure <- ggplot(intro_plot, aes(method, bias, colour = method, shape = method)) +
  geom_hline(yintercept = 0, colour = "#23373B", linetype = "dashed") +
  geom_errorbar(aes(ymin = bias - 1.96 * mcse, ymax = bias + 1.96 * mcse,
    linetype = "エラーバー：Monte Carlo誤差に基づく95%信頼区間"),
    width = .10, show.legend = c(linetype = TRUE, colour = FALSE)) +
  geom_point(size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = c("#E87722", "#007C91")) +
  scale_shape_manual(values = c(16, 17)) + scale_linetype_manual(values = 1) +
  guides(colour = "none", shape = "none",
    linetype = guide_legend(order = 2, override.aes = list(colour = "#23373B", shape = NA))) +
  labs(x = NULL, y = "推定値の平均\n− 真値", colour = NULL, shape = NULL, linetype = NULL) +
  theme_minimal(base_size = 17, base_family = "Yu Gothic") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    axis.text = element_text(colour = "#23373B"), legend.position = "bottom",
    legend.box = "vertical", legend.spacing.y = grid::unit(0, "pt"),
    plot.margin = margin(5, 15, 5, 10))
ggsave("results/figures/intro_bias.png", intro_figure,
  width = 10, height = 3.4, dpi = 180, device = ragg::agg_png, bg = "white")
