# Reuse saved fits; no models or simulation scenarios are rerun.
Sys.setlocale("LC_CTYPE", "Japanese_Japan.utf8")
source("R/analysis_helpers.R", local = TRUE)
library(ggplot2)
tables <- read.csv("data/tables.csv", stringsAsFactors = FALSE)
z <- tables[tables$design_id == "E050R1", ]
fits <- readRDS("results/model_fits.rds")
get_fit <- function(key) {
  out <- fits[fits$key == key & fits$config_id %in% z$config_id, ]
  stopifnot(!anyDuplicated(out$config_id), nrow(out) == 674L)
  out[match(z$config_id, out$config_id), ]
}
sas <- get_fit("sas_proc|fisher_default|Wald")
sas_pl <- get_fit("sas_proc|fisher_default|PL")
lf <- get_fit("logistf|lf_nr_default|PL")
lf_wald <- get_fit("logistf|lf_nr_default|Wald")
bg <- get_fit("brglm2|bg_asmean_default|Wald")
z$ci_available <- ci_ok(sas) & ci_ok(lf)
z$ci_disagree <- z$ci_available & xor(sas$low > 0 | sas$high < 0,
                                     lf$low > 0 | lf$high < 0)
z$ci_delta <- ifelse(z$ci_available,
  pmax(abs(sas$low - lf$low), abs(sas$high - lf$high)), NA_real_)
z$sas_point_converged <- sas$converged
z$logistf_point_converged <- lf$converged
z$sas_low <- sas$low; z$sas_high <- sas$high
z$logistf_low <- lf$low; z$logistf_high <- lf$high
z$aligned_ci_available <- ci_ok(sas_pl) & ci_ok(lf)
z$aligned_ci_disagree <- z$aligned_ci_available & xor(
  sas_pl$low > 0 | sas_pl$high < 0, lf$low > 0 | lf$high < 0)
z$aligned_decision <- ifelse(z$data_state != "ok", "全例イベントあり/なし",
  ifelse(!z$aligned_ci_available, "信頼区間取得不可",
    ifelse(z$aligned_ci_disagree, "判定差あり", "判定一致")))
z$aligned_decision <- factor(z$aligned_decision,
  c("判定一致", "判定差あり", "信頼区間取得不可", "全例イベントあり/なし"))
z$weight <- rep(1/nrow(z), nrow(z))
z$decision <- ifelse(z$data_state != "ok", "全例イベントあり/なし",
  ifelse(!z$ci_available, "信頼区間取得不可", ifelse(z$ci_disagree, "判定差あり", "判定一致")))
z$decision <- factor(z$decision, c("判定一致", "判定差あり", "信頼区間取得不可", "全例イベントあり/なし"))
available <- sum(z$weight[z$ci_available])
joint <- sum(z$weight[z$ci_disagree])
conditional <- joint / available
summary <- read.csv("results/tables/pair_summary.csv")
ref <- subset(summary, scenario_id == "E050R1" &
  comparison == "default_logistf_sas" & cohort == "all_outputs")
stopifnot(nrow(z) == 676L, !anyDuplicated(paste(z$a, z$c)),
  all(z$a + z$b == 25), all(z$c + z$d == 25),
  abs(sum(z$weight) - 1) < 1e-12, nrow(ref) == 1L,
  abs(available - ref$ci_available_probability) < 1e-12,
  abs(conditional - ref$ci_difference_conditional_probability) < 1e-12)
aligned_ref <- subset(summary, scenario_id == "E050R1" &
  comparison == "aligned_PL_logistf_sas" & cohort == "all_outputs")
stopifnot(nrow(aligned_ref) == 1L, sum(z$ci_disagree) == 36L,
  sum(z$aligned_ci_available) == 674L, sum(z$aligned_ci_disagree) == 0L,
  abs(sum(z$weight[z$aligned_ci_available]) - aligned_ref$ci_available_probability) < 1e-12,
  aligned_ref$ci_difference_conditional_probability == 0)

# Compare Wald intervals on the same 676 tables as the PL comparison.
wald_decision <- function(a, b) {
  available <- ci_ok(a) & ci_ok(b)
  different <- available & xor(a$low > 0 | a$high < 0, b$low > 0 | b$high < 0)
  factor(ifelse(z$data_state != "ok", "全例イベントあり/なし",
    ifelse(!available, "信頼区間取得不可", ifelse(different, "判定差あり", "判定一致"))),
    levels = levels(z$decision))
}
z$wald_logistf_decision <- wald_decision(sas, lf_wald)
z$wald_brglm2_decision <- wald_decision(sas, bg)
stopifnot(sum(z$wald_logistf_decision == "判定差あり") == 12L,
          sum(z$wald_brglm2_decision == "判定差あり") == 0L,
          identical(z$wald_logistf_decision, wald_decision(bg, lf_wald)))
wald_ref <- subset(summary, scenario_id == "E050R1" &
  comparison == "aligned_Wald_logistf_sas" & cohort == "all_outputs")
wald_probability <- sum(z$weight[z$wald_logistf_decision == "判定差あり"]) /
  sum(z$weight[ci_ok(sas) & ci_ok(lf_wald)])
stopifnot(nrow(wald_ref) == 1L,
          abs(wald_probability - wald_ref$ci_difference_conditional_probability) < 1e-12)

# Recalculate a common ordinary-information SE using the closed-form Firth
# group probabilities for this exposure-only 2x2 model. Keep each saved beta.
valid <- z$data_state == "ok"
p1 <- (z$a + 0.5) / (z$a + z$b + 1)
p0 <- (z$c + 0.5) / (z$c + z$d + 1)
z$wald_common_se <- ifelse(valid,
  sqrt(1 / (25 * p1 * (1 - p1)) + 1 / (25 * p0 * (1 - p0))), NA_real_)
common_sas <- data.frame(low = sas$beta - qnorm(.975) * z$wald_common_se,
                         high = sas$beta + qnorm(.975) * z$wald_common_se)
common_lf <- data.frame(low = lf_wald$beta - qnorm(.975) * z$wald_common_se,
                        high = lf_wald$beta + qnorm(.975) * z$wald_common_se)
z$wald_common_sas_low <- common_sas$low
z$wald_common_sas_high <- common_sas$high
z$wald_common_logistf_low <- common_lf$low
z$wald_common_logistf_high <- common_lf$high
z$wald_common_se_decision <- wald_decision(common_sas, common_lf)
se_reference <- read.csv("results/tables/wald_se.csv")
se_reference <- se_reference[match(z$config_id[valid], se_reference$config_id), ]
stopifnot(sum(ci_ok(common_sas) & ci_ok(common_lf)) == 674L,
          sum(z$wald_common_se_decision == "判定差あり") == 0L,
          sum(z$wald_common_se_decision == "判定一致") == 674L,
          sum(z$wald_common_se_decision == "全例イベントあり/なし") == 2L,
          max(abs(z$wald_common_se[valid] - se_reference$ordinary_se)) < 1e-12,
          identical(common_sas$low[valid] > 0 | common_sas$high[valid] < 0,
                    sas$low[valid] > 0 | sas$high[valid] < 0))
write.csv(z, "results/tables/exact_heatmap_n50.csv", row.names = FALSE)
theme_set(theme_minimal(base_size = 17, base_family = "Yu Gothic") +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "#23373B"),
        plot.title = element_text(face = "bold", size = 18),
        legend.position = "bottom", legend.title = element_blank(),
        plot.margin = margin(8, 10, 5, 8)))
base <- ggplot(z, aes(a, c)) +
  scale_x_continuous(breaks = seq(0, 25, 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 25, 5), expand = c(0, 0)) +
  coord_fixed() + labs(x = "曝露群のイベント数", y = "非曝露群のイベント数")
decision <- base + geom_tile(aes(fill = decision)) +
  scale_fill_manual(values = c("判定一致" = "#E3EBED", "判定差あり" = "#E87722",
    "信頼区間取得不可" = "#743D7B", "全例イベントあり/なし" = "#23373B"), drop = TRUE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "PROC LOGISTIC（Wald）\n対 logistf（PL）")
aligned <- base + geom_tile(aes(fill = aligned_decision), show.legend = TRUE) +
  scale_fill_manual(values = c("判定一致" = "#E3EBED", "判定差あり" = "#E87722",
    "信頼区間取得不可" = "#743D7B", "全例イベントあり/なし" = "#23373B"),
    breaks = c("判定一致", "判定差あり", "全例イベントあり/なし"), drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "PROC LOGISTIC（PL）\n対 logistf（PL）")
wald_plot <- function(column, other) {
  base + geom_tile(aes(fill = .data[[column]]), show.legend = TRUE) +
    scale_fill_manual(values = c("判定一致" = "#E3EBED", "判定差あり" = "#E87722",
      "信頼区間取得不可" = "#743D7B", "全例イベントあり/なし" = "#23373B"),
      breaks = c("判定一致", "判定差あり", "全例イベントあり/なし"), drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = paste0("PROC LOGISTIC（Wald）\n対 ", other, "（Wald）"))
}
wald_logistf <- wald_plot("wald_logistf_decision", "logistf")
wald_common_se <- wald_plot("wald_common_se_decision", "logistf") +
  labs(title = "PROC LOGISTIC 対 logistf\n標準誤差を揃えて再計算")
# Measure Japanese text on the same PNG device used for the final figures.
# The default PDF device lacks the Windows font database.
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(tempfile(fileext = ".png"), width = 5.5, height = 5.25,
                units = "in", res = 190)
} else {
  grDevices::png(tempfile(fileext = ".png"), width = 5.5, height = 5.25,
                 units = "in", res = 190, type = "cairo-png")
}
plot_grobs <- lapply(list(decision, aligned, wald_logistf, wald_common_se), ggplotGrob)
aligned_widths <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$widths))
aligned_heights <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$heights))
plot_grobs <- lapply(plot_grobs, function(g) {
  g$widths <- aligned_widths; g$heights <- aligned_heights; g
})
invisible(grDevices::dev.off())
ggsave("results/figures/exact_heatmap_decision.png", plot_grobs[[1]],
       width = 5.5, height = 5.25, dpi = 190, bg = "white")
ggsave("results/figures/exact_heatmap_aligned.png", plot_grobs[[2]],
       width = 5.5, height = 5.25, dpi = 190, bg = "white")
ggsave("results/figures/exact_heatmap_wald_logistf.png", plot_grobs[[3]],
       width = 5.5, height = 5.25, dpi = 190, bg = "white")
ggsave("results/figures/exact_heatmap_wald_common_se.png", plot_grobs[[4]],
       width = 5.5, height = 5.25, dpi = 190, bg = "white")
record <- c("Source: results/model_fits.rds; data/tables.csv",
  "Fixed group sizes: 25 exposed and 25 unexposed; a,c each range from 0 to 25",
  "Default settings; nonconverged numeric outputs retained as in pair_summary.csv",
  sprintf("Configurations: %d; available confidence interval pairs: %d; disagreements: %d",
    nrow(z), sum(z$ci_available), sum(z$ci_disagree)),
  sprintf("Weight sum: %.17g", sum(z$weight)),
  sprintf("Available-pair fraction of all configurations: %.17g", available),
  sprintf("Disagreement fraction of all configurations: %.17g", joint),
  sprintf("Disagreement fraction of available pairs: %.17g", conditional),
  "Aggregation: each table counted once; no event probability or true OR specified",
  sprintf("PL vs PL: available confidence interval pairs %d; disagreements %d",
    sum(z$aligned_ci_available), sum(z$aligned_ci_disagree)),
  sprintf("Wald: PROC vs logistf 12/674; PROC vs brglm2 0/674; table fraction %.17g", wald_probability),
  "Wald: brglm2 vs logistf has the same disagreement cells as PROC vs logistf.",
  "Recalculated Wald: each saved PROC/logistf beta retained; common ordinary-information SE at closed-form Firth group probabilities.",
  "Recalculated Wald: available pairs 674; disagreements 0; excluded all-event/no-event tables 2.",
  "Common SE agrees with wald_se.csv to 1e-12; recalculated decisions agree with saved PROC decisions.",
  "Validated against saved pair_summary.csv to 1e-12.")
writeLines(record, "results/checks/heatmaps.txt")
cat(paste(record, collapse = "\n"), "\n")

# Keep the companion Wald comparison table on the same equal-count basis.
wald_pairs <- list(c("sas_proc|fisher_default|Wald", "logistf|lf_nr_default|Wald"),
                   c("sas_proc|fisher_default|Wald", "brglm2|bg_asmean_default|Wald"),
                   c("logistf|lf_nr_default|Wald", "brglm2|bg_asmean_default|Wald"))
wald_counts <- do.call(rbind, lapply(wald_pairs, function(keys) {
 a <- get_fit(keys[1]); b <- get_fit(keys[2]); valid <- ci_ok(a) & ci_ok(b)
 different <- xor(a$low[valid] > 0 | a$high[valid] < 0,
                  b$low[valid] > 0 | b$high[valid] < 0)
 data.frame(scenario="E050R1", A=keys[1], B=keys[2], available=sum(valid),
            different=sum(different), proportion=mean(different))
}))
stopifnot(identical(wald_counts$different, c(12L, 0L, 12L)), all(wald_counts$available==674L))
write.csv(wald_counts, "results/tables/wald_aligned_examples.csv", row.names=FALSE)
