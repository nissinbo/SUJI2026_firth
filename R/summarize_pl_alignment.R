# Inspect the same tables before and after aligning PROC LOGISTIC to PL.
# Selection is based only on the original Wald-versus-PL decision difference.
source("R/analysis_helpers.R", local = TRUE)
tables <- read.csv("data/tables.csv")
tables <- tables[tables$design_id == "E050R1" & tables$data_state == "ok", ]
fits <- readRDS("results/model_fits.rds")
aligned_fit <- function(key) {
  result <- fits[fits$key == key, ]
  result <- result[match(tables$config_id, result$config_id), ]
  stopifnot(identical(result$config_id, tables$config_id), all(ci_ok(result)))
  result
}
sas_wald <- aligned_fit("sas_proc|fisher_default|Wald")
sas_pl <- aligned_fit("sas_proc|fisher_default|PL")
logistf_pl <- aligned_fit("logistf|lf_nr_default|PL")
excludes_zero <- function(result) result$low > 0 | result$high < 0
ci_difference <- function(a, b) pmax(abs(a$low - b$low), abs(a$high - b$high))
different <- xor(excludes_zero(sas_wald), excludes_zero(logistf_pl))
after_different <- xor(excludes_zero(sas_pl), excludes_zero(logistf_pl))
stopifnot(nrow(tables) == 674L, sum(different) == 36L, !any(after_different))

details <- data.frame(
  config_id = tables$config_id, exposed_events = tables$a, unexposed_events = tables$c,
  before_different = different, after_different = after_different,
  sas_wald_low = sas_wald$low, sas_wald_high = sas_wald$high,
  sas_pl_low = sas_pl$low, sas_pl_high = sas_pl$high,
  logistf_pl_low = logistf_pl$low, logistf_pl_high = logistf_pl$high,
  before_delta = ci_difference(sas_wald, logistf_pl),
  after_delta = ci_difference(sas_pl, logistf_pl)
)
summary <- do.call(rbind, lapply(c("original_disagreements", "all_available"), function(scope) {
  selected <- if (scope == "original_disagreements") different else rep(TRUE, nrow(tables))
  do.call(rbind, lapply(c("before", "after"), function(stage) {
    values <- details[[paste0(stage, "_delta")]][selected]
    data.frame(scope, stage, tables = sum(selected),
      disagreements = sum(details[[paste0(stage, "_different")]][selected]),
      median = median(values), maximum = max(values))
  }))
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(details, "results/tables/pl_alignment_details.csv", row.names = FALSE)
write.csv(summary, "results/tables/pl_alignment_summary.csv", row.names = FALSE)

# The illustrated 6-versus-1 table is close to the median original difference.
example_id <- "n50_6.19.1.24"
i <- match(example_id, tables$config_id)
stopifnot(different[i], all(c(sas_wald$converged[i], sas_pl$converged[i], logistf_pl$converged[i])),
          abs(details$before_delta[i] - median(details$before_delta[different])) < 1e-7)
example <- data.frame(
  config_id = example_id,
  method = c("PROC LOGISTIC (Wald)", "logistf (PL)", "PROC LOGISTIC (PL)"),
  low = c(sas_wald$low[i], logistf_pl$low[i], sas_pl$low[i]),
  high = c(sas_wald$high[i], logistf_pl$high[i], sas_pl$high[i])
)
example$includes_zero <- !excludes_zero(example)
write.csv(example, "results/tables/pl_alignment_example.csv", row.names = FALSE)
print(summary, row.names = FALSE, digits = 8)
print(example, row.names = FALSE, digits = 9)
