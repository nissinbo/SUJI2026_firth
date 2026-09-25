source("R/project_helpers.R", local = TRUE)
source("R/analysis_helpers.R", local = TRUE)
set_project_locale()

folder<-"results/tables"
z<-read.csv(file.path(folder,"pair_summary.csv"))
defs<-read.csv(file.path(folder,"contrast_definitions.csv"))
stopifnot(nrow(z)==nrow(defs)*16L*4L,!anyDuplicated(paste(z$scenario_id,z$comparison,z$cohort)))
for(type in c("ci","p")) {
 available<-z[[paste0(type,"_available_probability")]]
 joint<-z[[paste0(type,"_difference_joint_probability")]]
 conditional<-z[[paste0(type,"_difference_conditional_probability")]]
 stopifnot(all(available>=0&available<=z$cohort_probability+1e-12),all(joint>=0&joint<=available+1e-12),
  all(abs(conditional[available>0]-joint[available>0]/available[available>0])<1e-12),all(is.na(conditional[available==0])))
 counts<-z[[paste0(type,"_available_n")]]
 differences<-z[[paste0(type,"_difference_n")]]
 stopifnot(all(abs(conditional[counts>0]-differences[counts>0]/counts[counts>0])<1e-12),
  all(abs(available-counts/z$total_configurations)<1e-12))
}
stopifnot(all(z$p_delta_max[z$comparison=="SAS_CI"]==0),all(z$beta_delta_p99[z$comparison=="LF_inference"]==0))
scenarios<-read.csv("data/scenarios.csv")
diagnostics<-read.csv(file.path(folder,"scenario_diagnostics.csv"))
d<-merge(scenarios,diagnostics,by="scenario_id")
prob<-with(d,2/((n_exposed+1)*(n_unexposed+1)))
stopifnot(max(abs(prob-d$exclusion_probability))<1e-12)
mc<-z$kind=="simulation"&z$ci_available_n>0
expected_mcse<-with(z[mc,],sqrt(ci_difference_conditional_probability*(1-ci_difference_conditional_probability)/ci_available_n))
stopifnot(max(abs(z$ci_difference_mcse[mc]-expected_mcse))<1e-12)
writeLines(c(paste0("verified=",format(Sys.time())),"scenario_comparison_cohort_keys=complete and unique",
 "probability_accounting=joint <= available <= cohort; conditional=joint/available",
 "equal_count_proportions=disagreements / available pairs for every comparison and cohort",
 "one_outcome_exclusion_proportion=2 / ((n_exposed+1)*(n_unexposed+1))",
 "same_SAS_fit_Wald_vs_PL=p-values identical", "same_logistf_fit_Wald_vs_PL=coefficients identical",
 "simulation_MCSE=matches available-pair binomial formula"),"results/checks/summary.txt")
cat("Independent summary checks passed.\n")

# Independently check equal-count summaries against saved per-table fits.

tabdir <- "results/tables"
grid <- read.csv("data/tables.csv", stringsAsFactors=FALSE)
fits <- readRDS("results/model_fits.rds")
fits <- fits[fits$config_id %in% grid$config_id, ]
methods <- split(fits, fits$key)
rm(fits)
aligned <- lapply(methods, function(z) z[match(grid$config_id,z$config_id), ])
known_methods <- aligned[!grepl("^sas_fl", names(aligned))]
common <- Reduce(`&`, lapply(known_methods, function(z) !is.na(z$converged) & z$converged))
common_inference <- common & Reduce(`&`, lapply(known_methods, function(z) ci_ok(z) & p_ok(z)))
summary <- read.csv(file.path(tabdir,"pair_summary.csv"))
defs <- read.csv(file.path(tabdir,"contrast_definitions.csv"))
checked <- 0L
equal <- function(actual, expected) stopifnot(isTRUE(all.equal(actual, expected, tolerance=1e-12, check.attributes=FALSE)))
quantiles <- function(x) {
 x <- x[is.finite(x)]
 if(!length(x)) return(rep(NA_real_,3))
 c(quantile(x,c(.5,.99),type=7,names=FALSE),max(x))
}
for (i in seq_len(nrow(defs))) {
 a <- aligned[[defs$method_A[i]]]; b <- aligned[[defs$method_B[i]]]
 paired <- !is.na(a$converged) & a$converged & !is.na(b$converged) & b$converged
 for (design in unique(grid$design_id)) for (cohort in unique(summary$cohort)) {
  mask <- grid$design_id==design & switch(cohort, all_outputs=rep(TRUE,nrow(grid)),
    both_point_converged=paired,common_point_converged=common,common_inference_available=common_inference)
  row <- summary[summary$scenario_id==design & summary$comparison==defs$comparison[i] & summary$cohort==cohort, ]
  stopifnot(nrow(row)==1L)
  for (metric in c("ci","p")) {
   valid <- mask & if(metric=="ci") ci_ok(a)&ci_ok(b) else p_ok(a)&p_ok(b)
   disagree <- if(metric=="ci") xor(a$low[valid]>0|a$high[valid]<0,b$low[valid]>0|b$high[valid]<0) else
     xor(a$p[valid]<.05,b$p[valid]<.05)
   equal(row[[paste0(metric,"_available_n")]],sum(valid))
   equal(row[[paste0(metric,"_difference_n")]],sum(disagree))
   equal(row[[paste0(metric,"_difference_conditional_probability")]],if(sum(valid))mean(disagree)else NA_real_)
   delta <- if(metric=="ci") pmax(abs(a$low[valid]-b$low[valid]),abs(a$high[valid]-b$high[valid])) else abs(a$p[valid]-b$p[valid])
   equal(unlist(row[paste0(metric,c("_delta_median","_delta_p99","_delta_max"))]),quantiles(delta))
  }
  equal(unlist(row[c("beta_delta_median","beta_delta_p99")]),quantiles(abs(a$beta[mask]-b$beta[mask]))[1:2])
  checked <- checked+1L
 }
}
stopifnot(checked==26L*8L*4L)
record <- c(sprintf("raw_fit_validation=%d comparison/design/cohort combinations",checked),
 "exact_counts=available, disagreements and fractions independently matched",
 "exact_numeric_summaries=CI, p-value and coefficient quantiles matched R type 7; maxima matched")
writeLines(record,"results/checks/table_counts.txt")
cat(paste(record,collapse="\n"),"\n")

# Read saved summaries only; no data generation or model fitting.
design <- read.csv("data/scenarios.csv", stringsAsFactors = FALSE)
pairs <- read.csv("results/tables/pair_summary.csv", stringsAsFactors = FALSE)
selected <- design
z <- merge(selected, subset(pairs, comparison == "default_logistf_sas" & cohort == "all_outputs"), by = "scenario_id")
z <- z[order(z$n, z$ratio), ]
stopifnot(nrow(z) == 8L, !anyDuplicated(z$scenario_id),
          all(z$n_exposed + z$n_unexposed == z$n),
          all(z$n_unexposed / z$n_exposed == z$ratio),
          max(abs(z$ci_difference_conditional_probability - z$ci_difference_joint_probability / z$ci_available_probability)) < 1e-12,
          max(abs(z$ci_difference_conditional_probability - z$p_difference_conditional_probability)) < 1e-12)
other <- subset(pairs, scenario_id %in% selected$scenario_id & cohort == "all_outputs" &
                  comparison %in% c("default_brglm2_sas", "aligned_PL_logistf_sas"))
stopifnot(nrow(other) == 16L, all(other$ci_difference_conditional_probability == 0))
out <- z[c("scenario_id", "n", "n_exposed", "n_unexposed", "ratio",
           "comparison", "cohort", "ci_available_probability", "ci_difference_conditional_probability",
           "p_difference_conditional_probability", "ci_difference_n", "ci_available_n")]
out$ci_difference_percent <- 100 * out$ci_difference_conditional_probability
write.csv(out, "results/tables/exact_size_allocation.csv", row.names = FALSE)
print(out[c("n", "n_exposed", "n_unexposed", "ci_difference_percent")], row.names = FALSE)
cat("Verified eight designs, conditional denominators, equal default CI/p decision differences, and zero CI decision differences in the two other comparisons.\n")
