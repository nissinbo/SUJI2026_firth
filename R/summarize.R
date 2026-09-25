source("R/project_helpers.R", local = TRUE)
set_project_locale()
source("R/import_results.R", local = TRUE)
source("R/data_helpers.R", local = TRUE)
source("R/analysis_helpers.R", local = TRUE)
source("R/comparison_settings.R", local = TRUE)

output_dir <- "results/tables"
ensure_directories(output_dir)
write_result <- function(x,name) write.csv(x,file.path(output_dir,paste0(name,".csv")),row.names=FALSE,na="")
tables <- read.csv("data/tables.csv",stringsAsFactors=FALSE)
covmeta <- readRDS("results/raw/covariate_diagnostics.rds")
meta <- combine_metadata(tables, covmeta)
eligible_ids <- meta$config_id[meta$data_state=="ok"]
stopifnot(length(eligible_ids)==41877L,!anyDuplicated(meta$config_id))
all <- load_model_fits(eligible_ids, rebuild = "--rebuild" %in% commandArgs(TRUE))
by_method <- split(all,all$key);rm(all);gc()
get_method <- function(key) {
 z <- by_method[[key]];stopifnot(!is.null(z));z[match(meta$config_id,z$config_id),]
}
contrasts <- comparison_definitions()
write_result(contrasts,"contrast_definitions")
common <- meta$data_state=="ok"
common_inference <- common
# FL has no convergence flag, so it cannot define this common cohort.
known_methods <- names(by_method)[!grepl("^sas_fl", names(by_method))]
for (key in known_methods) {
  result <- get_method(key)
  converged <- confirmed_convergence(result$converged)
  common <- common & converged
  common_inference <- common_inference & converged & ci_ok(result) & p_ok(result)
}
write_result(data.frame(config_id=meta$config_id,common_point_converged=common,common_inference_available=common_inference),"common_cohort")
# Each distinct table counts once. Generate the eight design records from the
# saved tables; simulated data and fits are never modified here.
scenarios <- unique(tables[c("n", "n_exposed", "n_unexposed", "design_id")])
scenarios$ratio <- scenarios$n_unexposed / scenarios$n_exposed
scenarios$scenario_id <- scenarios$design_id
stopifnot(nrow(scenarios)==8L)
write.csv(scenarios,"data/scenarios.csv",row.names=FALSE)
# Retain the existing *_probability column names for compatibility. They now
# represent equal-count proportions for exact tables, as for simulated datasets.
covdesign <- read.csv("data/covariate_scenarios.csv",stringsAsFactors=FALSE)
weights <- lapply(seq_len(nrow(scenarios)),function(i) {
 s<-scenarios[i,];idx<-which(meta$design_id==s$design_id)
 w<-rep(1/length(idx),length(idx))
 stopifnot(abs(sum(w)-1)<1e-12)
 list(scenario_id=s$scenario_id,kind="exact",idx=idx,w=w,B=NA_integer_)
})
weights <- c(weights,lapply(seq_len(nrow(covdesign)),function(i) {
 s<-covdesign[i,];idx<-which(meta$design_id==s$scenario_id)
 stopifnot(length(idx)==s$repetitions)
 list(scenario_id=s$scenario_id,kind="simulation",idx=idx,w=rep(1/s$repetitions,s$repetitions),B=s$repetitions)
}))
write_result(do.call(rbind,lapply(weights,function(q) {
 z<-meta[q$idx,];data.frame(scenario_id=q$scenario_id,kind=q$kind,n=nrow(z),
  exclusion_probability=sum(q$w[z$data_state!="ok"]),
  complete_separation_probability=sum(q$w[z$separation=="complete_separation"]),
  quasi_separation_probability=sum(q$w[z$separation=="quasi_separation"]),
  common_point_converged_probability=sum(q$w[common[q$idx]]),expected_events=sum(q$w*z$events))
})),"scenario_diagnostics")
pair_rows <- list()
position <- 0L
for (i in seq_len(nrow(contrasts))) {
 sp <- contrasts[i, ]
 a <- get_method(sp$method_A)
 b <- get_method(sp$method_B)
 ci <- ci_ok(a) & ci_ok(b)
 pv <- p_ok(a) & p_ok(b)
 cd <- ci & xor(a$low > 0 | a$high < 0, b$low > 0 | b$high < 0)
 pd <- pv & xor(a$p < .05, b$p < .05)
 cd[is.na(cd)] <- FALSE
 pd[is.na(pd)] <- FALSE
 both_converged <- confirmed_convergence(a$converged) & confirmed_convergence(b$converged)
 any_bad <- confirmed_convergence(!a$converged) | confirmed_convergence(!b$converged)
 dc <- pmax(abs(a$low - b$low), abs(a$high - b$high))
 dp <- abs(a$p - b$p)
 db <- abs(a$beta - b$beta)
 cohorts <- list(
   all_outputs = rep(TRUE, nrow(meta)),
   both_point_converged = both_converged,
   common_point_converged = common,
   common_inference_available = common_inference
 )
 for (q in weights) for (cohort in names(cohorts)) {
  idx <- q$idx
  w <- q$w
  mask <- cohorts[[cohort]][idx]
  keepci <- mask & ci[idx]
  keepp <- mask & pv[idx]
  cden <- sum(w[keepci])
  pden <- sum(w[keepp])
  cnum <- sum(w[mask & cd[idx]])
  pnum <- sum(w[mask & pd[idx]])
  cp <- if (cden > 0) cnum / cden else NA_real_
  pp <- if (pden > 0) pnum / pden else NA_real_
  cmc <- if (q$kind == "simulation" && cden > 0) sqrt(cp * (1 - cp) / sum(keepci)) else NA_real_
  pmc <- if (q$kind == "simulation" && pden > 0) sqrt(pp * (1 - pp) / sum(keepp)) else NA_real_
  position<-position+1L
  pair_rows[[position]]<-data.frame(scenario_id=q$scenario_id,kind=q$kind,comparison=sp$comparison,factor=sp$factor,cohort=cohort,
   total_configurations=length(idx),cohort_probability=sum(w[mask]),
   excluded_probability=sum(w[meta$data_state[idx]!="ok"]),
   ci_available_probability=cden,p_available_probability=pden,
   ci_difference_joint_probability=cnum,p_difference_joint_probability=pnum,
   ci_difference_conditional_probability=cp,p_difference_conditional_probability=pp,
   ci_difference_mcse=cmc,p_difference_mcse=pmc,
   ci_difference_n=sum(mask&cd[idx]),p_difference_n=sum(mask&pd[idx]),
   ci_available_n=sum(keepci),p_available_n=sum(keepp),
   nonconverged_output_probability=sum(w[mask&any_bad[idx]&ci[idx]]),
   convergence_unknown_output_probability=sum(w[mask&ci[idx]&!any_bad[idx]&!both_converged[idx]]),
   ci_delta_median=summary_quantile(dc[idx][keepci],w[keepci],.5,q$kind),ci_delta_p99=summary_quantile(dc[idx][keepci],w[keepci],.99,q$kind),
   ci_delta_max=if(any(keepci))max(dc[idx][keepci])else NA_real_,
   p_delta_median=summary_quantile(dp[idx][keepp],w[keepp],.5,q$kind),p_delta_p99=summary_quantile(dp[idx][keepp],w[keepp],.99,q$kind),
   p_delta_max=if(any(keepp))max(dp[idx][keepp])else NA_real_,
   beta_delta_median=summary_quantile(db[idx][mask],w[mask],.5,q$kind),beta_delta_p99=summary_quantile(db[idx][mask],w[mask],.99,q$kind))
 }
 cat("Summarized",sp$comparison,"\n")
}
pair_summary<-do.call(rbind,pair_rows);write_result(pair_summary,"pair_summary")
method_rows<-list();position<-0L
for(key in names(by_method)) {
 z<-get_method(key)
 for(q in weights) {
  idx<-q$idx;w<-q$w;eligible<-meta$data_state[idx]=="ok";point<-is.finite(z$beta[idx]);ci<-ci_ok(z)[idx];pv<-p_ok(z)[idx]
  known<-!is.na(z$converged[idx]);conv<-known&z$converged[idx];fail<-known&!z$converged[idx]
  position<-position+1L
  method_rows[[position]]<-data.frame(scenario_id=q$scenario_id,kind=q$kind,method=key,
   excluded_probability=sum(w[!eligible]),point_missing_probability=sum(w[eligible&!point]),
   ci_missing_probability=sum(w[eligible&!ci]),p_missing_probability=sum(w[eligible&!pv]),
   point_converged_probability=sum(w[conv]),point_nonconverged_probability=sum(w[fail]),
   convergence_unknown_probability=sum(w[eligible&!known]),
   ci_diagnostic_failed_probability=sum(w[!is.na(z$ci_converged[idx])&!z$ci_converged[idx]&eligible]),
   plrt_iteration_limit_probability=sum(w[!is.na(z$plrt_limit[idx])&z$plrt_limit[idx]&eligible]),
   closed_form_error_over_1e8_probability=if(q$kind=="exact")sum(w[eligible&point&abs(z$beta[idx]-meta$beta_reference[idx])>1e-8])else NA_real_)
 }
}
write_result(do.call(rbind,method_rows),"method_status")
writeLines(c(paste0("completed=",format(Sys.time())),paste0("eligible_datasets=",length(eligible_ids)),
 "enumerated_tables=22437; eligible=22421; one_outcome=16",
 "covariate_datasets=20000; eligible=19456; one_outcome=544",
 paste0("integrated_rows=",length(eligible_ids)*length(by_method),"; unique methods=",length(by_method),"; target ID sets verified for every method"),
 paste0("scope=","all four implementations"),
 "exact_aggregation=each table counted once; 8 group-size designs; no binomial weighting",
 "exact_quantiles=R quantile type 7; simulation quantiles unchanged (empirical inverse CDF)",
 "decision_comparison=unrounded CI excludes zero; p < 0.05",
 "proportion_denominator=all configurations/replicates for availability; available pairs for disagreements; legacy column suffix probability retained",
 "MCSE=binomial SE over available independent replicates; zero observed events does not prove zero population probability",
 "common_cohort=all 13 R and 7 SAS PROC point fits converge; FL convergence not available",
 "common_inference_available=add finite ordered CIs and p in [0,1] for every R/SAS PROC variant; same subset for factor ranking"),"results/checks/aggregation.txt")
cat("Analysis summaries completed.\n")
