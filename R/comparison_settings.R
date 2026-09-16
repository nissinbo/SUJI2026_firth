# Definitions of the 26 comparisons; no data loading or model fitting.
comparison_definitions <- function() {
  K <- function(engine,variant,inference="Wald") paste(engine,variant,inference,sep="|")
  S <- function(variant="fisher_default",inference="Wald") K("sas_proc",variant,inference)
  L <- function(variant="lf_nr_default",inference="PL") K("logistf",variant,inference)
  G <- function(variant="bg_asmean_default") K("brglm2",variant)
  F <- function(variant="default") K("sas_fl",variant,"PL")
  spec <- list(
   c("default_logistf_sas",L(),S(),"default"),
   c("default_brglm2_sas",G(),S(),"default"),
   c("default_fl_sas",F(),S(),"default"),
   c("aligned_PL_logistf_sas",L(),S(inference="PL"),"alignment"),
   c("aligned_Wald_logistf_sas",L(inference="Wald"),S(),"alignment"),
   c("default_fl_logistf",F(),L(),"alignment"),
   c("tight_PL_logistf_sas",L("lf_nr_fit_pl_tight"),S("fisher_pl","PL"),"alignment"),
   c("tight_Wald_brglm2_sas",G("bg_asmean_epsilon12"),S("fisher_gradient"),"alignment"),
   c("tight_fl_logistf",F("tight"),L("lf_nr_fit_pl_tight"),"alignment"),
   c("LF_inference",L(inference="Wald"),L(),"inference"),
   c("SAS_CI",S(),S(inference="PL"),"inference"),
   c("LF_algorithm_default",L(),L("lf_irls_default"),"algorithm"),
   c("LF_algorithm_tight",L("lf_nr_fit_pl_tight"),L("lf_irls_fit_pl_tight"),"algorithm"),
   c("LF_maxit",L(),L("lf_nr_maxit5000"),"iterations"),
   c("LF_threshold",L("lf_nr_maxit5000"),L("lf_nr_fit_tight"),"threshold"),
   c("LF_PL_control",L(),L("lf_nr_pl_tight"),"profile"),
   c("LF_PL_control_tight_fit",L("lf_nr_fit_tight"),L("lf_nr_fit_pl_tight"),"profile"),
   c("BG_maxit",G(),G("bg_asmean_maxit5000"),"iterations"),
   c("BG_threshold",G("bg_asmean_maxit5000"),G("bg_asmean_epsilon12"),"threshold"),
   c("SAS_algorithm_default",S(),S("newton_default"),"algorithm"),
   c("SAS_algorithm_tight",S("fisher_gradient"),S("newton_gradient"),"algorithm"),
   c("SAS_maxit",S(),S("fisher_cap"),"iterations"),
   c("SAS_gradient",S("fisher_cap"),S("fisher_gradient"),"threshold"),
   c("SAS_parameter",S("fisher_cap"),S("fisher_parameter"),"threshold"),
   c("SAS_PL_control",S("fisher_gradient","PL"),S("fisher_pl","PL"),"profile"),
   c("FL_controls",F(),F("tight"),"combined_control")
  )
  contrasts <- as.data.frame(do.call(rbind,spec),stringsAsFactors=FALSE)
  names(contrasts)<-c("comparison","method_A","method_B","factor")
  contrasts
}
