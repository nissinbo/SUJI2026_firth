/* PROC LOGISTIC: seven settings, shared binary-exact inputs. */
options validvarname=any nomprint nomlogic nosymbolgen nonotes;
%let root=%sysget(FIRTH_PROJECT_ROOT);
%let phase=%sysget(FIRTH_PHASE);
%let engine=%sysget(FIRTH_ENGINE);
%let out=&root/results/raw/sas_&phase._&engine;
%macro export(ds,name);
 filename result "&out/&name..csv" encoding='utf-8';
 proc export data=&ds outfile=result dbms=csv replace; run;
 filename result clear;
%mend;
%macro procfit(variant=,opts=);
 proc datasets lib=work nolist; delete pe wc pc cs ih; quit;
 ods exclude all;
 ods output ParameterEstimates=pe ClparmWald=wc ClparmPL=pc ConvergenceStatus=cs IterHistory=ih;
 proc logistic data=batch_input;
  by config_id;
  model y(event='1')=&rhs / firth clparm=both itprint &opts;
 run;
 ods output close;
 data pe1; set pe; where upcase(Variable)='EXPOSURE';
  keep config_id Estimate StdErr ProbChiSq; rename Estimate=beta StdErr=se ProbChiSq=p;
 run;
 data wc1; set wc; where upcase(Parameter)='EXPOSURE';
  keep config_id LowerCL UpperCL; rename LowerCL=wald_low UpperCL=wald_high;
 run;
 data pc1; set pc; where upcase(Parameter)='EXPOSURE';
  keep config_id LowerCL UpperCL; rename LowerCL=pl_low UpperCL=pl_high;
 run;
 data cs1; set cs; where NullModel=0; keep config_id Status Reason; run;
 proc sql; create table ih1 as select config_id,max(Iteration) as iterations from ih where not missing(exposure) group by config_id; quit;
 data one;
  length variant $32 beta_hex se_hex p_hex wald_low_hex wald_high_hex pl_low_hex pl_high_hex $16;
  merge targets pe1 wc1 pc1 cs1 ih1; by config_id;
  variant="&variant";
  beta_hex=put(beta,hex16.);se_hex=put(se,hex16.);p_hex=put(p,hex16.);
  wald_low_hex=put(wald_low,hex16.);wald_high_hex=put(wald_high,hex16.);
  pl_low_hex=put(pl_low,hex16.);pl_high_hex=put(pl_high,hex16.);
 run;
 %export(one,&kind._&batch._&variant)
%mend;
%macro run_kind(kind);
 data shared;
  length config_id $32 x1_hex x2_hex $16;
  infile "&root/data/sas_&phase._&kind..csv" dsd firstobs=2 truncover encoding='utf-8';
  input batch config_id :$32. row_id y exposure x1_hex :$16. x2_hex :$16.;
  x1=input(x1_hex,hex16.);x2=input(x2_hex,hex16.);
 run;
 proc sort data=shared; by batch config_id row_id; run;
 proc sql noprint; select max(batch) into :nbatch trimmed from shared; quit;
 %if &kind=enum %then %let rhs=exposure;
 %else %let rhs=exposure x1 x2;
 %do batch=1 %to &nbatch;
  data batch_input;set shared;where batch=&batch;run;
  proc sort data=batch_input(keep=config_id) out=targets nodupkey;by config_id;run;

   %procfit(variant=fisher_default)
   %procfit(variant=newton_default,opts=technique=newton)
   %procfit(variant=fisher_cap,opts=maxiter=5000)
   %procfit(variant=fisher_gradient,opts=maxiter=5000 gconv=1e-12)
   %procfit(variant=fisher_parameter,opts=maxiter=5000 xconv=1e-12)
   %procfit(variant=newton_gradient,opts=technique=newton maxiter=5000 gconv=1e-12)
   %procfit(variant=fisher_pl,opts=maxiter=5000 gconv=1e-12 plconv=1e-12)
  %put ANALYSIS_PROGRESS &engine &phase &kind &batch / &nbatch;
 %end;
%mend;
%run_kind(enum)
%run_kind(cov)
data _null_; file "&out/status.txt"; put "completed=TRUE" / "SAS=&sysvlong4" / "phase=&phase" / "engine=&engine";run;
%put ANALYSIS_COMPLETE;
