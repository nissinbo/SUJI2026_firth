/* Unmodified %FL, with recursive isolation of groups that raise numerical errors.
   A failing individual remains a missing result, never silently omitted. */
options validvarname=any nomprint nomlogic nosymbolgen nonotes nosyntaxcheck;
%let root=%sysget(FIRTH_PROJECT_ROOT);
%let phase=%sysget(FIRTH_PHASE);
%let out=&root/results/raw/sas_&phase._fl;
%let range_from=%sysget(FIRTH_FROM);
%let range_to=%sysget(FIRTH_TO);
%if %length(&range_from)=0 %then %do;%let range_from=1;%end;
%if %length(&range_to)=0 %then %do;%let range_to=420;%end;
%let analysis_counter=0;
%include "&root/sas/vendor/fl.sas";
%macro isolate(data=,depth=0);
 %local total half dsid hascols rc found;
 data try_input;set &data;by config_id;retain _seq 0;
  if first.config_id then _seq+1;group_id=_seq;drop _seq;
 run;
 proc sort data=try_input(keep=config_id group_id) out=try_map nodupkey;by group_id;run;
 proc sql noprint;select count(*) into :total trimmed from try_map;quit;
 /* Preserve recursion inputs before the macro uses its own working datasets. */
 data level&depth;set try_input;run;
 proc datasets lib=work nolist;delete fltry flest;quit;
 %fl(data=try_input,y=y,varlist=&rhs,by=group_id,outtab=fltry,outest=flest,
     print=0,notes=0,pl=1,pl1=1,epsilon=&epsilon,maxit=&maxit,plmaxit=&plmaxit)
 options nosyntaxcheck obs=max firstobs=1;
 %let syscc=0;
 %let hascols=0;%let found=0;
 %if %sysfunc(exist(fltry)) %then %do;
  %let dsid=%sysfunc(open(fltry));
  %if &dsid>0 %then %do;
   %if %sysfunc(varnum(&dsid,_NAME_))>0 and %sysfunc(varnum(&dsid,BETA))>0 %then %let hascols=1;
   %let rc=%sysfunc(close(&dsid));
  %end;
 %end;
 %if &hascols=1 %then %do;
  data extracted;set fltry;where upcase(_NAME_)='EXPOSURE' and _VAR_=1;run;
  proc sql noprint;select count(*) into :found trimmed from extracted;quit;
 %end;
 %if &found=&total %then %do;
  data valid;
   merge try_map extracted;by group_id;
   length extraction_status $24;extraction_status='output_available';
   keep config_id BETA STDERR CI_LO CI_UP P_VALUE _ITER_ extraction_status;
  run;
  proc append base=collected data=valid force;run;
 %end;
 %else %if &total>1 %then %do;
  %let half=%sysfunc(floor(%sysevalf(&total/2)));
  data left&depth right&depth;set level&depth;
   if group_id<=&half then output left&depth;else output right&depth;
  run;
  %put FL_ISOLATE kind=&kind batch=&batch setting=&setting depth=&depth datasets=&total;
  %isolate(data=left&depth,depth=%eval(&depth+1))
  %isolate(data=right&depth,depth=%eval(&depth+1))
 %end;
 %else %do;
  data invalid;
   set level&depth(obs=1);length extraction_status $24;
   BETA=.;STDERR=.;CI_LO=.;CI_UP=.;P_VALUE=.;_ITER_=.;extraction_status='numeric_failure';
   keep config_id BETA STDERR CI_LO CI_UP P_VALUE _ITER_ extraction_status;
  run;
  proc append base=collected data=invalid force;run;
  data _null_;set invalid;put 'FL_INDIVIDUAL_FAILURE ' config_id " setting=&setting";run;
 %end;
%mend;
%macro fit(setting=,epsilon=,maxit=,plmaxit=);
 %if %sysfunc(fileexist(&out/&kind._&batch._&setting..csv)) %then %return;
 data collected;
  length config_id $32 extraction_status $24;
  BETA=.;STDERR=.;CI_LO=.;CI_UP=.;P_VALUE=.;_ITER_=.;stop;
 run;
 %isolate(data=batch_input)
 proc sort data=collected;by config_id;run;
 data one;
  length setting $8 beta_hex se_hex p_hex pl_low_hex pl_high_hex $16;
  set collected;setting="&setting";point_limit=&maxit;pl_limit=&plmaxit;
  beta_hex=put(BETA,hex16.);se_hex=put(STDERR,hex16.);p_hex=put(P_VALUE,hex16.);
  pl_low_hex=put(CI_LO,hex16.);pl_high_hex=put(CI_UP,hex16.);
 run;
 proc sql noprint;select count(*),count(distinct config_id) into :got trimmed,:unique_got trimmed from one;quit;
 %if &got ne &expected or &unique_got ne &expected %then %do;
  %put ERROR: FL wrapper lost or duplicated target rows.;%abort cancel;
 %end;
 filename result "&out/&kind._&batch._&setting..csv" encoding='utf-8';
 proc export data=one outfile=result dbms=csv replace;run;
 filename result clear;
%mend;
%macro run_kind(kind);
 data shared;
  length config_id $32 x1_hex x2_hex $16;
  infile "&root/data/sas_&phase._&kind..csv" dsd firstobs=2 truncover encoding='utf-8';
  input batch config_id :$32. row_id y exposure x1_hex :$16. x2_hex :$16.;
  x1=input(x1_hex,hex16.);x2=input(x2_hex,hex16.);
 run;
 proc sort data=shared;by batch config_id row_id;run;
 proc sql noprint;select max(batch) into :nbatch trimmed from shared;quit;
 %if &kind=enum %then %let rhs=exposure;
 %else %let rhs=exposure x1 x2;
 %do batch=1 %to &nbatch;
  %let analysis_counter=%eval(&analysis_counter+1);
  %if &analysis_counter>=&range_from and &analysis_counter<=&range_to %then %do;
  data batch_input;set shared;where batch=&batch;run;
  proc sql noprint;select count(distinct config_id) into :expected trimmed from batch_input;quit;
  %fit(setting=default,epsilon=1e-4,maxit=50,plmaxit=50)
  %fit(setting=tight,epsilon=1e-12,maxit=5000,plmaxit=500)
  %put ANALYSIS_PROGRESS fl &phase &kind &batch / &nbatch;
  %end;
 %end;
%mend;
%run_kind(enum)
%run_kind(cov)
data _null_;file "&out/status_&range_from._&range_to..txt";
 put 'completed=TRUE' / "SAS=&sysvlong4" / 'engine=fl' / 'numerical_failures=isolated individually and retained as missing outputs';
run;
%put ANALYSIS_COMPLETE;
