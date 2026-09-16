# Shared design for equally counted 2x2 configurations and paired simulations.
Sys.setlocale("LC_CTYPE", "Japanese_Japan.utf8")
source("R/data_helpers.R", local = TRUE)

dir.create("data", recursive=TRUE, showWarnings=FALSE)
dir.create("results", recursive=TRUE, showWarnings=FALSE)
dir.create("results/raw", recursive=TRUE, showWarnings=FALSE)
dir.create("results/checks", recursive=TRUE, showWarnings=FALSE)
B <- 2500L
design <- expand.grid(n=c(20L,50L,100L,200L), ratio=c(1L,4L))
design$n_exposed <- design$n/(1+design$ratio)
design$n_unexposed <- design$n-design$n_exposed
design$design_id <- sprintf("E%03dR%d",design$n,design$ratio)
scenarios <- design
scenarios$scenario_id <- scenarios$design_id
write.csv(scenarios,"data/scenarios.csv",row.names=FALSE)
tables <- do.call(rbind,lapply(seq_len(nrow(design)),function(i) {
  s <- design[i,]; z <- expand.grid(a=0:s$n_exposed,c=0:s$n_unexposed)
  z$b <- s$n_exposed-z$a; z$d <- s$n_unexposed-z$c
  z$config_id <- sprintf("n%d_%d.%d.%d.%d",s$n,z$a,z$b,z$c,z$d)
  z$design_id <- s$design_id; z$scenario <- "ENUM"
  z$n <- s$n; z$n_exposed <- s$n_exposed; z$n_unexposed <- s$n_unexposed
  z$events <- z$a+z$c
  z$data_state <- ifelse(z$events %in% c(0,s$n),"one_outcome","ok")
  z$separation <- ifelse((z$a==0 & z$d==0)|(z$b==0 & z$c==0),"complete_separation",
    ifelse(z$a*z$b*z$c*z$d==0,"quasi_separation","none"))
  z$separation[z$data_state!="ok"] <- "one_outcome"
  z$sparsity_class <- z$separation
  z$beta_reference <- with(z,log((a+.5)*(d+.5)/((b+.5)*(c+.5))))
  z
}))
stopifnot(nrow(tables)==22437L,nrow(scenarios)==8L,!anyDuplicated(tables$config_id))
write.csv(tables,"data/tables.csv",row.names=FALSE)
# Two fixed backgrounds; correlation contrasts are within each background.
covdesign <- expand.grid(background=c("N40_balanced","N100_unbalanced"),rho=c(0,.5,.9,.99),stringsAsFactors=FALSE)
covdesign$n <- ifelse(covdesign$background=="N40_balanced",40L,100L)
covdesign$n_exposed <- 20L
covdesign$target_p0 <- .05; covdesign$OR <- 2
covdesign$gamma1 <- .5; covdesign$gamma2 <- .5
covdesign$scenario_id <- sprintf("C%dR%02d",covdesign$n,round(covdesign$rho*100))
covdesign$repetitions <- B
# Integrate the normal linear predictor to fix the marginal nonexposed risk.
covdesign$intercept <- vapply(seq_len(nrow(covdesign)),function(i) {
  s <- covdesign[i,]; sd_eta <- sqrt(s$gamma1^2+s$gamma2^2+2*s$rho*s$gamma1*s$gamma2)
  uniroot(function(b) integrate(function(z) plogis(b+sd_eta*z)*dnorm(z),-10,10,rel.tol=1e-10)$value-s$target_p0,c(-15,5),tol=1e-12)$root
},numeric(1))
write.csv(covdesign,"data/covariate_scenarios.csv",row.names=FALSE)
RNGkind("Mersenne-Twister","Inversion","Rejection")
metas <- vector("list",nrow(covdesign)*B)
for (i in seq_len(nrow(covdesign))) {
  s <- covdesign[i,]
  datasets <- vector("list", B)
  for (replicate in seq_len(B)) {
    d <- generate_covariate_data(s, replicate)
    id <- sprintf("%s_%04d",s$scenario_id,replicate)
    mm <- model.matrix(~ exposure+x1+x2,d)
    metas[[(i-1L)*B+replicate]] <- data.frame(config_id=id,design_id=s$scenario_id,scenario="COVARIATE",
      n=s$n,n_exposed=s$n_exposed,n_unexposed=s$n-s$n_exposed,a=NA,b=NA,c=NA,d=NA,
      events=sum(d$y),data_state=ifelse(sum(d$y)%in%c(0,s$n),"one_outcome","ok"),
      separation="not_classified",sparsity_class="covariate",beta_reference=NA_real_,
      replicate=replicate,rho=s$rho,sample_correlation=cor(d$x1,d$x2),condition_number=kappa(mm,exact=TRUE),
      model_rank=qr(mm)$rank)
    datasets[[replicate]] <- cbind(config_id=id,d)
  }
  dat <- do.call(rbind,datasets)
  saveRDS(dat,sprintf("data/%s.rds",s$scenario_id))
}
covmeta <- do.call(rbind,metas)
stopifnot(nrow(covmeta)==20000L,all(covmeta$model_rank==4L),!anyDuplicated(covmeta$config_id))
saveRDS(covmeta,"data/covariate_metadata.rds")
writeLines(c("design_version=2026-09-09",paste0("repetitions=",B),
 "aggregation_revision=2026-09-15; enumeration=8 designs, each table counted once, no binomial weighting",
 "one_outcome=recorded as design exclusion; retained in availability denominator",
 "covariate_p0=marginal nonexposed probability 0.05; intercept calibrated per rho",
 "covariate_OR=conditional exposure odds ratio 2; fixed gamma1=gamma2=0.5",
 "common_random_numbers=same x1,z,uniform outcome draws within n and replicate",
 "covariate_backgrounds=40 with 20:20; 100 with 20:80; not a size-only contrast"),"data/design_record.txt")
cat("Design created:",nrow(tables),"tables,",nrow(scenarios),"enumeration designs,",nrow(covmeta),"covariate datasets\n")
