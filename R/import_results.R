# Preserve binary SAS values; decimal PROC EXPORT display formats are not used.
hex_double <- function(x) {
  missing <- is.na(x) | nchar(x)!=16L | x=="2E00000000000000"
  x[missing] <- "0000000000000000"
  bytes <- vapply(seq(1,15,2),function(i) strtoi(substr(x,i,i+1),16L),integer(length(x)))
  ans <- readBin(as.raw(t(bytes)),"double",n=length(x),size=8L,endian="big")
  ans[missing] <- NA_real_; ans
}
read_sas_csv <- function(path) {
  z <- read.csv(path,colClasses="character",encoding="UTF-8",check.names=FALSE)
  names(z)[1] <- sub("^\ufeff","",names(z)[1])
  for (name in grep("_hex$",names(z),value=TRUE)) z[[sub("_hex$","",name)]] <- hex_double(z[[name]])
  z
}
standard_r <- function(z) {
  data.frame(config_id=z$config_id,engine=z$implementation,variant=z$variant_id,inference=z$inference,
    beta=z$beta,se=z$standard_error,low=z$ci_low,high=z$ci_high,p=z$p_value,
    converged=z$software_converged,ci_converged=z$ci_converged,
    point_limit=z$full_iteration_limit,plrt_limit=z$plrt_iteration_limit,iterations=z$iterations,
    status=z$fit_status)
}
standard_sas <- function(path,engine) {
  z <- read_sas_csv(path)
  if(engine=="proc") do.call(rbind,lapply(c("Wald","PL"),function(inference) {
    data.frame(config_id=z$config_id,engine="sas_proc",variant=z$variant,inference=inference,
      beta=z$beta,se=z$se,low=z[[if(inference=="Wald") "wald_low" else "pl_low"]],
      high=z[[if(inference=="Wald") "wald_high" else "pl_high"]],p=z$p,
      converged=as.numeric(z$Status)==0,ci_converged=if(inference=="Wald") as.numeric(z$Status)==0 else NA,
      point_limit=NA,plrt_limit=NA,iterations=as.numeric(z$iterations),status=z$Reason)
  })) else data.frame(config_id=z$config_id,engine="sas_fl",variant=z$setting,inference="PL",
    beta=z$beta,se=z$se,low=z$pl_low,high=z$pl_high,p=z$p,
    converged=NA,ci_converged=NA,point_limit=as.numeric(z$`_ITER_`)>=as.numeric(z$point_limit),
    plrt_limit=NA,iterations=as.numeric(z$`_ITER_`),status=ifelse(is.finite(z$beta),"output_available_convergence_unknown","missing"))
}

# Import full runs only; pilot outputs are deliberately excluded.
load_model_fits <- function(eligible_ids, rebuild = FALSE, cache = "results/model_fits.rds") {
  if (rebuild || !file.exists(cache)) {
   rfiles <- unlist(lapply(1:2,function(k)list.files(sprintf("results/raw/r_full_%d",k),pattern="^chunk_[0-9]+.rds$",full.names=TRUE)))
   stopifnot(length(rfiles)==420L)
   message("Importing 420 R result batches")
   rparts <- lapply(rfiles,function(f)standard_r(readRDS(f)$results))
   r <- do.call(rbind,rparts);rm(rparts);gc()
   stopifnot(nrow(r)==length(eligible_ids)*20L)
   procfiles <- list.files("results/raw/sas_full_proc",pattern="^[a-z]+_[0-9]+_.*csv$",full.names=TRUE)
   flfiles <- list.files("results/raw/sas_full_fl",pattern="^[a-z]+_[0-9]+_.*csv$",full.names=TRUE)
   stopifnot(length(procfiles)==420L*7L)
   stopifnot(length(flfiles)==420L*2L)
   message("Importing 2,940 PROC LOGISTIC and 840 FL result batches")
   sp <- do.call(rbind,lapply(procfiles,standard_sas,engine="proc"))
   sf <- do.call(rbind,lapply(flfiles,standard_sas,engine="fl"))
   stopifnot(nrow(sp)==length(eligible_ids)*14L)
   stopifnot(nrow(sf)==length(eligible_ids)*2L)
   all <- rbind(r,sp,sf);rm(r,sp,sf);gc()
   all$key <- paste(all$engine,all$variant,all$inference,sep="|")
   stopifnot(!anyDuplicated(paste(all$config_id,all$key)),length(unique(all$key))==36L)
   for(key in unique(all$key)) stopifnot(setequal(all$config_id[all$key==key],eligible_ids))
   saveRDS(all,cache)
  } else all <- readRDS(cache)
  all
}
