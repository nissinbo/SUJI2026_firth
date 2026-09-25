# Fit the presentation model settings on shared inputs; resume saved chunks.
source("R/project_helpers.R", local = TRUE)
set_project_locale()
source("R/utils.R", local = TRUE)
source("R/data_helpers.R", local = TRUE)
source("R/fit_models.R", local = TRUE)
args <- commandArgs(TRUE)
phase <- if ("--pilot" %in% args) "pilot" else "full"
shard <- integer_argument(args, "shard", 1L)
nshards <- if (phase=="pilot") 1L else 2L
stopifnot(shard %in% seq_len(nshards))
tables <- read.csv("data/tables.csv",stringsAsFactors=FALSE)
covmeta <- readRDS("data/covariate_metadata.rds")
metas <- combine_metadata(tables, covmeta)
metas <- metas[metas$data_state=="ok",]
if (phase=="pilot") {
  ids <- pilot_ids(tables, covmeta)
  metas <- metas[metas$config_id%in%ids,]
} else {
  metas <- metas[(seq_len(nrow(metas))-1L)%%nshards+1L==shard,]
}
covariates <- read_covariate_inputs()
outdir <- sprintf("results/raw/r_%s_%d",phase,shard)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
started <- Sys.time()
chunks <- split(seq_len(nrow(metas)),ceiling(seq_len(nrow(metas))/100L))
start_arg <- integer_argument(args, "from", 1L)
end_arg <- integer_argument(args, "to", length(chunks))
stopifnot(start_arg >= 1L, start_arg <= length(chunks), end_arg >= start_arg)
for (chunk in seq.int(start_arg,min(end_arg,length(chunks)))) {
  outpath <- file.path(outdir,sprintf("chunk_%04d.rds",chunk))
  if(file.exists(outpath)) next
  audits <- outputs <- list(); pos <- 0L
  for (i in chunks[[chunk]]) {
    dat <- model_data(metas[i, ], covariates)
    for (v in variants) {
      got <- fit_one(metas[i, ], v, dat)
      pos <- pos + 1L
      audits[[pos]] <- got$audit
      outputs[[pos]] <- got$result
    }
  }
  saveRDS(list(audit=do.call(rbind,audits),results=do.call(rbind,outputs)),outpath)
  cat(sprintf("%s shard=%d chunk=%d/%d datasets=%d elapsed=%.1fs\n",format(Sys.time()),shard,chunk,length(chunks),max(chunks[[chunk]]),as.numeric(difftime(Sys.time(),started,units="secs"))))
  flush.console()
}
writeLines(c(paste0("completed=",format(Sys.time())),paste0("datasets=",nrow(metas)),paste0("variants=",length(variants)),
  paste0("logistf=",packageVersion("logistf")),paste0("brglm2=",packageVersion("brglm2")),
  paste0("factor_source_md5=",tools::md5sum("R/fit_models.R"))),file.path(outdir,sprintf("status_%03d_%03d.txt",start_arg,end_arg)))
