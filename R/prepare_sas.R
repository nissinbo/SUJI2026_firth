Sys.setlocale("LC_CTYPE","Japanese_Japan.utf8")
source("R/utils.R", local = TRUE)
source("R/data_helpers.R", local = TRUE)
args <- commandArgs(TRUE)
phase <- if("--pilot" %in% args) "pilot" else "full"
tables <- read.csv("data/tables.csv")
covmeta <- readRDS("data/covariate_metadata.rds")
tables <- tables[tables$data_state=="ok",]
covmeta <- covmeta[covmeta$data_state=="ok",]
if(phase=="pilot") {
  ids <- pilot_ids(tables, covmeta)
  tables <- tables[tables$config_id%in%ids,]
  covmeta <- covmeta[covmeta$config_id %in% ids, ]
}
covariates <- read_covariate_inputs()
for (kind in c("enum","cov")) {
  meta <- if(kind=="enum") tables else covmeta
  meta <- meta[order(meta$config_id),]
  meta$batch <- ceiling(seq_len(nrow(meta))/100L)
  outfile <- sprintf("data/sas_%s_%s.csv",phase,kind)
  con <- file(outfile,"wt",encoding="UTF-8")
  writeLines("batch,config_id,row_id,y,exposure,x1_hex,x2_hex",con)
  for (chunk in split(seq_len(nrow(meta)),meta$batch)) {
    dat <- do.call(rbind,lapply(chunk,function(i) {
      z <- meta[i,]
      d <- model_data(z, covariates)
      data.frame(batch=z$batch,config_id=z$config_id,row_id=d$row_id,y=d$y,exposure=d$exposure,
                 x1_hex=double_to_hex(d$x1),x2_hex=double_to_hex(d$x2))
    }))
    write.table(dat,con,sep=",",row.names=FALSE,col.names=FALSE,quote=TRUE)
  }
  close(con)
  cat(phase,kind,nrow(meta),"datasets written\n")
}
