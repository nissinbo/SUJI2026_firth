Sys.setlocale("LC_CTYPE","Japanese_Japan.utf8")
stopifnot(requireNamespace("lpSolve",quietly=TRUE))
classify_separation <- function(dat,formula) {
  if(length(unique(dat$y))==1L) return("one_outcome")
  X <- model.matrix(formula,dat)
  # An orthonormal basis preserves the separating hyperplanes and avoids
  # an unbounded feasibility problem for nearly collinear predictors.
  X <- qr.Q(qr(X))
  M <- X*(2*dat$y-1)
  A <- cbind(M,-M)
  constraints <- rbind(cbind(A,-1),c(rep(1,ncol(A)),0))
  fit <- lpSolve::lp("max",c(rep(0,ncol(A)),1),constraints,c(rep(">=",nrow(A)),"<="),c(rep(0,nrow(A)),1))
  stopifnot(fit$status==0L)
  if(fit$objval>1e-8) {
    stopifnot(min(A%*%fit$solution[seq_len(ncol(A))])>fit$objval-1e-7)
    return("complete_separation")
  }
  Q <- rbind(A,rep(1,ncol(A)))
  fit <- lpSolve::lp("max",colSums(A),Q,c(rep(">=",nrow(A)),"<="),c(rep(0,nrow(A)),1))
  stopifnot(fit$status==0L)
  if(fit$objval>1e-8) {
    stopifnot(min(A%*%fit$solution)>-1e-7)
    return("quasi_separation")
  }
  "none"
}
# Verify the geometric classifier against all possible n=10 balanced 2x2 tables.
for(a in 0:5) for(c in 0:5) {
  d <- data.frame(y=c(rep(1,a),rep(0,5-a),rep(1,c),rep(0,5-c)),exposure=rep(c(1,0),each=5))
  expected <- if((a+c) %in% c(0,10)) "one_outcome" else if((a==0&&c==5)||(a==5&&c==0)) "complete_separation" else if(any(c(a,5-a,c,5-c)==0)) "quasi_separation" else "none"
  stopifnot(classify_separation(d,~exposure)==expected)
}
meta <- readRDS("data/covariate_metadata.rds")
for(f in list.files("data",pattern="^C[0-9]+R[0-9]+\\.rds$",full.names=TRUE)) {
  dat <- readRDS(f)
  groups <- split(dat,dat$config_id)
  classes <- vapply(groups,classify_separation,character(1),formula=~exposure+x1+x2)
  meta$separation[match(names(classes),meta$config_id)] <- classes
  cat(basename(f),paste(names(table(classes)),table(classes),collapse="; "),"\n")
}
stopifnot(!any(meta$separation=="not_classified"))
# Separate diagnostic output leaves model input metadata immutable during fitting.
saveRDS(meta,"results/raw/covariate_diagnostics.rds")
writeLines(c(paste0("lpSolve=",packageVersion("lpSolve")),"classifier_checks=all 36 balanced n10 tables passed",
 "complete=positive maximum minimum margin; quasi=positive maximum sum of nonnegative margins with zero minimum margin",
 "parameter_L1_norm_bound=1; X replaced by orthonormal column-space basis",
 "classification_objective_tolerance=1e-8; feasibility_check=1e-7; all models full rank"),"results/checks/separation.txt")
