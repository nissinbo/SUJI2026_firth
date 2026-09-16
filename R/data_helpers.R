# Shared input selection and reconstruction. Sourcing this file performs no I/O.
combine_metadata <- function(tables, covariates) {
  columns <- intersect(names(tables), names(covariates))
  metadata <- rbind(tables[, columns], covariates[, columns])
  rownames(metadata) <- NULL
  metadata
}

pilot_ids <- function(tables, covariates) {
  tables <- tables[tables$data_state == "ok", ]
  table_ids <- unlist(lapply(split(tables, tables$design_id), function(group) {
    group$config_id[unique(round(seq(1, nrow(group), length.out = 3)))]
  }), use.names = FALSE)
  c(table_ids, covariates$config_id[
    covariates$replicate <= 5 & covariates$data_state == "ok"
  ])
}

table_data <- function(meta) {
  data.frame(
    row_id = seq_len(meta$n),
    y = c(rep(1L, meta$a), rep(0L, meta$b), rep(1L, meta$c), rep(0L, meta$d)),
    exposure = c(rep(1L, meta$n_exposed), rep(0L, meta$n_unexposed)),
    x1 = 0, x2 = 0
  )
}

read_covariate_inputs <- function(directory = "data") {
  files <- list.files(directory, pattern = "^C[0-9]+R[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No covariate inputs found in ", directory)
  data <- do.call(rbind, lapply(files, readRDS))
  list(data = data, index = split(seq_len(nrow(data)), data$config_id))
}

model_data <- function(meta, covariates) {
  if (meta$scenario != "COVARIATE") return(table_data(meta))
  index <- covariates$index[[meta$config_id]]
  if (is.null(index)) stop("Missing covariate dataset: ", meta$config_id)
  covariates$data[index, ]
}

model_formula <- function(meta) {
  if (meta$scenario == "COVARIATE") y ~ exposure + x1 + x2 else y ~ exposure
}

# Match the saved design's random-number sequence across correlation settings.
generate_covariate_data <- function(scenario, replicate) {
  n <- scenario$n
  set.seed(20260909L + n * 10000L + replicate)
  x1 <- rnorm(n)
  z <- rnorm(n)
  u <- runif(n)
  exposure <- c(rep(1L, scenario$n_exposed), rep(0L, n - scenario$n_exposed))
  x2 <- scenario$rho * x1 + sqrt(1 - scenario$rho^2) * z
  probability <- plogis(scenario$intercept + log(scenario$OR) * exposure +
                         scenario$gamma1 * x1 + scenario$gamma2 * x2)
  data.frame(row_id = seq_len(n), y = as.integer(u < probability), exposure, x1, x2)
}
