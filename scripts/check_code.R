# Read-only regression checks for shared inputs and model wrappers.
# Run from the project root; full experiments are not rerun.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Japanese_Japan.utf8")
source("R/utils.R", local = TRUE)
source("R/data_helpers.R", local = TRUE)
source("R/analysis_helpers.R", local = TRUE)
source("R/import_results.R", local = TRUE)
source("R/comparison_settings.R", local = TRUE)
source("R/fit_models.R", local = TRUE)

for (file in c(list.files("R", "\\.R$", full.names = TRUE),
               list.files("scripts", "\\.R$", full.names = TRUE))) {
  invisible(parse(file, encoding = "UTF-8"))
}

tables <- read.csv("data/tables.csv")
covmeta <- readRDS("data/covariate_metadata.rds")
design <- read.csv("data/covariate_scenarios.csv")
covariates <- read_covariate_inputs()
metadata <- combine_metadata(tables, covmeta)
ids <- pilot_ids(tables, covmeta)
stopifnot(!anyDuplicated(ids), all(ids %in% metadata$config_id),
          all(metadata$data_state[match(ids, metadata$config_id)] == "ok"))

# R reconstruction must agree with the binary-exact SAS input, including row order.
for (kind in c("enum", "cov")) {
  input <- read_sas_csv(paste0("data/sas_full_", kind, ".csv"))
  index <- split(seq_len(nrow(input)), input$config_id)
  scenario <- if (kind == "enum") "ENUM" else "COVARIATE"
  kind_ids <- ids[metadata$scenario[match(ids, metadata$config_id)] == scenario]
  stopifnot(length(kind_ids) > 0L, all(kind_ids %in% names(index)))
  for (id in kind_ids) {
    meta <- metadata[match(id, metadata$config_id), ]
    expected <- model_data(meta, covariates)
    actual <- input[index[[id]], ]
    for (column in c("row_id", "y", "exposure", "x1", "x2")) {
      stopifnot(identical(as.numeric(actual[[column]]), as.numeric(expected[[column]])))
    }
  }
}
message("R and SAS input reconstruction agrees for every pilot dataset.")

# Same saved random inputs in all eight conditions, including high correlation.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (i in seq_len(nrow(design))) for (replicate in c(1L, 1250L, 2500L)) {
  scenario <- design[i, ]
  generated <- generate_covariate_data(scenario, replicate)
  id <- sprintf("%s_%04d", scenario$scenario_id, replicate)
  saved <- model_data(covmeta[match(id, covmeta$config_id), ], covariates)
  stopifnot(isTRUE(all.equal(generated, saved[names(generated)],
                            tolerance = 0, check.attributes = FALSE)))
}
message("Generated inputs match saved inputs in 24 scenario/replicate combinations.")

# A missing or inverted interval must not enter an available-output denominator.
edge <- data.frame(low = c(-1, 0, NA, -Inf, 2), high = c(1, 0, 1, Inf, 1),
                   p = c(0, 1, NA, -0.01, 1.01))
stopifnot(identical(ci_ok(edge), c(TRUE, TRUE, FALSE, FALSE, FALSE)),
          identical(p_ok(edge), c(TRUE, TRUE, FALSE, FALSE, FALSE)),
          identical(confirmed_convergence(c(TRUE, FALSE, NA)), c(TRUE, FALSE, FALSE)))
stopifnot(integer_argument(character(), "from", 1L) == 1L,
          integer_argument("--from=2", "from", 1L) == 2L)
for (args in list("--from=x", "--from=-1", c("--from=1", "--from=2"))) {
  stopifnot(inherits(try(integer_argument(args, "from", 1L), silent = TRUE), "try-error"))
}

definitions <- comparison_definitions()
saved_definitions <- read.csv("results/tables/contrast_definitions.csv")
stopifnot(identical(definitions, saved_definitions))

# Exercise all 13 settings on no separation, quasi/complete separation, and
# a covariate model. Compare numeric output to the saved full experiment.
test_ids <- c("n50_7.18.1.24", "n50_0.25.1.24", "n50_25.0.0.25",
              covmeta$config_id[which(covmeta$data_state == "ok")[1]])
saved <- readRDS("results/model_fits.rds")
saved <- saved[saved$config_id %in% test_ids & saved$engine %in% c("logistf", "brglm2"), ]
for (id in test_ids) {
  meta <- metadata[match(id, metadata$config_id), ]
  dat <- model_data(meta, covariates)
  for (variant in variants) {
    actual <- standard_r(fit_one(meta, variant, dat)$result)
    expected <- saved[saved$config_id == id & saved$variant == variant$id, ]
    expected <- expected[match(actual$inference, expected$inference), ]
    stopifnot(nrow(actual) == nrow(expected), !anyNA(expected$config_id))
    for (column in c("beta", "se", "low", "high", "p", "converged", "ci_converged")) {
      stopifnot(isTRUE(all.equal(actual[[column]], expected[[column]],
                                tolerance = 1e-10, check.attributes = FALSE)))
    }
  }
}
message("All 52 fits agree with saved coefficients, intervals, p-values and convergence flags.")
