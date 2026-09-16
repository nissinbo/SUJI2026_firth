# Run from the project root. --rebuild also re-imports the full raw model outputs.
# Model fitting and input generation are separate, explicit operations.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Japanese_Japan.utf8")
if (!file.exists("presentation/slides.qmd")) stop("Run from the project root")
for (folder in c("results/tables", "results/figures", "results/checks")) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
steps <- c("summarize", "validate_results", "check_wald_se", "summarize_pl_alignment",
           "plot_results", "plot_heatmaps", "coverage", "plot_intro_bias")
for (step in steps) {
  message("Running R/", step, ".R")
  source(paste0("R/", step, ".R"), local = new.env(parent = globalenv()), encoding = "UTF-8")
  gc()
}
writeLines(capture.output(sessionInfo()), "results/checks/session.txt")
message("Updated and checked all presentation tables and figures.")
