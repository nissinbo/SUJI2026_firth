# Optional PDF quality checks. Requires pdftools; does not modify the PDF.
args <- commandArgs(trailingOnly = TRUE)
check_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
check_root <- if (length(check_args)) dirname(dirname(normalizePath(sub("^--file=", "", check_args[1]), winslash = "/"))) else normalizePath(".", winslash = "/")
if (!requireNamespace("pdftools", quietly = TRUE)) {
  # Optional existing library used by this workstation; no package installation.
  prior_library <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5")
  if (dir.exists(prior_library)) .libPaths(c(.libPaths(), prior_library))
}
if (!requireNamespace("pdftools", quietly = TRUE)) stop("Install pdftools for optional PDF checks")
pdf <- file.path(check_root, "presentation/slides.pdf")
qmd <- readLines(file.path(check_root, "presentation/slides.qmd"), encoding = "UTF-8")
info <- pdftools::pdf_info(pdf)
# R comments inside executable chunks are not slide headings.
in_code <- cumsum(grepl("^```", qmd)) %% 2L == 1L
expected <- sum(grepl("^#{1,2} ", qmd) & !in_code) + 1L
cat("Expected pages:", expected, "; PDF pages:", info$pages, "\n")
count_ok <- info$pages == expected
sizes <- pdftools::pdf_pagesize(pdf)
stopifnot(all(abs(sizes$width / sizes$height - 16 / 9) < 1e-5))
fonts <- pdftools::pdf_fonts(pdf)
stopifnot(all(fonts$embedded))
out <- file.path(check_root, "tmp/slides_qa")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
texts <- pdftools::pdf_text(pdf)
if (!requireNamespace("yaml", quietly = TRUE)) stop("Install yaml for optional PDF checks")
frontmatter_end <- which(qmd == "---")[2]
metadata <- yaml::yaml.load(paste(qmd[2:(frontmatter_end - 1)], collapse = "\n"))
project_metadata <- yaml::read_yaml(file.path(check_root, "_quarto.yml"))
footer <- if ("footer" %in% names(metadata)) metadata$footer else project_metadata$footer
if (!is.null(footer) && !identical(footer, FALSE) && nzchar(as.character(footer))) {
  footer_label <- gsub("[[:space:]]", "", as.character(footer))
  stopifnot(all(grepl(footer_label, gsub("[[:space:]]", "", texts), fixed = TRUE)))
  cat("Configured YAML footer found on every page.\n")
} else {
  cat("No YAML footer configured; footer check skipped.\n")
}
writeLines(unlist(lapply(seq_along(texts), function(i) c(paste0("PAGE ", i), texts[i]))), file.path(out, "extracted.txt"), useBytes = TRUE)
write.csv(fonts, file.path(out, "fonts.csv"), row.names = FALSE)
data <- pdftools::pdf_data(pdf)
for (i in seq_along(data)) {
  z <- data[[i]]
  if (any(z$x < -1 | z$y < -1 | z$x + z$width > sizes$width[i] + 1 | z$y + z$height > sizes$height[i] + 1)) {
    stop(paste("Text outside page", i))
  }
}
if (!"--no-render" %in% args) {
  invisible(pdftools::pdf_convert(pdf, dpi = 110,
    filenames = file.path(out, sprintf("page-%02d.png", seq_len(info$pages))), verbose = FALSE))
}
if (!count_ok) stop("Page count differs from slide count; inspect rendered pages for overflow")
cat("Page count, aspect ratio, embedded fonts, and page bounds passed.\n")
