# Shared project-level helpers for executable analysis and presentation scripts.

set_project_locale <- function() {
  if (.Platform$OS.type != "windows") return(invisible(NULL))
  locale <- suppressWarnings(Sys.setlocale("LC_CTYPE", "Japanese_Japan.utf8"))
  if (!nzchar(locale)) {
    warning("Could not set Japanese UTF-8 locale; continuing with the current locale.", call. = FALSE)
  }
  invisible(locale)
}

ensure_directories <- function(paths) {
  for (path in paths) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(paths)
}

assert_project_root <- function(marker = "presentation/slides.qmd") {
  if (!file.exists(marker)) {
    stop("Run from the project root: missing ", marker, call. = FALSE)
  }
  invisible(TRUE)
}
