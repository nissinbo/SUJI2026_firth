collapse_messages <- function(x) {
  x <- gsub("[\r\n]+", "\\\\n", trimws(x))
  x <- unique(x[nzchar(x)])
  if (length(x)) paste(x, collapse = " | ") else ""
}

with_conditions <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(message = conditionMessage(e)), class = "captured_error")
  )
  list(value = value, warnings = warnings, error = inherits(value, "captured_error"))
}

double_to_hex <- function(x) {
  vapply(as.double(x), function(one) {
    raw <- writeBin(one, raw(), size = 8L, endian = "big")
    paste(sprintf("%02X", as.integer(raw)), collapse = "")
  }, character(1))
}

# Strict parsing avoids silently resuming a different range after a typo.
integer_argument <- function(args, name, default) {
  prefix <- paste0("--", name, "=")
  values <- substring(args[startsWith(args, prefix)], nchar(prefix) + 1L)
  if (!length(values)) return(default)
  if (length(values) != 1L || !grepl("^[0-9]+$", values)) {
    stop("Expected one nonnegative integer for --", name)
  }
  value <- suppressWarnings(as.integer(values))
  if (is.na(value)) stop("Integer out of range for --", name)
  value
}
