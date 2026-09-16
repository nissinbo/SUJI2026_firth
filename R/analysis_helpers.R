# Shared output validity rules. Unknown convergence is never treated as success.
ci_ok <- function(z) {
  is.finite(z$low) & is.finite(z$high) & z$low <= z$high
}

p_ok <- function(z) {
  is.finite(z$p) & z$p >= 0 & z$p <= 1
}

confirmed_convergence <- function(flag) {
  !is.na(flag) & flag
}

# Preserve the existing quantile conventions: interpolated quantiles for
# enumerated tables, inverse empirical CDF for simulation summaries.
summary_quantile <- function(x, w, p, kind) {
  if (kind == "exact") {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    return(as.numeric(quantile(x, p, type = 7, names = FALSE)))
  }
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  x <- x[keep]
  w <- w[keep]
  order <- order(x)
  x <- x[order]
  w <- w[order]
  x[which(cumsum(w) >= p * sum(w))[1]]
}
