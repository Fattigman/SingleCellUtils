#' Normalize Single-Cell Data
#'
#' This function normalizes single-cell RNA-seq data using a simple method.
#'
#' @param counts A matrix of raw counts.
#' @param method A string specifying the normalization method ("log" or "size_factor").
#' @return A matrix of normalized counts.
#' @export
normalize_counts <- function(counts, method = "log") {
  if (method == "log") {
    return(log1p(counts))
  } else if (method == "size_factor") {
    size_factors <- colSums(counts)
    return(t(t(counts) / size_factors) * mean(size_factors))
  } else {
    stop("Unknown normalization method")
  }
}
