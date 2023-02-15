#' Calculate rarefaction curve data with multiple rarefactions
#'
#' @description Wrapper around calculate_alpha_diversity to generate rarefaction curve data with multiple rarefactions of each sample
#' @param pobject A \code{phyloseq} object.
#' @param ntables Number or rarefactions to perform. Default 10
#' @param step Resolution of sequencing depth (interval between rarefactions). Default 250]
#' @param maxdepth Maximum sequencing depth to include Default sample depth of 90th percentile
#' @param methods Which alpha diversity metrics to calculate. Default c("Observed","Chao1","FaithPD","Shannon")
#' @param seedstart Random seed to use. Default 500
#' @param verbose Logical, indicating if information should be reported during calculation. Default FALSE
#'
#' @return Data frame with average alpha diversity per sample at each sequencing depth (from 1 to maxdepth)
#' @export

Rcurve_data <- function(pobject, ntables=10, step=250,maxdepth = round(unname(quantile(sample_sums(pobject),0.9))), methods=c("Observed","Chao1","FaithPD","Shannon"), seedstart=500, verbose=FALSE) {
  require("vegan")
  loadNamespace("phyloseq")

  # prep list of
  step.seq <- seq(from = 1, to = maxdepth, by = step)

  # Calculate alpha diversity
  rare_tab <- lapply(step.seq,function(k) calculate_alpha_diversity(pobject = pobject, ntables = ntables, depth = k, INDECES = methods, seedstart = seedstart, verbose = verbose))

  # Format table
  rare_tab <- do.call(rbind, rare_tab)

  return(rare_tab)
}
