#' Filter by variance
#'
#' @description Filter taxa with low coefficient of variation in a phyloseq object and group them as "Others". Can be run for all samples or by groups
#' @param pobject A \code{phyloseq} object to filter
#' @param taxa.lim Numeric, if < 1 indicates the ratio of all taxa to keep and if > 1 indicates the number of taxa to keep
#' @return A \code{phyloseq} object with the excluded taxa grouped as "Others"
#' @export

filter_variance <- function(pobject, taxa.lim = 0.9){
  loadNamespace("phyloseq")

  # Normalise to ratios
  p.rel <- transform_sample_counts(pobject, function(x) x/sum(x))

  # # Orientate and export otu table
  if (taxa_are_rows(p.rel)) otu.table <- as.data.frame(otu_table(p.rel)) else otu.table <- as.data.frame(t(otu_table(p.rel)))

  # Update number of taxa to keep
  filt.tax <- ifelse(taxa.lim < 1, round(nrow(otu.table)*taxa.lim), taxa.lim)

  # Calculate cov and rank
  columns <- colnames(otu.table)
  sum.dat <- otu.table %>%
    mutate(ASV = row.names(otu.table)) %>%
    rowwise() %>%
    transmute(ASV = ASV,
              cov = abs(sd(c_across(all_of(columns)))/mean(c_across(all_of(columns))))) %>%
    ungroup() %>%
    mutate(Rank = rank(-cov, ties.method = "random"))

  # Create list of ASVs to merge
  otu.list <- sum.dat %>%
    filter(Rank > filt.tax) %>%
    pull(ASV)

  if(any("Others" %in% taxa_names(pobject))) otu.list <- c(otu.list,"Others")

  # Output count of filtered taxa
  message(paste(length(otu.list),"features grouped as 'Others' in the output"))

  # perform the merger with the original sample counts
  merged <- merge_taxa(pobject, otu.list, 1)

  # change the taxa name
  taxa_names(merged)[taxa_names(merged) %in% otu.list[1]] <- "Others"

  # change taxa to "Other" for all levels not being NAs
  for (i in 1:ncol(tax_table(merged))){
    if (sum(!is.na(tax_table(merged)[,i]))){
      tax_table(merged)["Others",i] <- "Others"
    }
  }
  return(merged)
}
