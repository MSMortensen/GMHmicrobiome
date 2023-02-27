#' Filter by abundance
#'
#' @description Filter low abundance taxa in a phyloseq object and group them as "Others". Can be run for all samples together or per group
#' @param pobject Phyloseq object to filter
#' @param group Sample data variable to filter by Default NA
#' @param min.abundance cut-off ratio (taxa with mean abundance below this ratio (0.1 = 10 \%) will be group in "Others") Default 0.001
#' @param includes If using groups, should taxa be included if they are above the min.abundance in any, all, or a specific number of groups Default "any"
#' @return A filtered phyloseq object
#' @export

filter_abundance <- function(pobject, group = NA, min.abundance=0.01, includes = "any"){
  loadNamespace("phyloseq")

  # transform to sample counts to relative values
  transformed <- transform_sample_counts(pobject, function(x) x/sum(x))

  # Orientate and export otu table
  if (taxa_are_rows(transformed)) otu.table <- as.data.frame(otu_table(transformed)) else otu.table <- as.data.frame(t(otu_table(transformed)))

  # remove Others from otu table
  otu.table <- otu.table[row.names(otu.table) != "Others",]

  # Calculate abundance
  if (is.na(group)) {
    counts <- data.frame(row.names = row.names(otu.table),
                         all = rowMeans(otu.table))
    counts$keep <- ifelse(counts$all > min.abundance, TRUE, FALSE)
  } else {
    # extract metadata
    dat <- data.frame(sample_data(pobject))

    # set groups
    vgroup <- as.character(unique(dat[,group]))

    # create count table
    counts <- data.frame(matrix(nrow = nrow(otu.table), ncol = length(vgroup)))
    row.names(counts) <- row.names(otu.table)
    colnames(counts) <- vgroup

    for (i in seq(length(vgroup))) counts[,vgroup[i]] <- rowMeans(otu.table[,dat[,group]==vgroup[i]])

    if (includes == "any") {
      counts$keep <- ifelse(rowSums(counts > min.abundance) > 0, TRUE, FALSE)
    } else if (includes == "all") {
      counts$keep <- ifelse(rowSums(counts > min.abundance) == length(vgroup), TRUE, FALSE)
    } else if (is.numeric(includes)) counts$keep <- ifelse(rowSums(counts > min.abundance) >= includes, TRUE, FALSE)
  }

  # Create list of ASVs to merge
  otu.list <- row.names(counts)[!counts$keep]
  if(any("Others" %in% taxa_names(pobject))) otu.list <- c(otu.list,"Others")

  # Output count of filtered taxa
  if(length(otu.list) == length(taxa_names(pobject))) stop("All features would be grouped as 'Others'!")
  if(length(otu.list) == 0) {
    merged <- pobject
    message(paste("No features were grouped!"))
  } else {
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
  }
  return(merged)
}
