#' Filter taxa by rank
#'
#' @description Filter taxa by rank and group low rank taxa as "Others". Can be run for all samples together or per group
#' @param pobject Phyloseq object to filter
#' @param group Sample data variable to filter by Default NA
#' @param min.rank cut-off ratio (taxa with mean abundance below this ratio (0.1 = 10 \%) will be group in "Others") Default 10
#' @return A filtered phyloseq object
#' @export

filter_rank <- function(pobject, group = NA, min.rank = 10, includes = "any"){
  loadNamespace("phyloseq")

  # transform to sample counts to relative values
  transformed <- transform_sample_counts(pobject, function(x) x/sum(x))

  # Orientate and export otu table
  if (taxa_are_rows(transformed)) otu.table <- as.data.frame(otu_table(transformed)) else otu.table <- as.data.frame(t(otu_table(transformed)))

  # remove Others from otu table
  otu.table <- otu.table[row.names(otu.table) != "Others",]

  # Calculate rank
  if (is.na(group)) {
    counts <- data.frame(row.names = row.names(otu.table),
                         all = order(rowMeans(otu.table), decreasing = TRUE))
    counts$keep <- ifelse(counts$all <= min.rank, TRUE, FALSE)
  } else {
    # extract metadata
    dat <- data.frame(sample_data(pobject))

    # set groups
    vgroup <- as.character(unique(dat[,group]))

    # create count table
    counts <- data.frame(row.names = row.names(otu.table))
    for (i in seq(length(vgroup))) counts[,vgroup[i]] <- order(rowMeans(otu.table[,dat[,group]==vgroup[i]]), decreasing = TRUE)

    if (includes == "any") {
      counts$keep <- ifelse(rowSums(counts <= min.rank) > 0, TRUE, FALSE)
    } else if (includes == "all") {
      counts$keep <- ifelse(rowSums(counts <= min.rank) == length(vgroup), TRUE, FALSE)
    } else if (is.numeric(includes)) counts$keep <- ifelse(rowSums(counts <= min.rank) >= includes, TRUE, FALSE)
  }
  # Create list of ASVs to merge
  otu.list <- row.names(counts)[!counts$keep]
  if(any("Others" %in% taxa_names(pobject))) otu.list <- c(otu.list,"Others")

  # Output count of filtered taxa
  if(length(otu.list) == length(taxa_names(pobject))) stop("All features would be grouped as 'Others'!")
  if(length(otu.list) == 0) {
    merged <- pobject
    message(paste("Only",length(taxa_names(transformed)), "in input. No features were grouped!", sep = " "))
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
