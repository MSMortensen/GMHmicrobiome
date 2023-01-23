#' Filter by prevalence
#'
#' @description Filter low prevalence taxa in a phyloseq object and group them as "Others". Can be run for all samples or by groups
#' @param pobject A \code{phyloseq} object to filter
#' @param group Sample data variable to filter by Default NA
#' @param min.samples cut-off value per group (taxa below this value will be group in "Others") Default 10% of samples
#' @param includes If using groups, should taxa be if they are prevalent in more than min.samples in any, all, or a specific number of groups Default "all"
#'
#' @return A \code{phyloseq} object with the excluded taxa grouped as "Others"
#' @export
#'
#' @examples
filter_prevalence <- function(pobject, group = NA, min.samples = (length(sample_names(pobject))/10), includes = "any"){
  loadNamespace("phyloseq")

  # # Orientate and export otu table
  if (taxa_are_rows(pobject)) otu.table <- as.data.frame(otu_table(pobject)) else otu.table <- as.data.frame(t(otu_table(pobject)))

  # make binary
  otu.table[otu.table != 0] <- 1

  # Calculate prevalence
  if (is.na(group)) {
    counts <- data.frame(row.names = row.names(otu.table),
                         all = rowSums(otu.table))
    counts$keep <- ifelse(counts$all > min.samples, TRUE, FALSE)
  } else {
    # extract metadata
    dat <- data.frame(sample_data(pobject))

    # set groups
    vgroup <- as.character(unique(dat[,group]))

    # create count table
    counts <- data.frame(row.names = row.names(otu.table))
    for (i in seq(length(vgroup))) counts[,vgroup[i]] <- rowSums(otu.table[,dat[,group]==vgroup[i]])

    if (includes == "any") {
      counts$keep <- ifelse(rowSums(counts > min.samples) > 0, TRUE, FALSE)
    } else if (includes == "all") {
      counts$keep <- ifelse(rowSums(counts > min.samples) == length(vgroup), TRUE, FALSE)
    } else if (is.numeric(includes)) counts$keep <- ifelse(rowSums(counts > min.samples) >= includes, TRUE, FALSE)
  }

  # Create list of ASVs to merge
  otu.list <- row.names(counts[!counts$keep,])
  if(any("Others" %in% rownames(otu.table))) otu.list <- c(otu.list,"Others")

  # Output count of filtered taxa
  if(length(otu.list) == 0) stop("No features to group!")
  if(length(otu.list) == nrow(count_table)) stop("All features would be grouped as 'Others'!")
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
