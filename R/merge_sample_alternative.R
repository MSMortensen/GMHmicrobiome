#' Merge samples alternative
#'
#' @description Modified version of the phyloseq merge samples function where it is now possible to use the mean instead of sum for otu table.
#' @param x A \code{phyloseq} object to filter
#' @param group The variable that should be used to merge the samples
#' @param fun Function to be used for merging otu_table. Default mean
#'
#' @return The provided \code{phyloseq} object with an updated otu_table() and sample_data()
#' @export

  merge_samples_alternative <- function(x, group, fun=mean){

  # Check if phyloseq object has a sample_data
    # Check class of group and modify if single "character" (column name)
    if( class(group)=="character" & length(group)==1 ){
      x1 <- data.frame(sample_data(x))
      if( !group %in% colnames(x1) ){stop("group not found among sample variable names.")}
      group <- x1[, group]
    }
    # coerce to factor
    if( class(group)!="factor" ){ group <- factor(group) }
    # Perform merges.
    ## Sample data
    x1    <- data.frame(sample_data(x))

    # Remove any non-coercable columns.
    # Philosophy is to keep as much as possible. If it is coercable at all, keep.
    # Coerce all columns to numeric matrix
    coercable    <- sapply(x1, canCoerce, "numeric")
    x2           <- sapply(x1[, coercable], as, "numeric")
    rownames(x2) <- rownames(x1)

    # Perform the aggregation.
    outdf <- aggregate(x2, list(group), fun)
    # get rownames from the "group" column (always first)
    # rownames(outdf) <- as.character(outdf[, 1])
    rownames(outdf) <- levels(group)
    # "pop" the first column
    newSM <- outdf[, -1, drop=FALSE]

    ## OTU table
    if( taxa_are_rows(x) ){ x1 <- t(otu_table(x)) } else { x1 <- otu_table(x)}
    # coerce to matrix, x2
    x2 <- as(x1, "matrix")

    # Aggregate samples
    out <- as.matrix(data.frame(stats::aggregate(x2, list(group), fun),row.names = 1))

    # convert back to otu_table, and return
    newOT <- otu_table(out, taxa_are_rows=FALSE)

    phyloseqList <- list(newOT, newSM)

  ### Add to build-call-list the remaining components, if present in x.
  ### NULL is returned by accessor if object lacks requested component/slot.
  ### Order of objects in list doesn't matter for phyloseq.
  ### The list should not be named.
  if( !is.null(access(x, "tax_table")) ){ phyloseqList <- c(phyloseqList, list(tax_table(x))) }
  if( !is.null(access(x, "phy_tree"))    ){ phyloseqList <- c(phyloseqList, list(phy_tree(x))) }

  return( do.call("phyloseq", phyloseqList) )
}
