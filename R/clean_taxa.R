#' Clean phyloseq tax table
#'
#' @param pobject A \code{phyloseq} object to filter
#' @param tax_remove The required level of classification. Any taxa not classified at least to this level will be removed [DEFAULTS = "Phylum"]
#' @param verbose Logical, indicating if summary of removed taxa should be reported. [DEFAULT = TRUE]
#'
#' @return The provided \code{phyloseq} object with an updated tax_table() and poorly classified taxa removed
#' @export
#'
#' @examples
clean_taxa <- function(pobject, tax_remove = "Phylum", verbose = TRUE) {
  loadNamespace("phyloseq")

  # Extract tax table
  tax <- data.frame(tax_table(pobject))

  # list ASVs that should be removed
  remove <- is.na(tax[,tax_remove])

  # remove ASVs
  phy.out <- prune_taxa(!remove, pobject)

  # Calculate and print statistics
  if (verbose) {
    # Calculate sample sums of original and cleaned
    output <- data.frame(row.names = sample_names(pobject),
                         org = sample_sums(pobject),
                         cleaned = sample_sums(phy.out))
    output$removed <- output$org - output$cleaned
    output$prc_removed <- output$removed*100/output$org

    # Print output
    cat("OVERVIEW OF ASVs REMOVED:\n",
        "Removed ASVs (%):\t",
        sum(remove),
        " (",
        round(sum(remove)*100/nrow(tax), digits = 3),
        ")\n",
        "Removed reads (%):\t",
        sum(output$removed),
        " (",
        round(sum(output$removed)*100/sum(output$org), digits = 3),
        ")\n",
        "Mean abundance removed:\t",
        round(mean(output$prc_removed), digits = 3),"\n",
        "Max abundance removed:\t",
        round(max(output$prc_removed), digits = 3),"\n", sep = "")
  }

  # Remove NA from tax table
  tax <- data.frame(tax_table(phy.out))

  for (i in seq(nrow(tax))) {
    if (is.na(tax[i,1])) {tax[i,1:7] <- "Unknown"
    } else if (is.na(tax[i,2])) {tax[i,2:7] <- paste(colnames(tax)[1],tax[i,1], sep = "_")
    } else if (is.na(tax[i,3])) {tax[i,3:7] <- paste(colnames(tax)[2],tax[i,2], sep = "_")
    } else if (is.na(tax[i,4])) {tax[i,4:7] <- paste(colnames(tax)[3],tax[i,3], sep = "_")
    } else if (is.na(tax[i,5])) {tax[i,5:7] <- paste(colnames(tax)[4],tax[i,4], sep = "_")
    } else if (is.na(tax[i,6])) {tax[i,6:7] <- paste(colnames(tax)[5],tax[i,5], sep = "_")
    } else if (is.na(tax[i,7])) {tax[i,7] <- paste(colnames(tax)[6],tax[i,6], sep = "_")
    }
  }

  # Insert modified tax_table in phyloseq object
  tax_table(phy.out) <- as.matrix(tax)

  # return the clean phyloseq object
  return(phy.out)
}
