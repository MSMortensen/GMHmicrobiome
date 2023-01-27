#' Sort taxa by full taxonomy
#'
#' @description Sorts a data.frame based on all available taxonomic levels, and ensures that the factor levels of the last reflects this order
#' @param dat Dataframe
#'
#' @return Dataframe ordered by taxonomy and where the lowest taxonomic level is a factor with levels in the present order.
#' @export
#'
#' @examples
sort_taxa <- function(dat) {
  # set columns
  tax.levels <- c("Kingdom", "Phylum", "Class","Order","Family","Genus","Species", "Taxa")
  tax.use <- tax.levels[tax.levels %in% colnames(dat)]

  # Sort data
  dat.ordered <- dat %>% arrange_at(tax.use)

  dat.ordered <- as.data.frame(dat.ordered)

  # extract sorted levels
  new.levels <- as.character(unique(dat.ordered[,last(tax.use)]))


  dat.ordered[,last(tax.use)] <- factor(dat.ordered[,last(tax.use)], levels = new.levels)

  return(dat.ordered)
}
