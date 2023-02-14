#' Multiple rarefication
#'
#' @description Create multiple rarefications of each sample and select most representative iteration.
#' @param pobject A \code{phyloseq} object.
#' @param ntables Number or rarefactions to perform. Default 100
#' @param depth Rarefaction depth at which the alpha diversity is calculated. Default 90 \% of the lowest sample depth
#' @param distmethod Distance/dissimilarity metric used to asses iteration similarity. Must be one of: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison". Default "bray"
#' @param summarymeasure Function to compare iterations. Default mean
#' @param seedstart Random seed to use. Default 500
#' @param verbose Logical, indicating if the current sample being rarefied should be reported during calculation. Default TRUE
#'
#' @return A rarefied \code{phyloseq} object with rarefied samples
#' @export

multiple_rarefy <- function(pobject, ntables=100, depth = round(min(sample_sums(pobject))*0.9), distmethod="bray", summarymeasure=mean, seedstart=500, verbose=TRUE) {
  require("vegan")
  loadNamespace("phyloseq")

  # Orientate the OTU correctly
  if (taxa_are_rows(pobject)){rawtab<-unclass(t(otu_table(pobject)))} else rawtab <- unclass(otu_table(pobject))

  # Ignore samples below rarefaction depth
  ind <- (rowSums(rawtab) < depth)
  sam.discard <- rownames(rawtab)[ind]
  otu.tab <- rawtab[!ind, ]

  # Rarefaction function
  rarefy <- function(x, depth) {
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    j <- numeric(length(x))
    j[as.numeric(names(y.tab))] <- y.tab
    j
  }

  # Table to output rarefied data
  final_tab = c()

  # Run each sample separately
  for (z in 1:nrow(otu.tab)) {
    if (verbose==TRUE) {
      print(paste("Rarefaction sample number", z, sep=" "))
    }
    numbers <- otu.tab[z,]

    # Rarefy the sample ntables times
    set.seed(seedstart + z)
    rare_tab <- lapply(1:ntables,function(k) rarefy(numbers,depth))

    rare_tab <- do.call(rbind, rare_tab)
    # # Remove columns with no reads
    # rare_tab_no_zero <- rare_tab[,colSums(rare_tab) != 0]
    # # distance across reps for subject z
    distmat = as.matrix(vegdist(rare_tab, method=distmethod))
    # calculate mean distance for each rep
    distsummary = apply(distmat, 2, summarymeasure)
    # the best rep is the one with the mean distance to all other reps. (in case of ties, just select the first)
    whichbestrep = which(distsummary == min(distsummary))[1]
    # select that rep only for subject z
    bestrep = rare_tab[whichbestrep,]
    # build that rep for subject y into final table
    final_tab = rbind(final_tab, bestrep)
  }

  # Remove samples with too few reads
  pobject <- prune_samples(!sample_names(pobject) %in% sam.discard, pobject)
  # Reformat final tab and return to the pobject object
  rownames(final_tab) = rownames(otu.tab)
  colnames(final_tab) = colnames(otu.tab)
  otu_table(pobject) <- otu_table(t(final_tab), taxa_are_rows = T)

  # Return pobject to the environment
  return(pobject)
}
