#' Calculate alpha diversity using multiple rarefaction
#'
#' @description
#' @param pobject A \code{phyloseq} object.
#' @param ntables Number or rarefactions to perform. Default 100
#' @param depth Rarefaction depth at which the alpha diversity is calculated. Default 90 \% of the lowest sample depth
#' @param methods Which alpha diversity metrics to calculate. Default c("Observed","Chao1","FaithPD","Shannon")
#' @param seedstart Random seed to use. Default 500
#' @param verbose Logical, indicating if the current sample being rarefied should be reported during calculation. Default FALSE
#'
#' @return Data frame with mean and sd alpha diversity for each alpha diveristy metric added to the sample data from phyloseq.
#' @export
#'
#' @examples
calculate_alpha_diversity <- function(pobject, ntables=100, depth = round(min(sample_sums(pobject))*0.9), methods=c("Observed","Chao1","FaithPD","Shannon"), seedstart=500, verbose=FALSE) {
  require("vegan")
  loadNamespace("phyloseq")

  # remove samples below depth
  phy.use <- prune_samples(sample_sums(pobject) >= depth, pobject )

  # Orientate the OTU correctly
  if (taxa_are_rows(phy.use)){otu.tab<-unclass(t(otu_table(phy.use)))} else otu.tab <- unclass(otu_table(phy.use))

  # Rarefaction function
  rarefy <- function(x, depth) {
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    j <- numeric(length(x))
    j[as.numeric(names(y.tab))] <- y.tab
    j
  }

  # Table to output alpha diversity table
  Alpha_diversity = data.frame(row.names = row.names(otu.tab))

  for (i in seq(length(methods))){
    Alpha_diversity[,methods[i]] <- numeric(length = nrow(otu.tab))
    Alpha_diversity[,paste0(methods[i],"_sd")] <- numeric(length = nrow(otu.tab))
  }

  # Run each sample separately
  for (z in 1:nrow(otu.tab)) {
    if (verbose==TRUE) {
      print(paste("Rarefaction sample number", z, sep=" "))
    }
    numbers <- otu.tab[z,]

    # Rarefy the sample ntables times
    set.seed(seedstart + z)
    rare_tab <- lapply(1:ntables,function(k) rarefy(numbers,depth))

    # Format table
    rare_tab <- do.call(rbind, rare_tab)

    # Calculate Observed richness, Chao1, and ACE.
    adiv <- data.frame(t(estimateR(rare_tab)))

    if ("Observed" %in% methods){
      # Save mean and sd of observed richness
      Alpha_diversity$Observed[z] <- mean(adiv$S.obs)
      Alpha_diversity$Observed_sd[z] <- sd(adiv$S.obs)
    }

    if ("Chao1" %in% methods){
      # Save mean and sd of observed richness
      Alpha_diversity$Chao1[z] <- mean(adiv$S.chao1)
      Alpha_diversity$Chao1_sd[z] <- sd(adiv$S.chao1)
    }

    if ("ACE" %in% methods){
      # Save mean and sd of observed richness
      Alpha_diversity$ACE[z] <- mean(adiv$se.ACE)
      Alpha_diversity$ACE_sd[z] <- sd(adiv$se.ACE)
    }

    if ("Shannon" %in% methods){
      # Calculate observed richness for each rep of sample z
      adiv <- vegan::diversity(rare_tab, index = "shannon")

      # Save mean and sd of observed richness
      Alpha_diversity$Shannon[z] <- mean(adiv)
      Alpha_diversity$Shannon_sd[z] <- sd(adiv)
    }

    if ("Simpson" %in% methods){
      # Calculate observed richness for each rep of sample z
      adiv <- diversity(rare_tab, index = "simpson")
      # Save mean and sd of observed richness
      Alpha_diversity$Simpson[z] <- mean(adiv)
      Alpha_diversity$Simpson_sd[z] <- sd(adiv)
    }

    if ("Evenness" %in% methods){
      # Calculate observed richness for each rep of sample z
      sha <- diversity(rare_tab, index = "shannon")
      obs <- rowSums(rare_tab != 0)
      adiv <- sha/log(obs)
      # Save mean and sd of observed richness
      Alpha_diversity$Evenness[z] <- mean(adiv)
      Alpha_diversity$Evenness_sd[z] <- sd(adiv)
    }

    if ("FaithPD" %in% methods){
      colnames(rare_tab) <- taxa_names(pobject)
      # Calculate Faith Phylogenetic distance for each rep of sample z
      tmp <- pd(rare_tab, phy_tree(pobject), include.root = T)
      Alpha_diversity$FaithPD[z] <- mean(tmp$PD)
      Alpha_diversity$FaithPD_sd[z] <- sd(tmp$PD)
    }

  }

  # Add alpha diversity to sample data
  output <- cbind(sample_data(phy.use),Alpha_diversity)
  output$depth = depth
  # Return pobject to the environment
  return(output)
}
