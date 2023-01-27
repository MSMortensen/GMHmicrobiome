#' Annotation ZymoBIOMICS™ Microbial Community DNA Standard
#'
#' The ZymoBIOMICS™ Microbial Community DNA Standard (Catalog Nos. D6305 (200ng) and D6306 (2000ng)) genomes processed through the DF_GMH_pipeline to determine the expected abundance of the mock community when used in the DF_GMH_pipeline. Begining from full genomes, the 16S sequences were extracted with Barnap, and then trimmed with cutadapt (standard pipeline settings). Following annotation data was aggregated to Species level. It is important to note that this is an ~~in silico~~ representation of the mock community annotated with our pipeline, so this will differ from the actual content.
#'
#' @format ## `who`
#' A data frame with 6 rows and 9 columns:
#' \describe{
#'   \item{Kingdom}{Assigned Kingdom}
#'   \item{Phylum}{Assigned Phylum}
#'   \item{Class}{Assigned Class}
#'   \item{Order}{Assigned Order}
#'   \item{Family}{Assigned Family}
#'   \item{Genus}{Assigned Genus}
#'   \item{Species}{Assigned Species}
#'   \item{SampleID}{Sample name for merger with actual samples}
#'   \item{Abundance}{Calculated relative abundance}
#'   ...
#' }
#' @source https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip>
"who"
