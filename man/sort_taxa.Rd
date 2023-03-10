% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sort_taxa.R
\name{sort_taxa}
\alias{sort_taxa}
\title{Sort taxa by full taxonomy}
\usage{
sort_taxa(dat)
}
\arguments{
\item{dat}{Dataframe}
}
\value{
Dataframe ordered by taxonomy and where the lowest taxonomic level is a factor with levels in the present order.
}
\description{
Sorts a data.frame based on all available taxonomic levels, and ensures that the factor levels of the last reflects this order
}
