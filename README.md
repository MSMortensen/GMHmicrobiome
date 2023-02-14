GMH microbiome analalysis package
================

## Intro

This package is build for support analyses of microbiome data at [The
Research Group for Gut, Microbes, and
Health](https://www.food.dtu.dk/english/research/gut-microbes-and-health "Group page")
at [DTU National Food
Institute](https://www.food.dtu.dk/english "Institute page"), Denmark.

The package consists of templates for data import, microbiome
description, analysis of individual variables, beta diversity analysis,
and differential abundance testing. Additionally, some custom functions
used in the templates are included.

## Content

The package is structured with Rmarkdown templates that together makes
up a fundamental microbiome analysis as well as some functions that
simplify parts of the analyses.

### Templates

As the package is still a work in progress, this is a complete list of
all templates that will be included. Not all of the templates have been
finalised, but this will be updated when they are added.

- **GMH_import**  
  *Input*: DADA2 generated phyloseq object and a metadata csv file.  
  *Description*: Import, cleaning, QC, decontamination, and calculates
  alpha diversities.  
  *Output*: Two phyloseq objects (with different levels of
  decontamination) and a structured phyloseq object summarizing the
  process.

- **GMH_description**  
  *Input*: One phyloseq objects from **GMH_import.**  
  *Description*: Generates a general description of the microbiomes in
  the samples. This does not include any statistical analysis.  
  *Output*: A structured .html document with all tables and plots.

- **GMH_test_variables** (Lacks visualization part)  
  *Input*: One phyloseq objects from **GMH_import.**  
  *Description*: Guide to identify and apply the correct statistical
  test differences and correlations between alpha diversity and project
  variables.  
  *Output*: Relevant plots and a structured .html document showing the
  statistical analyses.

- **GMH_beta_diversity**  
  *Input*: One phyloseq objects from **GMH_import.**  
  *Description*: Fundamental analysis of beta diversity, identifying
  differences between groups and correlations with quantitative
  variables.  
  *Output*: Relevant plots and a structured .html document showing the
  beta diversity analysis

- **GMH_differential_abundance** (In progress, only placeholder
  template)  
  *Input*: One phyloseq objects from **GMH_import.**  
  *Description*: Identification of differentially expressed taxa. This
  is structured around the DAtest package and supplemented with plots.  
  *Output*: Relevant plots and a structured .html document showing the
  differential abundance analysis.

### Functions

The package contains the following functions:

- **Project initiation**  
  `initiate_GMH_project()`: Creates the folder structure used by the
  script in the the current working directory. It can also create copies
  of the five templates with generic names.

- **Rarefaction and alpha diversity**  
  `calculate_alpha_diversity()`: Calculates alpha diversity as the mean
  of multiple independent rarefactions.  
  `Rcurve_data()`: Wraps around `calculate_alpha_diversity()` to
  calculate data for a rarefaction curve plot.  
  `multiple_rarefy()`: Rarefies each sample multiple times and select
  the best representation for each sample based on within sample
  betadiversity comparison.

- **Filtering**  
  `filter_abundance()`: Filters taxa based on average abundance and
  merges low abundant taxa into “Others”.  
  `filter_prevalence()`: Filters taxa based on prevalence
  (presence/absence) and merges low prevalence taxa into “Others”.  
  `filter_rank()`: Filters taxa based on ranked average abundance and
  merges lower ranked taxa into “Others”.  
  `filter_variance()`: Filters taxa based on coefficient of variation
  (cov) and merges less variating taxa into “Others”. Low abundant taxa
  more likely to have high cov so I recommend filtering by abundance or
  prevalence first.

- **Taxa cleaning**  
  `clean_taxa()`: Replaces NA values in `tax_table()` with the lowest
  identified classification for that taxa.  
  `sort_taxa()`: Sort taxa by full taxonomy for plotting (e.i. taxa from
  the same phylum is plottet next to each other).

## Installation

``` r
if (!requireNamespace("remotes")) install.packages("remotes")

remotes::install_github("MSMortensen/GMHmicrobiome")
```

## How to use

I recommend starting with a clean RStudio project and then follow the
steps below:

- Run `initiate_GMH_project(files = FALSE)` to create the necessary
  folders.

- Copy input data to the input folder

- Import data

  - Create a new Rmarkdown file from the **GMH_import** template  
    (File \> New File \> R Markdown … \> From Template \> GMH_import).

  - Modify the template as necessary and when all works, knit the
    Rmarkdown.

- Describe the microbiome

  - Create a new Rmarkdown file from the **GMH_import** template  
    (File \> New File \> R Markdown … \> From Template \>
    GMH_description).

  - Modify the template as necessary and when all works, knit the
    Rmarkdown.

- Run any of the three analysis templates.
