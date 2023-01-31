#' Initiate GMH project
#'
#' @description Creates the folders and/or import the base Rmarkdown templates needed to run basic microbiome analysis using the \code{GMHmicrobiome} pipeline.
#' @param folders Logical indicator if folder should be created. Default = TRUE
#' @param files Logical indicator if Rmarkdown templates should be created. Default = TRUE
#' @param overwrite Logical indicator if files should be overwritten if present. Default = FALSE
#' @return
#' @export
#'
#' @examples
initiate_GMH_project <- function(folders = TRUE, files = TRUE, overwrite = FALSE){
  # Create used folders if missing
  if (isTRUE(folders)) {
    if (!file.exists("R_objects")) dir.create(file.path(getwd(), "R_objects"))
    if (!file.exists("plots")) dir.create(file.path(getwd(), "plots"))
    if (!file.exists("tables")) dir.create(file.path(getwd(), "tables"))
    if (!file.exists("output")) dir.create(file.path(getwd(), "output"))
  }

  # Create the initial files if not found
  if (isTRUE(files) & isFALSE(overwrite)) {
    file.copy(system.file("rmarkdown/templates/gmh_import/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_1_import.Rmd")
    file.copy(system.file("rmarkdown/templates/gmh_description/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_2_description.Rmd")
    file.copy(system.file("rmarkdown/templates/gmh_test_variables/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_3_test_variables.Rmd")
    file.copy(system.file("rmarkdown/templates/gmh_beta_diversity/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_4_betadiversity.Rmd")
    file.copy(system.file("rmarkdown/templates/gmh_differential_abundance/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_5_differential_abundance.Rmd")
  }

  # Create and/or overwrite the initial files
  if (isTRUE(files) & isFALSE(overwrite)) {
    if (!file.exists("GMH_1_import.Rmd")) file.copy(system.file("rmarkdown/templates/gmh_import/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_1_import.Rmd")
    if (!file.exists("GMH_2_description.Rmd")) file.copy(system.file("rmarkdown/templates/gmh_description/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_2_description.Rmd")
    if (!file.exists("GMH_3_test_variables.Rmd")) file.copy(system.file("rmarkdown/templates/gmh_test_variables/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_3_test_variables.Rmd")
    if (!file.exists("GMH_4_betadiversity.Rmd")) file.copy(system.file("rmarkdown/templates/gmh_betadiversity/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_4_betadiversity.Rmd")
    if (!file.exists("GMH_5_differential_abundance.Rmd")) file.copy(system.file("rmarkdown/templates/gmh_differential_abundance/skeleton","skeleton.Rmd", package = "GMHmicrobiome"), "GMH_5_differential_abundance.Rmd")
  }
}

