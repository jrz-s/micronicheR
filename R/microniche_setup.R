#' Install and load required packages for microniche
#'
#' @title Install and load microniche dependencies
#' @description Checks, installs, and loads CRAN and GitHub packages required by the microniche package workflow. This function is especially useful for automatically installing `NicheMapR` from GitHub when it is not available locally.
#' @author
#' \itemize{
#'   \item Zárate-Salazar, J. Rafael, PhD
#' }
#' @param cran Character vector. Names of CRAN packages to install and load.
#' @param github Character vector. GitHub repositories in the format `"user/repository"` to install and load.
#' @return Invisibly returns a list of loaded packages.
#' @examples
#' \dontrun{
#' microniche_setup(
#'   cran = c("terra", "dplyr", "tidyr"),
#'   github = c("mrke/NicheMapR")
#' )
#' }
#' @importFrom utils install.packages installed.packages
#' @importFrom remotes install_github
#' @export
microniche_setup <- function(cran = NULL, github = NULL) {

  message("Checking CRAN packages...")

  # Identify missing CRAN packages
  missing_cran <- cran[!cran %in% rownames(installed.packages())]

  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, dependencies = TRUE)
  } else {
    message("All CRAN packages already installed.")
  }

  # Install GitHub packages if provided
  if (!is.null(github)) {

    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }

    for (repo in github) {

      pkg <- basename(repo)

      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing GitHub package: ", repo)
        remotes::install_github(repo, dependencies = TRUE)
      } else {
        message(pkg, " already installed.")
      }
    }
  }

  # Combine package names
  all_packages <- c(cran, basename(github))

  message("Loading packages...")

  invisible(
    lapply(
      all_packages,
      function(pkg) requireNamespace(pkg, quietly = TRUE)
    )
  )

  message("All packages ready.")

}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
