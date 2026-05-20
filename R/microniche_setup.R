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

  if (!is.null(cran)) {

    missing_cran <- cran[!vapply(
      cran,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]

    if (length(missing_cran) > 0) {
      message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
      utils::install.packages(missing_cran, dependencies = TRUE)
    } else {
      message("All CRAN packages are already installed.")
    }

    failed_cran <- cran[!vapply(
      cran,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]

    if (length(failed_cran) > 0) {
      stop(
        "The following CRAN packages could not be installed or loaded: ",
        paste(failed_cran, collapse = ", ")
      )
    }
  }

  if (!is.null(github)) {

    if (!requireNamespace("remotes", quietly = TRUE)) {
      utils::install.packages("remotes")
    }

    for (repo in github) {

      pkg <- basename(repo)

      if (!requireNamespace(pkg, quietly = TRUE)) {

        message("Installing GitHub package: ", repo)

        remotes::install_github(
          repo,
          dependencies = TRUE,
          upgrade = "never"
        )
      } else {
        message(pkg, " is already installed.")
      }

      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
          "GitHub package '", pkg, "' could not be installed or loaded.\n",
          "Repository: ", repo, "\n\n",
          "On Windows, this may require Rtools. For NicheMapR, install Rtools44, ",
          "restart RStudio, and try again."
        )
      }
    }
  }

  message("All requested packages are ready.")

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
