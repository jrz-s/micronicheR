#' Download and register NicheMapR global climate data
#'
#' @title Set up NicheMapR global climate data
#' @description Downloads and registers the global climate data required by `NicheMapR::micro_global()`.
#' @author
#' \itemize{
#'   \item Zárate-Salazar, J. Rafael, PhD
#' }
#' @param folder Character. Folder where the global climate data should be downloaded. Default is the current working directory.
#' @param timeout Numeric. Download timeout in seconds. Default is 600.
#' @return Invisibly returns TRUE if the setup process finishes.
#' @export
microniche_setup_global_climate <- function(
    folder = getwd()
    , timeout = 600
) {

  if (!requireNamespace("NicheMapR", quietly = TRUE)) {
    stop(
      "Package 'NicheMapR' is required.\n",
      "Install it first using:\n",
      "microniche_setup(github = 'mrke/NicheMapR')"
    )
  }

  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  message(
    "Downloading and registering NicheMapR global climate data in:\n",
    folder
  )

  NicheMapR::get.global.climate(folder = folder)

  message(
    "\nGlobal climate data setup finished.\n",
    "Please restart R before running endonicheR() for the first time."
  )

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
