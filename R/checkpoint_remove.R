#' Remove micronicheR checkpoint files
#'
#' Deletes a checkpoint directory created by
#' \code{\link{endonicheR}}.
#'
#' @param checkpoint_dir Character. Path to the checkpoint directory.
#'
#' @return Invisibly returns TRUE if the directory was removed.
#'
#' @examples
#' \dontrun{
#'
#' micronicheR::checkpoint_remove(
#'   checkpoint_dir = "checkpoint_test"
#' )
#'
#' }
#'
#' @export

checkpoint_remove <- function(
    checkpoint_dir
){

  if (missing(checkpoint_dir)) {

    stop(
      "'checkpoint_dir' must be provided."
    )

  }

  if (!dir.exists(checkpoint_dir)) {

    stop(
      "Directory does not exist:\n",
      checkpoint_dir
    )

  }

  unlink(
    checkpoint_dir,
    recursive = TRUE,
    force = TRUE
  )

  message(
    "Checkpoint directory removed:\n",
    normalizePath(
      checkpoint_dir,
      winslash = "/",
      mustWork = FALSE
    )
  )

  invisible(TRUE)

}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
