#' Convert rasters or coordinates into a pixel-level data frame
#' @title Convert raster inputs to data frame
#' @description Converts vegetation cover and elevation rasters into a pixel-level data frame using either a study raster or a set of coordinates. The function aligns, projects, resamples, crops, masks, and extracts raster values as needed.
#' @author
#' \itemize{
#'   \item Zárate-Salazar, J. Rafael, PhD
#' }
#' @param rcover SpatRaster. Raster containing vegetation cover or shade values.
#' @param rtop SpatRaster. Raster containing elevation or topographic values.
#' @param rast_or_coord SpatRaster, matrix, or data.frame. Study raster used to crop and align environmental rasters, or coordinates used to extract raster values. If coordinates are supplied, the first two columns must represent longitude (`x`) and latitude (`y`).
#' @param method_cover Character. Resampling or extraction method for the cover raster. Default is `"near"`.
#' @param method_top Character. Resampling or extraction method for the elevation raster. Default is `"bilinear"`.
#' @param return Character. Output format. Options are `"data"`, `"rasters"`, or `"both"`. Default is `"both"`.
#' @param use_disk Logical. If TRUE, intermediate raster operations are written to disk using temporary files. Default is FALSE.
#' @param resample_res Numeric or NULL. Target resolution used when coordinates are supplied. If NULL, the coarsest resolution among input rasters is used.
#' @return If `return = "data"`, a data frame with pixel ID, longitude (`x`), latitude (`y`), vegetation cover (`cov`), and elevation (`elv`). If `return = "rasters"`, a list of processed rasters. If `return = "both"`, a list containing processed rasters, the data frame, and CRS information.
#' @examples
#' \dontrun{
#' df <- micronicheR::rast_to_df(
#'   rcover = cover_raster,
#'   rtop = elevation_raster,
#'   rast_or_coord = study_raster,
#'   return = "data"
#' )
#' }
#' @importFrom terra rast vect minmax crs same.crs project resample crop mask extract as.data.frame ext res app
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @export
rast_to_df <- function(
    rcover
    , rtop
    , rast_or_coord
    , method_cover = "near"
    , method_top = "bilinear"
    , return = "both"
    , use_disk = FALSE
    , resample_res = NULL
) {

  # -------------------------------
  # 0. Defining the CRS (always for "WGS"|"WGS84"|"EPSG:4326")
  # -------------------------------

  target_crs <- "EPSG:4326"

  # -------------------------------
  # 1. Detect method
  # -------------------------------

  detect_method <- function(r, user_method = NULL) {

    if (!is.null(user_method)) {
      return(user_method)
    }

    vals <- terra::minmax(r)

    if (all(vals == round(vals), na.rm = TRUE)) {
      return("near")
    } else {
      return("bilinear")
    }
  }

  method_cover <- detect_method(rcover, method_cover)
  method_top   <- detect_method(rtop, method_top)

  # =========================================================
  # COORDINATE MODE
  # =========================================================

  is_coords <- is.matrix(rast_or_coord) || is.data.frame(rast_or_coord)

  if (is_coords) {

    # -------------------------------
    # 2. Prepare coordinates
    # -------------------------------

    coords <- as.data.frame(rast_or_coord)

    if (ncol(coords) < 2) {
      stop("A matrix de coordenadas deve ter pelo menos duas colunas (x, y).")
    }

    names(coords)[1:2] <- c("x", "y")

    # -------------------------------
    # 3. Create points (WGS84)
    # -------------------------------

    pts <- terra::vect(coords, geom = c("x", "y"), crs = "EPSG:4326")

    # -------------------------------
    # 4. Force rasters to WGS84
    # -------------------------------

    force_wgs84 <- function(r) {

      crs_wgs <- "EPSG:4326"

      if (is.na(terra::crs(r))) {
        stop("Raster without a defined CRS. Reprojection is not possible.")
      }

      if (!terra::same.crs(r, crs_wgs)) {
        r <- terra::project(r, crs_wgs)
      }

      return(r)
    }

    rcover <- force_wgs84(rcover)
    rtop   <- force_wgs84(rtop)

    # -------------------------------
    # 4.1 Resolution adjustment
    # -------------------------------

    adjust_resolution <- function(r, res_target, method) {

      if (is.null(res_target)) return(r)

      r_template <- terra::rast(
        extent = terra::ext(r),
        resolution = res_target,
        crs = terra::crs(r)
      )

      terra::resample(r, r_template, method = method)
    }

    # coarser resolution (if NULL)
    if (is.null(resample_res)) {

      res_cover <- terra::res(rcover)
      res_top   <- terra::res(rtop)

      resample_res <- max(c(res_cover, res_top), na.rm = TRUE)
    }

    rcover <- adjust_resolution(rcover, resample_res, method_cover)
    rtop   <- adjust_resolution(rtop, resample_res, method_top)

    # -------------------------------
    # 5. Extraction
    # -------------------------------

    cov_vals <- terra::extract(rcover, pts, method = method_cover)
    top_vals <- terra::extract(rtop, pts, method = method_top)

    df <- data.frame(
      ID  = seq_len(nrow(coords)),
      x   = coords$x,
      y   = coords$y,
      cov = cov_vals[,2],
      elv = top_vals[,2]
    )

    df <- df[stats::complete.cases(df), ]

    # -------------------------------
    # 6. Output
    # -------------------------------

    if (return == "data") {
      return(df)
    }

    if (return == "both") {
      return(list(
        rasters = list(
          cover = rcover,
          top   = rtop
        ),
        data = df,
        crs = "EPSG:4326"
      ))
    }

    if (return == "rasters") {
      return(list(
        cover = rcover,
        top   = rtop,
        crs   = "EPSG:4326"
      ))
    }

    stop("Invalid 'return' argument.")
  }

  # =========================================================
  # ORIGINAL MODE (RASTER)
  # =========================================================

  crs_list <- c(
    terra::crs(rcover),
    terra::crs(rtop),
    terra::crs(rast_or_coord)
  )

  target_crs_name <- paste0(target_crs," | WGS | WGS84")

  if (length(unique(crs_list)) > 1) {
    warning("Different CRS values detected. Reprojecting to: ", target_crs_name)
  }

  reproject_if_needed <- function(r) {
    if (!terra::same.crs(r, target_crs)) {
      r <- terra::project(r, target_crs)
    }
    return(r)
  }

  rcover <- reproject_if_needed(rcover)
  rtop   <- reproject_if_needed(rtop)
  rast_or_coord <- reproject_if_needed(rast_or_coord)

  align_to_study <- function(r, rast_or_coord, method) {

    filename <- if (use_disk) tempfile(fileext = ".tif") else ""

    terra::mask(
      terra::resample(
        terra::crop(r, rast_or_coord, filename = filename),
        rast_or_coord,
        method = method,
        filename = filename
      ),
      rast_or_coord,
      filename = filename
    )
  }

  cover <- align_to_study(rcover, rast_or_coord, method_cover)
  top   <- align_to_study(rtop, rast_or_coord, method_top)
  study <- rast_or_coord

  clean_raster <- function(r, method, tol = 1e-10) {

    if (method == "near") {

      filename <- if (use_disk) tempfile(fileext = ".tif") else ""

      r <- terra::app(r, fun = function(x) {
        x[abs(x) < tol] <- 0
        round(x)
      }, filename = filename)
    }

    return(r)
  }

  cover <- clean_raster(cover, method_cover)
  top   <- clean_raster(top, method_top)

  names(cover) <- "cov"
  names(top)   <- "elv"
  names(study) <- "study"

  if (return == "rasters") {
    return(list(
      cover = cover,
      top   = top,
      crs   = target_crs
    ))
  }

  stack <- c(cover, top)

  df <- terra::as.data.frame(stack, xy = TRUE, na.rm = TRUE) |>
    tibble::rownames_to_column("ID") #%>% tibble::as_tibble()

  if (return == "data") {
    return(df)
  }

  if (return == "both") {
    return(list(
      rasters = list(
        study = study,
        cover = cover,
        top   = top
      ),
      data = df,
      crs = target_crs
    ))
  }

  stop("Invalid 'return' argument.")
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
