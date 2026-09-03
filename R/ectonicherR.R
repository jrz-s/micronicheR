#' Run microclimate-based ectotherm physiological niche simulations
#' @title Microclimate-based physiological niche modeling
#' @description Runs microclimate and ectotherm simulations for one or more species. The function first converts raster or coordinate inputs into a pixel-level environmental data frame using `rast_to_df()`, then runs `NicheMapR::micro_global()` once per pixel and reuses the resulting microclimate for all species in `traits_df`.
#'
#' @author
#' \itemize{
#'   \item Zárate-Salazar, J. Rafael, PhD
#'   \item Valença, Saulo E. S., MSc
#'   \item Castro, Luis Miguel Senzano, PhD
#'   \item Rubalcaba, Juan G., PhD
#'   \item Gouveia, Sidney Feitosa, PhD
#' }
#'
#' @param rcover SpatRaster. Raster containing vegetation cover or shade values.
#' @param rtop SpatRaster. Raster containing elevation or topographic values.
#' @param rast_or_coord SpatRaster, matrix, or data.frame. Study raster or coordinates used to generate the environmental input data frame.
#' @param traits_df data.frame or NULL. Table containing species-level physiological, morphological, and behavioural parameters.
#'
#' @param method_cover Character. Resampling or extraction method for the cover raster passed to `rast_to_df()`. Default is `"near"`.
#' @param method_top Character. Resampling or extraction method for the elevation raster passed to `rast_to_df()`. Default is `"bilinear"`.
#' @param use_disk Logical. If TRUE, intermediate raster operations are written to disk using temporary files. Default is FALSE.
#' @param resample_res Numeric or NULL. Target resolution used by `rast_to_df()` when coordinates are supplied.
#'
#' @param sample_n Integer or NULL. Optional number of pixels or locations to sample before running the microclimate model.
#' @param sample_frac Numeric or NULL. Optional fraction of pixels or locations to sample before running the microclimate model.
#' @param seed Integer or NULL. Optional random seed used when `sample_n` or `sample_frac` is supplied.
#'
#' @param maxshade Numeric. Maximum shade value passed to `NicheMapR::micro_global()`. Default is 100.
#' @param runshade Numeric. Shade mode passed to `NicheMapR::micro_global()`. Default is 1.
#' @param Usrhyt Numeric. User-specified reference height passed to `NicheMapR::micro_global()`. Default is 0.01.
#'
#' @param Ww_g Numeric. Wet body mass (g).
#' @param pct_wet Numeric. Fraction of body surface acting as a free-water exchanger.
#' @param alpha_min Numeric. Minimum solar absorptivity.
#' @param alpha_max Numeric. Maximum solar absorptivity.
#' @param shape Numeric. Body shape code used by `NicheMapR::ectotherm()`.
#' @param T_RB_min Numeric. Minimum body temperature at which the animal attempts to leave its retreat.
#' @param T_B_min Numeric. Minimum body temperature required for basking.
#' @param T_F_min Numeric. Minimum body temperature allowing activity.
#' @param T_F_max Numeric. Maximum body temperature allowing activity.
#' @param T_pref Numeric. Preferred body temperature.
#' @param CT_max Numeric. Critical thermal maximum.
#' @param CT_min Numeric. Critical thermal minimum.
#' @param mindepth Numeric. Minimum burrow depth (soil node).
#' @param maxdepth Numeric. Maximum burrow depth (soil node).
#' @param shade_seek Numeric. Indicates whether shade seeking behaviour is enabled.
#' @param burrow Numeric. Indicates whether burrowing behaviour is enabled.
#' @param climb Numeric. Indicates whether climbing behaviour is enabled.
#' @param nocturn Numeric. Indicates nocturnal activity.
#' @param crepus Numeric. Indicates crepuscular activity.
#' @param diurn Numeric. Indicates diurnal activity.
#'
#' @param break_n Numeric or NULL. Optional chunk size used to split large pixel/location data frames into smaller blocks before running the microclimate and physiological models. This option is useful for reducing memory usage during large simulations.
#'
#' @param list_format Logical. Controls the format of the returned results.
#' If TRUE, returns a nested list organized by species and pixel, with
#' `energy` and `evap` components. Each component contains `DAY`, `time`,
#' `TC` (body temperature), `TA` (air temperature), and the corresponding
#' energy balance (`enbal`) or evaporative water loss (`masbal`).
#' If FALSE, returns a tidy data frame containing, in order, species (`sp`),
#' output type (`type`), pixel ID (`ID`), coordinates (`x`, `y`), day (`DAY`),
#' time (`time`), body temperature (`TC`), air temperature (`TA`), energy
#' balance (`enbal`), and evaporative water loss (`masbal`).
#'
#' @param summary Logical. If TRUE and `list_format = FALSE`, returns daily
#' summary statistics for each species and pixel. The output contains, in
#' order, species (`sp`), pixel ID (`ID`), coordinates (`x`, `y`), day (`DAY`),
#' mean and standard deviation of air temperature (`temp.day.mean`, `temp.day.sd`),
#' mean and standard deviation of body temperature (`body.temp.day.mean`, `body.temp.day.sd`),
#' daily summed energy balance (`enbal.day.sum`),
#' and daily summed evaporative water loss (`masbal.day.sum`).
#'
#' @param checkpoint_dir Character. Directory used to store checkpoint files during block processing. If NULL (default), checkpoint files are not created.
#'
#' @param resume Logical. If TRUE and checkpoint files already exist, previously completed blocks are loaded and skipped. Default is TRUE.
#'
#' @return
#' If `list_format = TRUE`, returns a nested list with one element per species, each containing `energy` (`enbal`) and `evap` (`masbal`) lists for every pixel.
#'
#' If `list_format = FALSE`, returns a tidy data frame containing species, pixel coordinates, day, time, body temperature, air temperature, energy balance (`enbal`), and evaporative water loss (`masbal`).
#'
#' If `summary = TRUE`, returns daily summary statistics for each species and pixel.
#'
#' @examples
#' \dontrun{
#' out <- micronicheR::ectonicheR(
#'     rcover = cover_raster
#'   , rtop = elevation_raster
#'   , rast_or_coord = study_raster
#'   , traits_df = traits_table
#'   , list_format = FALSE
#'   , summary = TRUE
#' )
#' }
#' @importFrom dplyr mutate select rename group_by summarise slice_sample coalesce
#' @importFrom tidyr unnest_longer unnest pivot_wider
#' @importFrom tibble enframe as_tibble
#' @importFrom purrr map
#' @importFrom here here
#' @importFrom stats setNames
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @importFrom rlang .data
#' @export
ectonicheR <- function(

  # ----------------------- #
  # Inputs provided by the researcher
  # ----------------------- #

  # Tree cover raster
    rcover = NULL

  # Elevation raster
  , rtop = NULL

  # Raster or coordinates
  , rast_or_coord = NULL

  # Species traits
  , traits_df = NULL

  # Arguments passed to rast_to_df()
  , method_cover = "near"
  , method_top = "bilinear"
  , use_disk = FALSE
  , resample_res = NULL

  # Optional sampling
  , sample_n = NULL
  , sample_frac = NULL
  , seed = NULL

  # ----------------------- #
  # Microclimate arguments
  # ----------------------- #

  , maxshade = 100
  , runshade = 1
  , Usrhyt = 0.01

  # ----------------------- #
  # Species default parameters
  # ----------------------- #

  , Ww_g       = 40
  , pct_wet    = 0.2
  , alpha_min  = 0.85
  , alpha_max  = 0.85
  , shape      = 3

  , T_RB_min   = 17.5
  , T_B_min    = 17.5
  , T_F_min    = 24
  , T_F_max    = 34
  , T_pref     = 30

  , CT_max     = 40
  , CT_min     = 6

  , mindepth   = 2
  , maxdepth   = 10

  , shade_seek = 1
  , burrow     = 1
  , climb      = 0

  , nocturn    = 0
  , crepus     = 0
  , diurn      = 1

  # ----------------------- #
  # Output and execution
  # ----------------------- #

  , list_format = TRUE
  , summary = FALSE
  , break_n = NULL

  , checkpoint_dir = NULL
  , resume = TRUE

){ # start

  # Install NicheMapR dependencies

  if (!requireNamespace("NicheMapR", quietly = TRUE)) {

    stop(
      "Package 'NicheMapR' is required.\n",
      "Run microniche_setup(github = 'mrke/NicheMapR')."
    )

  }

  # ================================================================
  # 0. CALL raster_to_df()
  # ------------------------------------------------
  # This step calls the external raster_to_df() function to transform
  # cover, elevation, and study rasters or coordinates into the input
  # data frame used internally by ectonicheR().
  # ================================================================

  df <- rast_to_df(
    rcover = rcover
    , rtop = rtop
    , rast_or_coord = rast_or_coord
    , method_cover = method_cover
    , method_top = method_top
    , return = "data"
    , use_disk = use_disk
    , resample_res = resample_res
  )

  # Optional sampling of pixels/locations

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(sample_n) && !is.null(sample_frac)) {
    stop("Use only one of 'sample_n' or 'sample_frac', not both.")
  }

  if (!is.null(sample_n)) {
    df <- dplyr::slice_sample(df, n = sample_n)
  }

  if (!is.null(sample_frac)) {
    df <- dplyr::slice_sample(df, prop = sample_frac)
  }

  # Optional chunking for large runs (28/05/2026)

  if (!is.null(break_n)) {

    if (!is.numeric(break_n) || length(break_n) != 1 || break_n < 1) {
      stop("'break_n' must be NULL or a single positive integer.")
    }

    break_n <- as.integer(break_n)

  }

  # ================================================================
  # 1. INTERNAL FUNCTION: micro_by_pixel()
  # ------------------------------------------------
  # This internal function runs NicheMapR::micro_global()
  # ONLY ONCE for each pixel/location.
  #
  # Input:
  #   df containing ID, x, y, cov, and elv
  #
  # Output:
  #   micro_pixels, a list containing the stored microclimate
  #   for each pixel.
  #
  # Important:
  #   The user does NOT need to run micro_by_pixel()
  #   outside ectonicheR().
  # ================================================================

  micro_by_pixel <- function(df){

    micro_pixels <- vector("list", nrow(df))

    for(l in seq_len(nrow(df))){

      message("Running microclimate for pixel ID: ", df$ID[l])

      # ------------------------------------------------------------
      # Ensure minshade < maxshade
      # ------------------------------------------------------------

      if(df$cov[l] >= maxshade){

        warning(
            "Pixel "
          , df$ID[l]
          , ": minshade (", df$cov[l]
          , ") >= maxshade (", maxshade
          , "). Adjusting minshade automatically."
        )

        minshade_i <- maxshade - 1e-6

      } else {

        minshade_i <- df$cov[l]

      }

      # ------------------------------------------------------------
      # Run micro_global()
      # ------------------------------------------------------------

      micro <- tryCatch(

        NicheMapR::micro_global(

            loc      = c(df$x[l], df$y[l])
          , dem      = df$elv[l]

          , minshade = minshade_i
          , maxshade = maxshade

          , runshade = runshade
          , Usrhyt   = Usrhyt

        ),

        error = function(e){

          warning(
            "micro_global() failed for pixel ",
            df$ID[l],
            ": ",
            conditionMessage(e)
          )

          NULL

        }

      )

      micro_pixels[[l]] <- list(

        ID    = df$ID[l],
        x     = df$x[l],
        y     = df$y[l],
        cov   = df$cov[l],
        elv   = df$elv[l],
        micro = micro

      )

    }

    micro_pixels

  }

  # ================================================================
  # 2. INTERNAL FUNCTION: ecto_by_species()
  # ------------------------------------------------
  # Runs NicheMapR::ectotherm() once for each species
  # and each pixel.
  # ================================================================

  ecto_by_species <- function(
    micro_pixels
    , traits_df
  ){

    niche_list <- list()

    for(s in seq_len(nrow(traits_df))){

      pars <- traits_df[s, ]

      sp_name <- as.character(pars$sp)

      if(is.na(sp_name) || sp_name == ""){
        stop(
          "Column 'sp' in traits_df contains missing or empty species names."
        )
      }

      message(
        "Running ectotherm for species: "
        , sp_name
      )

      mod.xy.eb <- vector(
        "list"
        , length(micro_pixels)
      )

      mod.xy.mb <- vector(
        "list"
        , length(micro_pixels)
      )

      for(l in seq_along(micro_pixels)){

        pix <- micro_pixels[[l]]

        # Skip pixels whose microclimate simulation failed
        if(is.null(pix$micro)){

          warning(
              "Skipping pixel "
            , pix$ID
            , " because microclimate simulation failed."
          )

          next

        }

        micro <- pix$micro

        if(
          is.null(micro) ||
          !is.list(micro) ||
          is.null(micro$metout)
        ){
          next
        }

        # ----------------------------------------------------------
        # Build argument list for NicheMapR::ectotherm()
        # ----------------------------------------------------------

        args_ecto <- list(

          # --------------------------------------------------------
          # Morphological parameters
          # --------------------------------------------------------

          Ww_g      = pars$Ww_g,
          shape     = pars$shape,
          alpha_max = pars$alpha_max,
          alpha_min = pars$alpha_min,
          pct_wet   = pars$pct_wet,

          # --------------------------------------------------------
          # Thermal physiology
          # --------------------------------------------------------

          T_F_max   = pars$T_F_max,
          T_F_min   = pars$T_F_min,
          T_B_min   = pars$T_B_min,
          T_RB_min  = pars$T_RB_min,
          T_pref    = pars$T_pref,
          CT_max    = pars$CT_max,
          CT_min    = pars$CT_min,

          # --------------------------------------------------------
          # Behavioural parameters
          # --------------------------------------------------------

          shade_seek = pars$shade_seek,
          burrow     = pars$burrow,
          climb      = pars$climb,
          nocturn    = pars$nocturn,
          diurn      = pars$diurn,
          crepus     = pars$crepus,
          mindepth   = pars$mindepth,
          maxdepth   = pars$maxdepth,
          postur     = 1,

          # --------------------------------------------------------
          # Output control
          # --------------------------------------------------------

          write_input = 0,
          write_csv   = 0,

          # --------------------------------------------------------
          # Microclimate data
          # --------------------------------------------------------

          nyears     = micro$nyears,
          metout     = micro$metout,
          shadmet    = micro$shadmet,
          soil       = micro$soil,
          shadsoil   = micro$shadsoil,
          soilmoist  = micro$soilmoist,
          shadmoist  = micro$shadmoist,
          humid      = micro$humid,
          shadhumid  = micro$shadhumid,
          soilpot    = micro$soilpot,
          shadpot    = micro$shadpot,
          tcond      = micro$tcond,
          shadtcond  = micro$shadtcond,
          rainfall   = micro$RAINFALL,

          # --------------------------------------------------------
          # Atmospheric conditions
          # --------------------------------------------------------

          rainhr = rep(
            -1,
            nrow(micro$metout)
          ),

          preshr = rep(

            101325 *
              (
                (1 - (0.0065 * as.numeric(micro$elev) / 288))
                ^
                  (1 / 0.190284)
              ),

            nrow(micro$metout)

          ),

          # --------------------------------------------------------
          # Soil properties
          # --------------------------------------------------------

          DEP = micro$DEP,

          KS  = micro$KS[
            seq(1, length(micro$KS), 2)
          ],

          b   = micro$BB[
            seq(1, length(micro$BB), 2)
          ],

          PE  = -abs(

            micro$PE[
              seq(1, length(micro$PE), 2)
            ]

          ),

          # --------------------------------------------------------
          # Site characteristics
          # --------------------------------------------------------

          alpha_sub = 1 - micro$REFL,
          elev      = micro$elev,
          longitude = micro$longlat[1],
          latitude  = micro$longlat[2],

          # --------------------------------------------------------
          # Shade regime
          # --------------------------------------------------------

          minshades = rep(
            micro$minshade,
            each = 24
          ),

          maxshades = rep(
            micro$maxshade,
            each = 24
          )

        )

        mod.t <- tryCatch(

          do.call(

            NicheMapR::ectotherm,
            args_ecto

          ),

          error = function(e){

            warning(

              "ectotherm() failed for species '",
              sp_name,
              "', pixel ",
              pix$ID,
              ": ",
              conditionMessage(e)

            )

            NULL

          }

        )

        if(
          !is.list(mod.t) ||
          is.null(mod.t$metout) ||
          is.null(mod.t$enbal) ||
          is.null(mod.t$masbal)
        ){
          next
        }

        # ----------------------------------------------------------
        # Extract hourly outputs
        # ----------------------------------------------------------

        met <- mod.t$metout
        env <- mod.t$environ

        DAY  <- as.numeric(met[, "DOY"])
        TIME <- as.numeric(met[, "TIME"])
        TC   <- as.numeric(env[, "TC"])
        TA   <- as.numeric(met[, "TALOC"])

        qgens <- as.numeric(
          mod.t$enbal[, "ENB"]
        )

        if(
          !all(
            c("H2OResp_g", "H2OCut_g") %in%
            colnames(mod.t$masbal)
          )
        ){
          next
        }

        evap <- as.numeric(

          rowSums(

            mod.t$masbal[
              ,
              c("H2OResp_g", "H2OCut_g")
            ],

            na.rm = TRUE

          )

        )

        # ----------------------------------------------------------
        # Store outputs
        # ----------------------------------------------------------

        mod.xy.eb[[l]] <- data.frame(

            DAY         = DAY
          , time        = TIME
          , TC          = TC
          , TA          = TA
          , enbal       = qgens

          , check.names = FALSE

        )

        mod.xy.mb[[l]] <- data.frame(

            DAY         = DAY
          , time        = TIME
          , TC          = TC
          , TA          = TA
          , masbal      = evap

          , check.names = FALSE

        )

      }

      niche_list[[sp_name]] <- list(

          energy = mod.xy.eb
        , evap = mod.xy.mb

      )

    }

    return(niche_list)

  }

  # ================================================================
  # 3. INTERNAL FUNCTION: merge_niche_lists()
  # ------------------------------------------------
  # This internal function combines the outputs generated when
  # ectonicheR() is run in blocks using the break_n argument.
  #
  # Input:
  #   niche_chunks = a list where each element corresponds to one
  #                  processed block of pixels/locations.
  #
  # Each block contains the same hierarchical structure returned by
  # ecto_by_species():
  #
  #   species
  #   energy
  #   evap
  #
  # Output:
  #   merged, a single hierarchical list containing all species and
  #   all pixels from all processed blocks.
  #
  # Important:
  #   This function does NOT run any model.
  #   It only joins the results produced by previous block-level runs.
  #
  #   The order of pixels is preserved because the blocks are processed
  #   sequentially from the original df.
  # ================================================================

  merge_niche_lists <- function(niche_chunks) {

    sp_names <- unique(unlist(lapply(niche_chunks, names)))

    merged <- list()

    for (sp in sp_names) {

      merged[[sp]] <- list(
        energy = unlist(
          lapply(niche_chunks, function(x) x[[sp]]$energy),
          recursive = FALSE
        ),
        evap = unlist(
          lapply(niche_chunks, function(x) x[[sp]]$evap),
          recursive = FALSE
        )
      )
    }

    return(merged)
  }

  # ================================================================
  # 4. MAIN FLOW OF ectonicheR()
  # ------------------------------------------------
  # From this point onward, the main workflow of the
  # function begins.
  #
  # The internal functions above have already been
  # defined, but they have not been executed yet.
  # ================================================================

  # 4.1. Check global climate data

  gcfolder_file <- file.path(.libPaths()[1], "gcfolder.rda")

  if (!file.exists(gcfolder_file)) {

    stop(
        "NicheMapR global climate data were not found.\n\n"
      , "Please run this once before using ectonicheR():\n\n"
      , "micronicheR::microniche_setup_global_climate(folder = getwd())\n\n"
      , "After the download finishes, restart R and run ectonicheR() again."
    )

  } else {

    message(
      "NicheMapR global climate data are registered. Continuing process..."
    )

  }

  # 4.2. If traits_df is not provided, create a default species
  # Default single-species behavior.

  if (is.null(traits_df)) {

    message("traits_df not provided. Using default parameters for one species.")

    traits_df <- tibble::tibble(

      sp = "sp1"

      # ------------------------------------------------
      , Ww_g = 40        # wet weight of animal (g)
      , pct_wet = 0.2    # % of surface area acting as a free-water exchanger
      , alpha_min =0.85  # minimum solar absorbtivity (dec %)
      , alpha_max = 0.85 # maximum solar absorbtivity (dec %)
      , shape = 3        # lizard shape
      , T_RB_min = 17.5  # min Tb at which they will attempt to leave retreat
      , T_B_min = 17.5   # min Tb at which leaves retreat to bask
      , T_F_min = 24     # minimum Tb at which activity occurs
      , T_F_max = 34     # maximum Tb at which activity occurs
      , T_pref = 30      # preferred Tb (will try and regulate to this)
      , CT_max = 40      # critical thermal minimum (affects choice of retreat)
      , CT_min = 6       # critical thermal maximum (affects choice of retreat)
      , mindepth = 2     # min depth (node, 1-10) allowed
      , maxdepth = 10    # max depth (node, 1-10) allowed
      , shade_seek = 1   # shade seeking?
      , burrow = 1       # can it burrow?
      , climb = 0        # can it climb to thermoregulate?
      , nocturn = 0      # nocturnal activity
      , crepus = 0       # crepuscular activity
      , diurn = 1        # diurnal activity

    )

  } else {

    message("traits_df provided. Running ectotherm simulations for multiple species.")
  }

  # 4.3. Run microclimate only once per pixel
  # The internal function micro_by_pixel() is called here.

  if (!requireNamespace("NicheMapR", quietly = TRUE)) {
    stop(
        "Package 'NicheMapR' is required.\n"
      , "Please install it with:\n"
      , "micronicheR::microniche_setup(github = 'mrke/NicheMapR')"
    )
  }

  # Note:
  # NicheMapR must be attached because both
  # micro_global() and ectotherm()
  # rely on internal package objects
  # that are not available through
  # requireNamespace() alone.

  if (!"package:NicheMapR" %in% search()) {

    suppressPackageStartupMessages(
      require("NicheMapR", character.only = TRUE)
    )
  }

  # 4.4. Run the physiological model for each species
  # The internal function ecto_by_species() is called here.

  if (is.null(break_n)) {

    micro_pixels <- micro_by_pixel(df = df)

    niche_list <- ecto_by_species(
        micro_pixels = micro_pixels
      , traits_df = traits_df
    )

  } else {

    df_chunks <- split(
        df
      , ceiling(seq_len(nrow(df)) / break_n)
    )

    message(
        "Running ectonicheR in "
      , length(df_chunks)
      , " block(s) of up to "
      , break_n
      , " pixel(s)."
    )

    niche_chunks <- vector("list", length(df_chunks))

    # =====================================================================
    # >>> CHECKPOINT SYSTEM (Start)
    # =====================================================================

    micro_dir <- NULL
    ecto_dir  <- NULL
    micro_file <- NULL
    ecto_file  <- NULL

    if(!is.null(checkpoint_dir)){

      checkpoint_dir <- file.path(
          checkpoint_dir
        , "ectonicheR"
      )

      micro_dir <- file.path(
          checkpoint_dir
        , "micro_global"
      )

      ecto_dir <- file.path(
          checkpoint_dir
        , "ectotherm"
      )

      dir.create(
          micro_dir
        , recursive    = TRUE
        , showWarnings = FALSE
      )

      dir.create(
          ecto_dir
        , recursive    = TRUE
        , showWarnings = FALSE
      )

      message(

          "Checkpoint files will be stored in:\n"
        , normalizePath(
            checkpoint_dir
          , winslash = "/"
          , mustWork = FALSE
        ),
          "\n\n"
        , "The default checkpoint directory can be changed using the "
        , "'checkpoint_dir' argument."

      )

    }

    for (b in seq_along(df_chunks)) {

      message(
          "Running block "
        , b
        , " of "
        , length(df_chunks)
        , " | pixels: "
        , nrow(df_chunks[[b]])
      )

      if(is.null(checkpoint_dir)){

        micro_file <- NULL
        ecto_file  <- NULL

      } else {

        micro_file <- file.path(

          micro_dir

          , paste0(

              "micro_block_"
            , sprintf("%04d", b)
            , ".rds"

          )

        )

        ecto_file <- file.path(

          ecto_dir

          , paste0(

             "ecto_block_"
            , sprintf("%04d", b)
            , ".rds"

          )

        )

      }

      # ------------------------------------------------------------
      # Microclimate simulation
      # ------------------------------------------------------------

      if (
        resume &&
        !is.null(micro_file) &&
        file.exists(micro_file)
      ) {

        message(
            "Loading saved microclimate block "
          , b
        )

        micro_pixels <- readRDS(micro_file)

      } else {

        micro_pixels <- micro_by_pixel(
          df = df_chunks[[b]]
        )

        if (!is.null(micro_file)) {

          message("Saving microclimate block ", b)

          saveRDS(
              micro_pixels
            , micro_file
          )

        }

      }

      # ------------------------------------------------------------
      # Physiological simulation (ectotherm)
      # ------------------------------------------------------------

      if (
        resume &&
        !is.null(ecto_file) &&
        file.exists(ecto_file)
      ) {

        message(
            "Loading saved ectotherm block "
          , b
        )

        niche_chunks[[b]] <- readRDS(
          ecto_file
        )

      } else {

        niche_chunks[[b]] <- ecto_by_species(
            micro_pixels = micro_pixels
          , traits_df = traits_df
        )

        if (!is.null(ecto_file)) {

          message("Saving ectotherm block ", b)

          saveRDS(
              niche_chunks[[b]]
            , ecto_file
          )

        }

      }

    }

    message(
      "All blocks processed successfully. Merging outputs..."
    )

    # =====================================================================
    # >>> CHECKPOINT SYSTEM (End)
    # =====================================================================

    niche_list <- merge_niche_lists(niche_chunks)

    message(
      "ectonicheR finished successfully."
    )

  }

  # 4.5. Control the three output formats

  if(list_format == TRUE){

    return(niche_list)

  } else{

    niche_df <- niche_list %>%
      tibble::enframe(name = "sp", value = "data_sp") %>%
      tidyr::unnest_longer(.data$data_sp, indices_to = "type") %>%
      tidyr::unnest_longer(.data$data_sp, indices_to = "ID") %>%
      dplyr::mutate(
          x  = df$x[.data$ID]
        , y  = df$y[.data$ID]
        , ID = df$ID[.data$ID]
      ) %>%
      dplyr::mutate(data_sp = purrr::map(.data$data_sp, tibble::as_tibble)) %>%
      tidyr::unnest(.data$data_sp) %>%
      dplyr::mutate(
        flux = dplyr::coalesce(.data$enbal, .data$masbal)
      ) %>%
      dplyr::select(
          .data$sp, .data$type, .data$ID, .data$x, .data$y
        , .data$DAY, .data$time, .data$TC, .data$TA, .data$flux
      ) %>%
      tidyr::pivot_wider(
          names_from  = .data$type
        , values_from = .data$flux
      ) %>%
      dplyr::rename(
          enbal  = .data$energy
        , masbal = .data$evap
      )

    if(summary == FALSE){

      return(niche_df)

    } else{

      niche_df_summary <- niche_df %>%
        dplyr::group_by(.data$sp, .data$ID, .data$x, .data$y, .data$DAY) %>%
        dplyr::summarise(

            temp.day.mean       = mean(.data$TA, na.rm = TRUE)
          , temp.day.sd         = stats::sd(.data$TA, na.rm = TRUE)

          , body.temp.day.mean  = mean(.data$TC, na.rm = TRUE)
          , body.temp.day.sd    = stats::sd(.data$TC, na.rm = TRUE)

          , enbal.day.sum       = sum(.data$enbal, na.rm = TRUE)
          , masbal.day.sum      = sum(.data$masbal, na.rm = TRUE)
          , .groups             = "drop"
        )

      return(niche_df_summary)
    }
  }

} # end

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
