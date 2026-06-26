#' Run microclimate-based endotherm physiological niche simulations
#' @title Microclimate-based physiological niche modeling
#' @description Runs microclimate and endotherm physiological niche simulations for one or more species. The function first converts raster or coordinate inputs into a pixel-level environmental data frame using `rast_to_df()`, then runs `NicheMapR::micro_global()` once per pixel and reuses the resulting microclimate for all species in `traits_df`.
#' @author
#' \itemize{
#'   \item Zárate-Salazar, J. Rafael, PhD
#'   \item Valença, Saulo E. S., MSc
#'   \item Castro, Luis Miguel Senzano, PhD
#'   \item Rubalcaba, Juan G., PhD
#'   \item Gouveia, Sidney Feitosa, PhD
#' }
#' @param rcover SpatRaster. Raster containing vegetation cover or shade values.
#' @param rtop SpatRaster. Raster containing elevation or topographic values.
#' @param rast_or_coord SpatRaster, matrix, or data.frame. Study raster or coordinates used to generate the environmental input data frame.
#' @param traits_df data.frame or NULL. Table containing species-level physiological and morphological parameters. Must include a species identifier column named `sp`. If NULL, default parameters are used for a single species named `"sp1"`.
#' @param method_cover Character. Resampling or extraction method for the cover raster passed to `rast_to_df()`. Default is `"near"`.
#' @param method_top Character. Resampling or extraction method for the elevation raster passed to `rast_to_df()`. Default is `"bilinear"`.
#' @param use_disk Logical. If TRUE, intermediate raster operations are written to disk using temporary files. Default is FALSE.
#' @param resample_res Numeric or NULL. Target resolution used by `rast_to_df()` when coordinates are supplied.
#' @param sample_n Integer or NULL. Optional number of pixels or locations to sample before running the microclimate model.
#' @param sample_frac Numeric or NULL. Optional fraction of pixels or locations to sample before running the microclimate model.
#' @param seed Integer or NULL. Optional random seed used when `sample_n` or `sample_frac` is supplied.
#' @param maxshade Numeric. Maximum shade value passed to `NicheMapR::micro_global()`. Default is 100.
#' @param runshade Numeric. Shade mode passed to `NicheMapR::micro_global()`. Default is 1.
#' @param Usrhyt Numeric. User-specified reference height passed to `NicheMapR::micro_global()`. Default is 0.01.
#' @param Z Numeric. Solar zenith angle in degrees.
#' @param ABSSB Numeric. Substrate solar absorptivity.
#' @param FLTYPE Numeric. Fluid type: 0 = air, 1 = freshwater, 2 = saltwater.
#' @param KSUB Numeric. Substrate thermal conductivity.
#' @param BP Numeric. Barometric pressure. Negative values indicate that elevation is used.
#' @param O2GAS Numeric. Oxygen concentration in air.
#' @param N2GAS Numeric. Nitrogen concentration in air.
#' @param CO2GAS Numeric. Carbon dioxide concentration in air.
#' @param PDIF Numeric. Proportion of diffuse solar radiation.
#' @param FLYHR Numeric. Indicates whether flight is occurring at a given hour.
#' @param UNCURL Numeric. Increment controlling uncurling toward maximum body shape ratio.
#' @param TC_INC Numeric. Increment for increasing core temperature.
#' @param PCTWET_INC Numeric. Increment for increasing wet skin surface area.
#' @param PCTWET_MAX Numeric. Maximum wet skin surface area.
#' @param AK1_INC Numeric. Increment for increasing flesh thermal conductivity.
#' @param AK1_MAX Numeric. Maximum flesh thermal conductivity.
#' @param PANT Numeric. Respiratory frequency multiplier for panting.
#' @param PANT_INC Numeric. Increment for panting multiplier.
#' @param PANT_MULT Numeric. Basal metabolic rate multiplier at maximum panting.
#' @param AMASS Numeric. Body mass in kg.
#' @param ANDENS Numeric. Animal density in kg/m3.
#' @param SUBQFAT Numeric. Indicates whether subcutaneous fat is present.
#' @param FATPCT Numeric. Body fat percentage.
#' @param SHAPE Numeric. Body shape code.
#' @param SHAPE_B Numeric. Ratio of long to short body axis.
#' @param SHAPE_B_MAX Numeric. Maximum ratio of long to short body axis.
#' @param PVEN Numeric. Fraction of surface area covered by ventral fur.
#' @param PCOND Numeric. Fraction of surface area in contact with the substrate.
#' @param SAMODE Numeric. Surface area calculation mode.
#' @param ORIENT Numeric. Animal orientation relative to solar radiation.
#' @param FURTHRMK Numeric. User-specified fur thermal conductivity.
#' @param DHAIRD Numeric. Dorsal hair diameter.
#' @param DHAIRV Numeric. Ventral hair diameter.
#' @param LHAIRD Numeric. Dorsal hair length.
#' @param LHAIRV Numeric. Ventral hair length.
#' @param ZFURD Numeric. Dorsal fur depth.
#' @param ZFURV Numeric. Ventral fur depth.
#' @param RHOD Numeric. Dorsal hair density.
#' @param RHOV Numeric. Ventral hair density.
#' @param REFLD Numeric. Dorsal fur reflectivity.
#' @param REFLV Numeric. Ventral fur reflectivity.
#' @param KHAIR Numeric. Hair thermal conductivity.
#' @param XR Numeric. Fractional fur depth at which longwave radiation is exchanged.
#' @param EMISAN Numeric. Animal emissivity.
#' @param FABUSH Numeric. Vegetation factor below or around the animal.
#' @param FGDREF Numeric. Ground reference configuration factor.
#' @param FSKREF Numeric. Sky configuration factor.
#' @param TC Numeric. Core temperature.
#' @param TC_MAX Numeric. Maximum core temperature.
#' @param AK1 Numeric. Initial flesh thermal conductivity.
#' @param AK2 Numeric. Fat conductivity.
#' @param PCTWET Numeric. Percentage of wet skin surface.
#' @param FURWET Numeric. Percentage of fur or feathers wet after rain.
#' @param PCTBAREVAP Numeric. Bare skin surface area available for evaporation.
#' @param PCTEYES Numeric. Surface area composed of eyes.
#' @param DELTAR Numeric. Offset between air and respiratory temperature.
#' @param RELXIT Numeric. Relative humidity of exhaled air.
#' @param RQ Numeric. Respiratory quotient.
#' @param EXTREF Numeric. Oxygen extraction efficiency.
#' @param PANT_MAX Numeric. Maximum panting multiplier.
#' @param PZFUR Numeric. Fractional reduction in fur depth from the piloerected state.
#' @param Q10 Numeric. Q10 factor for metabolic adjustment.
#' @param TC_MIN Numeric. Minimum core temperature during torpor.
#' @param DIFTOL Numeric. Numerical tolerance for the solver.
#' @param THERMOREG Numeric. Indicates whether thermoregulatory responses are invoked.
#' @param RESPIRE Numeric. Indicates whether respiration and associated heat loss are calculated.
#' @param TREGMODE Numeric. Thermoregulation mode.
#' @param WRITE_INPUT Numeric. Argument passed to `NicheMapR::endoR()` controlling whether input is written.
#' @param break_n Numeric or NULL. Optional chunk size used to split large pixel/location data frames into smaller blocks before running the microclimate and physiological models. If NULL, all pixels are processed at once. For example, `break_n = 500` processes the input in blocks of up to 500 pixels.
#' @param list_format Logical. If TRUE, returns a nested list by species, energy balance, evaporation, and pixel. If FALSE, returns a tidy data frame.
#' @param summary Logical. If TRUE and `list_format = FALSE`, returns a daily summary by species and pixel.
#' @param checkpoint_dir
#' Character. Directory used to store checkpoint files generated
#' during long simulations. If NULL (default), checkpoint files
#' are not created.
#' @param resume
#' Logical. If TRUE and checkpoint files already exist,
#' previously completed blocks are loaded and skipped.
#' Default is TRUE.
#' @return If `list_format = TRUE`, a nested list with one element per species, each containing `energy` and `evap` lists by pixel. If `list_format = FALSE`, a data frame with species, pixel coordinates, day, time, air temperature, energy balance, and mass balance. If `summary = TRUE`, returns daily summary statistics.
#' @examples
#' \dontrun{
#' out <- micronicheR::endonicheR(
#'   rcover = cover_raster,
#'   rtop = elevation_raster,
#'   rast_or_coord = study_raster,
#'   traits_df = traits_table,
#'   list_format = FALSE,
#'   summary = TRUE
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
endonicheR <- function(

  # ----------------------- #
  # Inputs provided by the researcher
  # ----------------------- #

  # Coverage raster
  rcover = NULL
  # topographic raster
  , rtop = NULL
  # Raster or coordinates of the study
  , rast_or_coord = NULL
  # Dataframe with the traits provided by the researcher
  , traits_df = NULL

  # arguments passed to rast_to_df()
  , method_cover = "near"
  , method_top = "bilinear"
  , use_disk = FALSE
  , resample_res = NULL

  # conditional samples
  , sample_n = NULL
  , sample_frac = NULL
  , seed = NULL

  # ----------------------- #
  # Environmental arguments
  # ----------------------- #
  , maxshade = 100
  , runshade = 1
  , Usrhyt = 0.01

  # ----------------------- #
  # Arguments of the species
  # ----------------------- #

  , Z = 90 # solar zenith angle (degrees)
  , ABSSB = 0.8 # substrate solar absorptivity (fractional, 0-1)
  , FLTYPE = 0 # fluid type: 0 = air; 1 = freshwater; 2 = saltwater
  , KSUB = 2.79 # substrate thermal conductivity
  , BP = -1 # Pa; negative value means elevation is used
  , O2GAS = 20.95 # oxygen concentration in air, used to account for non-atmospheric concentrations, e.g., burrows
  , N2GAS = 79.02 # nitrogen concentration in air, used to account for non-atmospheric concentrations, e.g., burrows
  , CO2GAS = 0.0412 # carbon dioxide concentration in air, used to account for non-atmospheric concentrations, e.g., burrows
  , PDIF = 0.15 # proportion of solar radiation that is diffuse (fractional, 0-1)

  # ---- BEHAVIOR ----

  , FLYHR = 0 # is flight occurring at this hour? (imposes forced evaporative loss)
  , UNCURL = 0.1 # allows the animal to uncurl toward SHAPE_B_MAX; the value is the increment by which SHAPE_B is increased per iteration
  , TC_INC = 0.1 # enables core temperature elevation; the value is the increment by which TC is increased per iteration
  , PCTWET_INC = 0.1 # enables sweating; the value is the increment by which PCTWET is increased per iteration
  , PCTWET_MAX = 70 # maximum surface area that can be wet
  , AK1_INC = 0.1 # enables increased thermal conductivity (W/mK); the value is the increment by which AK1 is increased per iteration
  , AK1_MAX = 2.8 # maximum flesh conductivity
  , PANT = 1 # multiplier on respiratory frequency to simulate panting
  , PANT_INC = 0.1 # increment for the respiratory frequency multiplier to simulate panting
  , PANT_MULT = 1.05 # multiplier on basal metabolic rate at the maximum panting level

  # ---- MORPHOLOGY ----

  # Geometry
  , AMASS = 1.2 # kg (body mass) (Jerusalinsky, 2013)
  , ANDENS = 1000 # kg/m3 (density)
  , SUBQFAT = 0 # is subcutaneous fat present? (0 = no, 1 = yes)
  , FATPCT = 20 # body fat percentage
  , SHAPE = 4 # shape (1 = cylinder, 2 = sphere, 3 = plate, 4 = ellipsoid)
  , SHAPE_B = 7 # ratio of long to short axis, must be > 1 (-)
  , SHAPE_B_MAX = 7 # maximum possible ratio of long to short axis, must be > 1 (-)
  , PVEN = 0.3 # fraction of surface area covered by ventral fur (fractional, 0-1)
  , PCOND = 0 # fraction of surface area in contact with the substrate (fractional, 0-1)
  , SAMODE = 2 # if 0, use surface area based on SHAPE geometry; if 1, use bird skin surface area allometry from Walsberg & King 1978, JEB 76:185-189; if 2, use mammal surface area from Stahl 1967, J. Appl. Physiol. 22:453-460
  , ORIENT = 1 # if 1 = normal to solar rays (heat maximization), if 2 = parallel to solar rays (heat minimization), if 3 = vertical and varying with solar altitude, or if 0 = average

  # Fur properties
  , FURTHRMK = 0 # user-specified fur thermal conductivity (W/mK), not used if 0
  , DHAIRD = 30E-06 # hair diameter, dorsal (m)
  , DHAIRV = 30E-06 # hair diameter, ventral (m)
  , LHAIRD = 57E-03 # hair length, dorsal (m)
  , LHAIRV = 44E-03 # hair length, ventral (m)
  , ZFURD = 22E-03 # fur depth, dorsal (m)
  , ZFURV = 17E-03 # fur depth, ventral (m)
  , RHOD = 3000E+04 # hair density, dorsal (1/m2)
  , RHOV = 3000E+04 # hair density, ventral (1/m2)
  , REFLD = 0.2 # dorsal fur reflectivity (fractional, 0-1)
  , REFLV = 0.2 # ventral fur reflectivity (fractional, 0-1)
  , KHAIR = 0.209 # hair thermal conductivity
  , XR = 1 # fractional fur depth at which longwave radiation is exchanged (0-1)

  # Radiation exchange
  , EMISAN = 0.99 # animal emissivity (-)
  , FABUSH = 0 # for vegetation below/around the animal (in TALOC)
  , FGDREF = 0.5 # ground reference configuration factor
  , FSKREF = 0.5 # sky configuration factor

  # Thermal physiology
  , TC = 37 # core temperature
  , TC_MAX = 39 # maximum core temperature
  , AK1 = 0.9 # initial flesh thermal conductivity (0.412 - 2.8 W/m per degree Celsius)
  , AK2 = 0.230 # fat conductivity

  # Evaporation
  , PCTWET = 0.5 # portion of the skin surface that is wet (%)
  , FURWET = 0 # portion of fur/feathers that becomes wet after rain (%)
  , PCTBAREVAP = 0 # surface area available for evaporation that is bare skin, e.g., paw licking (%)
  , PCTEYES = 0 # surface area composed of eyes (%) - zero if sleeping
  , DELTAR = 0 # offset between air temperature and respiratory temperature (per degree Celsius)
  , RELXIT = 100 # relative humidity of exhaled air (%)

  # Metabolism | Respiration
  , RQ = 0.80 # respiratory quotient (fractional, 0-1)
  , EXTREF = 20 # O2 extraction efficiency (%)
  , PANT_MAX = 5 # maximum respiratory frequency multiplier to simulate panting (-)
  , PZFUR = 0 # fractional incremental reduction in ZFUR from the piloerected state (-); values greater than zero trigger a piloerection response
  , Q10 = 2 # Q10 factor for adjusting BMR to TC
  , TC_MIN = 19 # minimum core temperature during torpor (TREGMODE = 0)

  # Other model settings
  , DIFTOL = 0.001 # tolerance for SIMULSOL
  , THERMOREG = 1 # invoke thermoregulatory response
  , RESPIRE = 1 # calculate respiration and associated heat loss
  , TREGMODE = 1 # 0 = torpor, 1 = raise core temperature, then pant, then sweat; 2 = raise core temperature and pant simultaneously, then sweat

  , WRITE_INPUT = 0

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
  # data frame used internally by endonicheR().
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

    # =====================================================================
    # >>> V2 CHECKPOINT SYSTEM (BEGIN)
    # =====================================================================

    if (!is.null(checkpoint_dir)) {

      message(
        "Checkpoint directory: ",
        normalizePath(
          checkpoint_dir,
          winslash = "/",
          mustWork = FALSE
        )
      )

      if (!dir.exists(checkpoint_dir)) {

        dir.create(
          checkpoint_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )

      }

    }

    # =====================================================================
    # >>> V2 CHECKPOINT SYSTEM (END)
    # =====================================================================

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
  #   for each pixel
  #
  # Important:
  #   The user does NOT need to run micro_by_pixel()
  #   outside endonicheR().
  #   It exists only to organize the internal workflow
  #   of the function.
  # ================================================================

  micro_by_pixel <- function(df){

    micro_pixels <- vector("list", nrow(df))

    for(l in 1:nrow(df)){

      message("Running microclimate for pixel ID: ", df$ID[l])

      micro <- NicheMapR::micro_global(
        loc = c(df$x[l], df$y[l])
        , dem = df$elv[l]
        , minshade = df$cov[l]
        , maxshade = maxshade
        , runshade = runshade
        , Usrhyt = Usrhyt
      )

      micro_pixels[[l]] <- list(
        ID = df$ID[l]
        , x = df$x[l]
        , y = df$y[l]
        , cov = df$cov[l]
        , elv = df$elv[l]
        , tshade = as.data.frame(micro$shadmet)
        , shadsoil = as.data.frame(micro$shadsoil)
        , metout = as.data.frame(micro$metout)
      )
    }

    return(micro_pixels)
  }

  # ================================================================
  # 2. INTERNAL FUNCTION: endo_by_species()
  # ------------------------------------------------
  # This internal function runs NicheMapR::endoR()
  # for each species and each pixel.
  #
  # Input:
  #   micro_pixels = microclimate previously calculated
  #                  by micro_by_pixel()
  #
  #   traits_df    = table containing species parameters
  #
  # Output:
  #   niche_df, a data frame containing sp, ID, x, y,
  #   DAY, time, TA, enbal, and masbal.
  #
  # Important:
  #   endo_by_species() does NOT create traits_df.
  #   traits_df is either provided by the user or
  #   automatically created when traits_df = NULL.
  # ================================================================

  endo_by_species <- function(micro_pixels, traits_df){

    niche_list <- list()

    for(s in 1:nrow(traits_df)){

      pars <- traits_df[s, ]
      sp_name <- pars$sp

      message("Running endoR for species: ", sp_name)

      mod.xy.eb <- list()
      mod.xy.mb <- list()

      for(l in seq_along(micro_pixels)){

        pix <- micro_pixels[[l]]

        tshade   <- pix$tshade
        shadsoil <- pix$shadsoil
        metout   <- pix$metout

        TGRD  <- max(shadsoil$D0cm)
        TSKY  <- max(tshade$TSKYC)
        VEL   <- mean(tshade$VLOC)
        RH    <- mean(tshade$RH)
        QSOLR <- max(tshade$SOLR)

        TA <- metout$TALOC

        mod.t <- lapply(seq_along(TA), function(i){

          NicheMapR::endoR(

            # Group 01
            TA = TA[i]
            , TAREF = TA[i]
            , TGRD = TGRD
            , TSKY = TSKY
            , VEL = VEL
            , RH = RH
            , QSOLR = QSOLR

            # Group 02
            , Z = pars$Z
            , ELEV = pix$elv
            , ABSSB = pars$ABSSB
            , FLTYPE = pars$FLTYPE
            , TCONDSB = TGRD
            , KSUB = pars$KSUB
            , TBUSH = TA[i]
            , BP = pars$BP
            , O2GAS = pars$O2GAS
            , N2GAS = pars$N2GAS
            , CO2GAS = pars$CO2GAS
            , R_PCO2 = pars$CO2GAS / 100
            , PDIF = pars$PDIF

            # ---- BEHAVIOR ----

            # Group 03
            , SHADE = pix$cov
            , FLYHR = pars$FLYHR
            , UNCURL = pars$UNCURL
            , TC_INC = pars$TC_INC
            , PCTWET_INC = pars$PCTWET_INC
            , PCTWET_MAX = pars$PCTWET_MAX
            , AK1_INC = pars$AK1_INC
            , AK1_MAX = pars$AK1_MAX
            , PANT = pars$PANT
            , PANT_INC = pars$PANT_INC
            , PANT_MULT = pars$PANT_MULT

            # ---- MORPHOLOGY ----

            # Geometry
            # Group 04
            , AMASS = pars$AMASS
            , ANDENS = pars$ANDENS
            , SUBQFAT = pars$SUBQFAT
            , FATPCT = pars$FATPCT
            , SHAPE = pars$SHAPE
            , SHAPE_B = pars$SHAPE_B
            , SHAPE_B_MAX = pars$SHAPE_B_MAX
            , SHAPE_C = pars$SHAPE_B
            , PVEN = pars$PVEN
            , PCOND = pars$PCOND
            , SAMODE = pars$SAMODE

            # Fur/Skin properties
            # Group 05
            , FURTHRMK = pars$FURTHRMK
            , DHAIRD = pars$DHAIRD
            , DHAIRV = pars$DHAIRV
            , LHAIRD = pars$LHAIRD
            , LHAIRV = pars$LHAIRV
            , ZFURD_MAX = pars$LHAIRD
            , ZFURV_MAX = pars$LHAIRV
            , ZFURD = pars$ZFURD
            , ZFURV = pars$ZFURV
            , RHOD = pars$RHOD
            , RHOV = pars$RHOV
            , REFLD = pars$REFLD
            , REFLV = pars$REFLV
            , ZFURCOMP = pars$ZFURV
            , KHAIR = pars$KHAIR
            , XR = pars$XR

            # Radiation exchange
            # Group 06
            , EMISAN = pars$EMISAN
            , FABUSH = pars$FABUSH
            , FGDREF = pars$FGDREF
            , FSKREF = pars$FSKREF

            # Thermal physiology
            # Group 07
            , TC = pars$TC
            , TC_MAX = pars$TC_MAX
            , AK1 = pars$AK1
            , AK2 = pars$AK2

            # Evaporation
            # Group 08
            , PCTWET = pars$PCTWET
            , FURWET = pars$FURWET
            , PCTBAREVAP = pars$PCTBAREVAP
            , PCTEYES = pars$PCTEYES
            , DELTAR = pars$DELTAR
            , RELXIT = pars$RELXIT

            # Metabolism / Respiration
            # Group 09
            , QBASAL = (70 * pars$AMASS ^ 0.75) * (4.185 / (24 * 3.6))
            , RQ = pars$RQ
            , EXTREF = pars$EXTREF
            , PANT_MAX = pars$PANT_MAX
            , PZFUR = pars$PZFUR
            , Q10 = pars$Q10
            , TC_MIN = pars$TC_MIN

            # Initial conditions
            # Group 11
            , TS = pars$TC - 3
            , TFA = TA[i]

            # Other model settings
            # Group 12
            , DIFTOL = pars$DIFTOL
            , THERMOREG = pars$THERMOREG
            , RESPIRE = pars$RESPIRE
            , TREGMODE = pars$TREGMODE

            # Default
            , WRITE_INPUT = WRITE_INPUT
          )
        })

        # ------------------------------------------------------------------
        # Extract results using the same logic as the original function

        qgens <- vapply(
          mod.t,
          function(m) {
            as.numeric(m$enbal[1, "ENB"])
          },
          numeric(1)
        )

        evap <- vapply(
          mod.t,
          function(m) {
            sum(
              as.numeric(m$masbal[1, c("H2OResp_g", "H2OCut_g")]),
              na.rm = TRUE
            )
          },
          numeric(1)
        )

        names(qgens) <- NULL
        names(evap) <- NULL

        mod.xy.eb[[l]] <- cbind(
          metout$DOY
          , metout$TIME
          , TA
          , qgens
        )

        mod.xy.mb[[l]] <- cbind(
          metout$DOY
          , metout$TIME
          , TA
          , evap
        )

        colnames(mod.xy.eb[[l]]) <- c("DAY", "time", "TA", "enbal")
        colnames(mod.xy.mb[[l]]) <- c("DAY", "time", "TA", "masbal")
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
  # endonicheR() is run in blocks using the break_n argument.
  #
  # Input:
  #   niche_chunks = a list where each element corresponds to one
  #                  processed block of pixels/locations.
  #
  # Each block contains the same hierarchical structure returned by
  # endo_by_species():
  #
  #   species
  #     energy
  #     evap
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
  # 4. MAIN FLOW OF endonicheR()
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
      "NicheMapR global climate data were not found.\n\n",
      "Please run this once before using endonicheR():\n\n",
      "micronicheR::microniche_setup_global_climate(folder = getwd())\n\n",
      "After the download finishes, restart R and run endonicheR() again."
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
      # Environmental parameters

      , Z = Z
      , ABSSB = ABSSB
      , FLTYPE = FLTYPE
      , KSUB = KSUB
      , BP = BP
      , O2GAS = O2GAS
      , N2GAS = N2GAS
      , CO2GAS = CO2GAS
      , PDIF = PDIF

      # ------------------------------------------------
      # Behavioral parameters

      , FLYHR = FLYHR
      , UNCURL = UNCURL
      , TC_INC = TC_INC
      , PCTWET_INC = PCTWET_INC
      , PCTWET_MAX = PCTWET_MAX
      , AK1_INC = AK1_INC
      , AK1_MAX = AK1_MAX
      , PANT = PANT
      , PANT_INC = PANT_INC
      , PANT_MULT = PANT_MULT

      # ------------------------------------------------
      # Morphological parameters

      # Geometry

      , AMASS = AMASS
      , ANDENS = ANDENS
      , SUBQFAT = SUBQFAT
      , FATPCT = FATPCT
      , SHAPE = SHAPE
      , SHAPE_B = SHAPE_B
      , SHAPE_B_MAX = SHAPE_B_MAX
      , PVEN = PVEN
      , PCOND = PCOND
      , SAMODE = SAMODE
      , ORIENT = ORIENT

      # Fur/Skin properties

      , FURTHRMK = FURTHRMK
      , DHAIRD = DHAIRD
      , DHAIRV = DHAIRV
      , LHAIRD = LHAIRD
      , LHAIRV = LHAIRV
      , ZFURD = ZFURD
      , ZFURV = ZFURV
      , RHOD = RHOD
      , RHOV = RHOV
      , REFLD = REFLD
      , REFLV = REFLV
      , KHAIR = KHAIR
      , XR = XR

      # ------------------------------------------------
      # Radiation exchange

      , EMISAN = EMISAN
      , FABUSH = FABUSH
      , FGDREF = FGDREF
      , FSKREF = FSKREF

      # ------------------------------------------------
      # Thermal physiology

      , TC = TC
      , TC_MAX = TC_MAX
      , AK1 = AK1
      , AK2 = AK2

      # ------------------------------------------------
      # Evaporation

      , PCTWET = PCTWET
      , FURWET = FURWET
      , PCTBAREVAP = PCTBAREVAP
      , PCTEYES = PCTEYES
      , DELTAR = DELTAR
      , RELXIT = RELXIT

      # ------------------------------------------------
      # Metabolism / respiration

      , RQ = RQ
      , EXTREF = EXTREF
      , PANT_MAX = PANT_MAX
      , PZFUR = PZFUR
      , Q10 = Q10
      , TC_MIN = TC_MIN

      # ------------------------------------------------
      # Model settings

      , DIFTOL = DIFTOL
      , THERMOREG = THERMOREG
      , RESPIRE = RESPIRE
      , TREGMODE = TREGMODE

    )

  } else {

    message("traits_df provided. Running endotherm simulations for multiple species.")
  }

  # 4.3. Run microclimate only once per pixel
  # The internal function micro_by_pixel() is called here.

  if (!requireNamespace("NicheMapR", quietly = TRUE)) {
    stop(
      "Package 'NicheMapR' is required.\n",
      "Please install it with:\n",
      "micronicheR::microniche_setup(github = 'mrke/NicheMapR')"
    )
  }

  # Note:
  # NicheMapR must be attached because micro_global()
  # requires internal objects (e.g. "CampNormTbl9_1")
  # that are not available through requireNamespace() alone.

  if (!"package:NicheMapR" %in% search()) {

    suppressPackageStartupMessages(
      require("NicheMapR", character.only = TRUE)
    )
  }

  # 4.4. Run the physiological model for each species
  # The internal function endo_by_species() is called here.

  if (is.null(break_n)) {

    micro_pixels <- micro_by_pixel(df = df)

    niche_list <- endo_by_species(
      micro_pixels = micro_pixels,
      traits_df = traits_df
    )

  } else {

    df_chunks <- split(
      df,
      ceiling(seq_len(nrow(df)) / break_n)
    )

    message(
      "Running endonicheR in ",
      length(df_chunks),
      " block(s) of up to ",
      break_n,
      " pixel(s)."
    )

    niche_chunks <- vector("list", length(df_chunks))

    # =====================================================================
    # >>> V2 CHECKPOINT SYSTEM (BEGIN) 19.06.2026
    # =====================================================================

    for (b in seq_along(df_chunks)) {

      message(
        "Running block ",
        b,
        " of ",
        length(df_chunks),
        " | pixels: ",
        nrow(df_chunks[[b]])
      )

      micro_file <- NULL
      endo_file  <- NULL

      if (!is.null(checkpoint_dir)) {

        micro_file <- file.path(
          checkpoint_dir,
          paste0(
            "micro_block_",
            sprintf("%04d", b),
            ".rds"
          )
        )

        endo_file <- file.path(
          checkpoint_dir,
          paste0(
            "endo_block_",
            sprintf("%04d", b),
            ".rds"
          )
        )

      }

      # ------------------------------------------------------------
      # MICROCLIMATE
      # ------------------------------------------------------------

      if (
        resume &&
        !is.null(micro_file) &&
        file.exists(micro_file)
      ) {

        message(
          "Loading saved microclimate block ",
          b
        )

        micro_pixels <- readRDS(micro_file)

      } else {

        micro_pixels <- micro_by_pixel(
          df = df_chunks[[b]]
        )

        if (!is.null(micro_file)) {

          message("Saving microclimate block ", b)

          saveRDS(
            micro_pixels,
            micro_file
          )

        }

      }

      # ------------------------------------------------------------
      # ENDOR SIMULATION
      # ------------------------------------------------------------

      if (
        resume &&
        !is.null(endo_file) &&
        file.exists(endo_file)
      ) {

        message(
          "Loading saved endoR block ",
          b
        )

        niche_chunks[[b]] <- readRDS(
          endo_file
        )

      } else {

        niche_chunks[[b]] <- endo_by_species(
          micro_pixels = micro_pixels,
          traits_df = traits_df
        )

        if (!is.null(endo_file)) {

          message("Saving endoR block ", b)

          saveRDS(
            niche_chunks[[b]],
            endo_file
          )

        }

      }

    }

    # =====================================================================
    # >>> V2 CHECKPOINT SYSTEM (END)
    # =====================================================================

    niche_list <- merge_niche_lists(niche_chunks)

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
        x = df$x[.data$ID],
        y = df$y[.data$ID],
        ID = df$ID[.data$ID]
      ) %>%
      dplyr::mutate(data_sp = purrr::map(.data$data_sp, tibble::as_tibble)) %>%
      tidyr::unnest(.data$data_sp) %>%
      dplyr::mutate(
        flux = dplyr::coalesce(.data$enbal, .data$masbal)
      ) %>%
      dplyr::select(
        .data$sp, .data$type, .data$ID, .data$x, .data$y,
        .data$DAY, .data$time, .data$TA, .data$flux
      ) %>%
      tidyr::pivot_wider(
        names_from = .data$type,
        values_from = .data$flux
      ) %>%
      dplyr::rename(
        enbal = .data$energy,
        masbal = .data$evap
      )

    if(summary == FALSE){

      return(niche_df)

    } else{

      niche_df_summary <- niche_df %>%
        dplyr::group_by(.data$sp, .data$ID, .data$x, .data$y, .data$DAY) %>%
        dplyr::summarise(
          temp.day.mean = mean(.data$TA, na.rm = TRUE),
          temp.day.sd = stats::sd(.data$TA, na.rm = TRUE),
          enbal.day.sum = sum(.data$enbal, na.rm = TRUE),
          masbal.day.sum = sum(.data$masbal, na.rm = TRUE),
          .groups = "drop"
        )

      return(niche_df_summary)
    }
  }

} # end

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
