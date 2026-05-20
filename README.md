micronicheR
================

`micronicheR` is an R package for microclimate-based physiological niche
modeling using environmental rasters, species traits, and `NicheMapR`.

The package was designed to:

- run microclimate simulations only once per pixel;
- reuse microclimate outputs across multiple species;
- integrate environmental rasters and physiological traits;
- support large-scale spatial simulations.

## Pre-instalation

The `micronicheR` package works correctly with `R version 4.4.3`, which
can be downloaded from the following link: [Download R version
4.4.3](https://cran.r-project.org/bin/windows/base/old/4.4.3/R-4.4.3-win.exe)

After downloading it, change the current R version in your RStudio, as
shown in the figure below:

<img src="man/figures/change-r-version-rstudio.png" width="700px">

## Installation

``` r
install.packages("remotes")

remotes::install_github("jrz-s/micronicheR")
```

## Install dependencies

``` r
microniche_setup(
  cran = c(
    "terra",
    "dplyr",
    "tidyr",
    "tibble",
    "purrr",
    "here"
  ),
  github = "mrke/NicheMapR"
)
```

## Load package

``` r
library(micronicheR)
library(terra)
```

## Example datasets

The package contains example rasters, coordinates, and species traits.
Check the sample data from the package `microniche`.

``` r
list.files(
  system.file("extdata", package = "micronicheR"))
```

    ## [1] "coord_example.csv"     "cover_example.tif"     "elevation_example.tif"
    ## [4] "study_example.tif"     "traits_example.csv"

The following five files should appear:

- “coord_example.csv”<br>
- “cover_example.tif”<br>
- “elevation_example.tif”<br>
- “study_example.tif”<br>
- “traits_example.csv”<br>

### Load example rasters

``` r
cover <- terra::rast(
  system.file(
    "extdata",
    "cover_example.tif",
    package = "micronicheR"
  )
)

elev <- terra::rast(
  system.file(
    "extdata",
    "elevation_example.tif",
    package = "micronicheR"
  )
)

study <- terra::rast(
  system.file(
    "extdata",
    "study_example.tif",
    package = "micronicheR"
  )
)
```

### Load species traits

``` r
traits <- read.csv(
  system.file(
    "extdata",
    "traits_example.csv",
    package = "micronicheR"
  )
)

head(traits)
```

### Load coordinates

``` r
coords <- read.csv(
  system.file(
    "extdata",
    "coord_example.csv",
    package = "micronicheR"
  )
)

head(coords)
```

## Convert rasters to data frame

``` r
df <- rast_to_df(
  rcover = cover,
  rtop = elev,
  rast_or_coord = study,
  return = "data"
)

head(df)
```

## Run microniche model using study raster

``` r
out_df_toraster <- microniche(
    rcover = cover
  , rtop = elev
  , rast_or_coord = study
  , traits_df = traits
  , sample_n = 2
  , list_format = FALSE
  , download_global_climate = TRUE
)

head(out_df_toraster)
```

## Run microniche model using coordinates

``` r
out_df_tocoord <- microniche(
    rcover = cover
  , rtop = elev
  , rast_or_coord = coords
  , traits_df = traits
  , sample_n = 2
  , list_format = FALSE
  , download_global_climate = TRUE
)

head(out_df_tocoord)
```

## Other Outputs

The package supports three output modes.

### 1. List

``` r
out_list <- microniche(
  rcover = cover,
  rtop = elev,
  rast_or_coord = study,
  traits_df = traits,
  sample_n = 2,
  list_format = TRUE,
  download_global_climate = TRUE
)

str(out_list, max.level = 3)
```

Returns hierarchical outputs:

- species
  - energy
  - evaporation
  - pixels

### 2. Data frame

`See` “Run microniche model using study raster” and “Run microniche
model using coordinates”

``` r
list_format = FALSE
```

Returns tidy tabular output.

### 3. Summary statistics

``` r
summary = TRUE
```

``` r
out_summary <- microniche(
    rcover = cover
  , rtop = elev
  , rast_or_coord = coords
  , traits_df = traits
  , sample_n = 2
  , list_format = FALSE
  , download_global_climate = TRUE
  , summary = TRUE
    
)

head(out_summary)
```

Returns daily summaries by species and pixel.

## Disclaimer

`micronicheR` is currently under active development and should be
considered experimental software. Functions, arguments, outputs, and
workflows may change in future versions.

Some functions rely on external resources from `NicheMapR`, especially
the global climate database required by `NicheMapR::micro_global()`.
Downloading these data may require a stable internet connection and can
occasionally fail due to timeout, interrupted downloads, GitHub
availability, or local permission issues.

If the automatic download fails, users should try:

``` r
options(timeout = 600)

NicheMapR::get.global.climate(
  folder = getwd()
)
```

## Author(s)

Zárate-Salazar, J. Rafael, PhD<br> Postdoctoral Researcher<br> Graduate
Program in Ecology and Conservation<br> Federal University of Sergipe
(UFS), São Cristóvão, Sergipe, Brazil<br>

Valença, Saulo E. S., MSc<br> Doctoral Researcher<br> Graduate Program
in Ecology and Conservation<br> Federal University of Sergipe (UFS), São
Cristóvão, Sergipe, Brazil

Senzano Castro, Luis Miguel, PhD<br> Postdoctoral Researcher<br> São
Paulo State University (UNESP), São Paulo, Brazil

Rubalcaba, Juan G., PhD<br> Ramón Y Cajal Researcher<br> Pyrenean
Institute of Ecology, Huesca, Spain

Gouveia, Sidney Feitosa, PhD<br> Professor, Department of Ecology<br>
Federal University of Sergipe, São Cristóvão, Sergipe, Brazil
