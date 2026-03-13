# ==============================================================================
# TITLE: OMI/MLS DU to ppb conversion
# DESCRIPTION: calculates ppbv from OMI/MLS Dobson Units (DU).
# ==============================================================================

# 1. LOAD LIBRARIES & CUSTOM FUNCTIONS
library(tidyverse)
library(ncdf4)
library(lubridate)

# Source the spatial functions file
source("spatial_functions.R")

# ==============================================================================
# 2. CONFIGURATION & PATHS
# ==============================================================================
DIR_BASE       <- "./project/OMI-MLSmonthly data"
DIR_PRESS      <- file.path(DIR_BASE, "L3_OMI_MLS_Pressure_1x1")
DIR_DU         <- "omi_mls_2004_2024"
DIR_OUT_PPBV   <- "omi_mls_ppbv"
DIR_SURF_PRESS <- "./project/TCR2/12pm to 2pm monthly surface 700hpa"

setwd(DIR_BASE)

if (!dir.exists(DIR_OUT_PPBV)) {
  dir.create(DIR_OUT_PPBV, recursive = TRUE)
}

# ==============================================================================
# 3. INITIALIZATION
# ==============================================================================
months_abbr <- month.abb

# Load example data to extract target coordinate grid
ex_file <- file.path(DIR_DU, paste0(months_abbr[5], "_2020.csv"))
ex_data <- read_csv(ex_file, show_col_types = FALSE)

target_lon <- sort(unique(ex_data$Longitude))
target_lat <- sort(unique(ex_data$Latitude))

# ==============================================================================
# 4. MAIN PROCESSING LOOP (2005 - 2022)
# ==============================================================================
cat("\n>>> Starting DU to PPBV conversion pipeline...\n")

for (year in 2005:2022) {
  cat(sprintf("\n--- Processing Year: %d ---\n", year))
  
  # A. Process Surface Pressure from NetCDF
  surface_file <- file.path(DIR_SURF_PRESS, paste0("ps_12-14_", year, ".nc"))
  
  surface_p <- nc_open(surface_file)
  lon_nc    <- ncvar_get(surface_p, "lon")
  lat_nc    <- ncvar_get(surface_p, "lat")
  pl_nc     <- ncvar_get(surface_p, "ps")
  nc_close(surface_p)
  
  data_to_convert <- list(array = pl_nc, longitudes = lon_nc, latitudes = lat_nc)
  
  # Shift longitudes and regrid to match target format
  converted_pl_data <- convert_lon_to_180(data_to_convert)
  surface_pl <- regrid_tp_data(converted_pl_data, target_lon, target_lat)
  
  # Fill missing surface pressure values with a default of 1000 hPa
  surface_pl$array[is.na(surface_pl$array)] <- 1000
  
  # B. Iterate over each month for OMI/MLS files
  for (m in 1:12) {
    month_abbr <- months_abbr[m]
    
    du_file    <- file.path(DIR_DU, paste0(month_abbr, "_", year, ".csv"))
    press_file <- file.path(DIR_PRESS, paste0("tp_", month_abbr, year, ".csv"))
    
    if (!file.exists(du_file) || !file.exists(press_file)) next
    
    du_data    <- read_csv(du_file, show_col_types = FALSE)
    press_data <- read_csv(press_file, show_col_types = FALSE)
    
    # Merge DU and Tropopause Pressure by Coordinates
    merged_data <- left_join(du_data, press_data, by = c("Latitude", "Longitude"))
    
    # Extract Surface Pressure for the specific month and format as dataframe
    pressure_matrix_for_month <- surface_pl$array[, , m]
    pressure_df_for_month <- as.data.frame.table(pressure_matrix_for_month, responseName = "SurfacePressure") %>%
      rename(Longitude = longitude, Latitude = latitude) %>%
      mutate(
        Longitude = as.numeric(as.character(Longitude)),
        Latitude  = as.numeric(as.character(Latitude))
      )
    
    # Merge Surface Pressure into the main dataset
    merged_data <- left_join(merged_data, pressure_df_for_month, by = c("Latitude", "Longitude"))
    
    if (!"TP" %in% names(merged_data)) {
      cat(sprintf("WARNING: TP column missing after merge for %s\n", basename(du_file)))
      next
    }
    
    # C. Calculate Volumetric Mixing Ratio (VMR) in ppbv
    # Formula: VMR = 1270 * Ozone_Column / (SurfacePressure - TP)
    merged_data <- merged_data %>%
      mutate(
        VMR = if_else(
          is.na(Ozone_Column) | is.na(SurfacePressure) | is.na(TP), 
          NA_real_,  
          1270 * Ozone_Column / (SurfacePressure - TP)
        )
      ) %>%
      rename(`Ozone_Column(DU)` = Ozone_Column, `VMR(ppbv)` = VMR) %>%
      dplyr::select(Latitude, Longitude, TP, `Ozone_Column(DU)`, `VMR(ppbv)`)
    
    # D. Save Monthly Output
    out_file <- file.path(DIR_OUT_PPBV, paste0(month_abbr, "_", year, "_omi_mls_ppbv.csv"))
    write_csv(merged_data, out_file)
    cat(sprintf("  Saved: %s\n", basename(out_file)))
  }
}
