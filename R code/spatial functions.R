# ==============================================================================
# TITLE: spatial data processing functions
# DESCRIPTION: A collection of custom functions for regridding spatial arrays, 
#              shifting coordinate systems, and constructing 3D data arrays.
# ==============================================================================

library(tidyverse)
library(akima)

# ------------------------------------------------------------------------------
# FUNCTION 1: Regrid 3D Spatial Arrays via Bilinear Interpolation
# ------------------------------------------------------------------------------
regrid_tp_data <- function(source_data, target_lon, target_lat) {
  
  # 1. Extract source data and dimensions
  source_array <- source_data$array
  source_lon   <- source_data$longitudes
  source_lat   <- source_data$latitudes
  n_months     <- dim(source_array)[3]
  
  # 2. Initialize an empty array for the regridded results
  regridded_array <- array(
    NA,
    dim = c(length(target_lon), length(target_lat), n_months),
    dimnames = list(longitude = as.character(target_lon),
                    latitude  = as.character(target_lat),
                    month     = dimnames(source_array)[[3]]) 
  )
  
  cat("Starting regridding process for", n_months, "months...\n")
  
  # 3. Loop through each month and interpolate the 2D slice
  for (m in 1:n_months) {
    cat("  - Regridding month:", m, "\n")
    
    source_slice_2d <- source_array[, , m]
    
    # Flatten the 2D matrix into a dataframe for coordinate matching
    source_df <- data.frame(
      lon = rep(source_lon, times = length(source_lat)),
      lat = rep(source_lat, each  = length(source_lon)),
      tp  = as.vector(source_slice_2d)
    )
    
    # Remove NA values to prevent akima::interp from failing
    source_df_no_na <- na.omit(source_df)
    
    if (nrow(source_df_no_na) < 4) {
      warning(sprintf("Skipping month %d due to insufficient valid data points.", m))
      next
    }
    
    # 4. Perform bilinear interpolation
    interp_result <- akima::interp(
      x = source_df_no_na$lon,
      y = source_df_no_na$lat,
      z = source_df_no_na$tp,
      xo = target_lon,
      yo = target_lat,
      linear = TRUE,   
      extrap = TRUE    
    )
    
    # 5. Store the interpolated matrix back into the array
    regridded_array[, , m] <- interp_result$z
  }
  
  cat("Regridding complete.\n")
  
  return(list(
    array      = regridded_array,
    longitudes = target_lon,
    latitudes  = target_lat
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION 2: Create 3D Array from Monthly CSV Files
# ------------------------------------------------------------------------------
create_tp_array <- function(yr, dir_press, lon_col = "Longitude", lat_col = "Latitude", tp_col = "TP") {
  
  # Create a list of all 12 file paths and filter existing ones
  file_paths <- map_chr(month.abb, ~ file.path(dir_press, paste0("tp_", .x, yr, ".csv")))
  existing_files <- file_paths %>% keep(file.exists)
  
  if (length(existing_files) == 0) {
    stop("No data files found for the year: ", yr)
  }
  
  # Read all data for the year into a single data frame
  all_year_data <- map_dfr(1:12, function(m) {
    file_path <- file.path(dir_press, paste0("tp_", month.abb[m], yr, ".csv"))
    if (file.exists(file_path)) {
      read_csv(file_path, show_col_types = FALSE) %>%
        mutate(month_index = m)
    }
  })
  
  if (nrow(all_year_data) == 0) stop("Data files were found, but they are all empty.")
  
  # Determine array dimensions and coordinate labels
  lons   <- sort(unique(all_year_data[[lon_col]]))
  lats   <- sort(unique(all_year_data[[lat_col]]))
  months <- 1:12
  
  # Create the empty 3D array with named dimensions
  tp_array <- array(
    NA, 
    dim = c(length(lons), length(lats), length(months)),
    dimnames = list(longitude = as.character(lons), 
                    latitude  = as.character(lats), 
                    month     = month.abb)
  )
  
  # Populate the array
  for (i in 1:nrow(all_year_data)) {
    row <- all_year_data[i, ]
    o_index <- match(row[[lon_col]], lons)
    a_index <- match(row[[lat_col]], lats)
    t_index <- row$month_index
    
    if (!is.na(o_index) && !is.na(a_index) && !is.na(t_index)) {
      tp_array[o_index, a_index, t_index] <- row[[tp_col]]
    }
  }
  
  return(list(
    array      = tp_array,
    longitudes = lons,
    latitudes  = lats
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION 3: Convert Longitudes from [0, 360] to [-180, 180]
# ------------------------------------------------------------------------------
convert_lon_to_180 <- function(data_list) {
  
  # 1. Shift longitudes > 180 into the negative range
  lons_180 <- ifelse(data_list$longitudes > 180, data_list$longitudes - 360, data_list$longitudes)
  
  # 2. Retrieve sorting indices to reorder from -180 to 180
  reorder_idx <- order(lons_180)
  
  # 3. Create the newly sorted longitude vector
  new_lons <- lons_180[reorder_idx]
  
  # 4. Reorder the array data along the longitude dimension (1st dimension)
  new_array <- data_list$array[reorder_idx, , , drop = FALSE]
  
  # 5. Update the dimension names of the new array
  dimnames(new_array)$longitude <- as.character(new_lons)
  
  return(list(
    array      = new_array,
    longitudes = new_lons,
    latitudes  = data_list$latitudes
  ))
}