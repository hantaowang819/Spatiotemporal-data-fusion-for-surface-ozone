# ==============================================================================
# TITLE: Parallelized GEOS-CF hourly to monthly MDA8 Ozone
# DESCRIPTION: Processes GEOS-CF NetCDF to calculate Daily 8-Hour Average ozone. 
# ==============================================================================

# Enable error handling
options(error = quote({dump.frames(to.file=TRUE); q()}))

# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(chron)
library(lattice)
library(ncdf4)
library(mgcv)
library(data.table)
library(fields)
library(tidyverse)
library(doSNOW)     
library(foreach)

# ==============================================================================
# 2. CONFIGURATION & PATHS
# ==============================================================================
MODEL_NAME <- "GEOS-CF"
START_YEAR <- 2021
END_YEAR   <- 2023
NUM_CORES  <- 10

DIR_INPUT  <- "./GEOS-CF/hourly"
DIR_OUTPUT <- "./GEOS-CF/monthly"

setwd(DIR_INPUT)

# ==============================================================================
# 3. INITIALIZATION & HELPER FUNCTIONS
# ==============================================================================
message(">>> Initializing GEOS-CF MDA8 Parallel Pipeline...")

# Read sample NetCDF file to extract static grid metadata
sample_ncfname <- "GEOS-CF.v01.rpl.aqc_tavg_1hr_g1440x721_v1.20171231_daily.nc4"
ncin <- nc_open(sample_ncfname)

variablename <- "O3" 
lon  <- ncvar_get(ncin, "lon")
lat  <- ncvar_get(ncin, "lat")
time <- ncvar_get(ncin, "time", verbose = FALSE)

nlon <- length(lon)  
nlat <- length(lat)  
fill_value <- ncatt_get(ncin, variablename, "_FillValue")$value

nc_close(ncin)

# Helper Function: Convert (Year, Day-of-Year) to "MMDD" string format
day_to_mmdd <- function(year, day) {
  days_in_year <- ifelse(year %% 4 == 0, 366, 365)
  if (day < 1 || day > days_in_year) {
    stop(sprintf("Invalid day_of_year: %d for year: %d", day, year))
  }
  date_obj <- as.Date(paste0(year, "-01-01")) + (day - 1)
  format(date_obj, "%m%d")
}

# ==============================================================================
# 4. PARALLEL CLUSTER SETUP
# ==============================================================================
message(sprintf(">>> Configuring parallel cluster with %d cores...", NUM_CORES))
cl <- makeCluster(NUM_CORES, type = "SOCK")
registerDoSNOW(cl)

# ==============================================================================
# 5. MAIN PROCESSING LOOP (Iterate over Years)
# ==============================================================================
for (year in START_YEAR:END_YEAR) {
  
  ndays <- ifelse(year %% 4 == 0, 366, 365)
  output_filename <- file.path(DIR_OUTPUT, paste0(MODEL_NAME, "-monthly-mda8-", year, ".csv"))
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("--- Processing Year: %d | Total Days: %d ---\n", year, ndays))
  cat(sprintf("==================================================\n"))
  
  # Initialize progress bar for the parallel workers
  pb <- txtProgressBar(max = ndays, style = 3)
  progress_fn <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress_fn)
  
  # --- PARALLEL FOREACH LOOP (Iterate over Days) ---
  day_results <- foreach(day = 1:ndays,
                         .packages = c("ncdf4", "chron"),
                         .options.snow = opts) %dopar% {
                           
                           date_str <- day_to_mmdd(year, day)
                           # Construct strings for cross-day boundary conditions
                           if (day == 1) {
                             pre_date <- "1231"    
                             pre_year <- year - 1
                           } else {
                             pre_date <- day_to_mmdd(year, day - 1)
                             pre_year <- year
                           }
                           
                           if (day == ndays) {
                             post_date <- "0101"   
                             post_year <- year + 1
                           } else {
                             post_date <- day_to_mmdd(year, day + 1)
                             post_year <- year
                           }
                           
                           # Construct dynamic filenames for current, previous, and next days
                           fname_curr <- sprintf("GEOS-CF.v01.rpl.aqc_tavg_1hr_g1440x721_v1.%d%s_daily.nc4", year, date_str)
                           fname_pre  <- sprintf("GEOS-CF.v01.rpl.aqc_tavg_1hr_g1440x721_v1.%d%s_daily.nc4", pre_year, pre_date)
                           fname_post <- sprintf("GEOS-CF.v01.rpl.aqc_tavg_1hr_g1440x721_v1.%d%s_daily.nc4", post_year, post_date)
                           
                           # Load daily arrays into memory
                           nc_curr <- nc_open(fname_curr); curr_data <- ncvar_get(nc_curr, variablename); nc_close(nc_curr)
                           nc_pre  <- nc_open(fname_pre);  pre_data  <- ncvar_get(nc_pre, variablename);  nc_close(nc_pre)
                           nc_post <- nc_open(fname_post); post_data <- ncvar_get(nc_post, variablename); nc_close(nc_post)
                           
                           # Initialize matrices for daily results (nlon x nlat)
                           world8hr_day <- matrix(0, nrow = nlon, ncol = nlat)  # MDA8 Ozone (ppb)
                           worldhr_day  <- matrix(0, nrow = nlon, ncol = nlat)  # Starting hour of the maximum window
                           
                           # Iterate across the spatial grid
                           for (i in 1:nlon) {
                             
                             # Calculate local time offset based on longitude
                             if (lon[i] >= 0) {
                               xloctm <- lon[i] / 15
                               offset <- -floor(xloctm)
                             } else {
                               xloctm <- (lon[i] + 360) / 15
                               offset <- 24 - floor(xloctm)
                             }
                             
                             for (j in 1:nlat) {
                               max_val <- -Inf
                               maxhr   <- NA
                               
                               # Scan candidate starting hours (7:00 to 23:00 local time)
                               for (hr in 7:23) {
                                 m <- offset + hr   
                                 window_vals <- numeric(8)
                                 
                                 # Construct the 8-hour rolling window
                                 for (h in 0:7) {
                                   target <- m + h
                                   
                                   # Handle cross-day boundaries
                                   if (target < 1) {
                                     tmp <- pre_data[i, j, target + 24]
                                   } else if (target > 24) {
                                     tmp <- post_data[i, j, target - 24]
                                   } else {
                                     tmp <- curr_data[i, j, target]
                                   }
                                   
                                   # Apply missing value mask
                                   if (!is.na(tmp) && tmp == fill_value) {
                                     tmp <- NA
                                   }
                                   window_vals[h + 1] <- tmp
                                 }
                                 
                                 # Calculate rolling average and update maximum
                                 avg_8hr <- mean(window_vals, na.rm = TRUE)
                                 if (!is.na(avg_8hr) && avg_8hr > max_val) {
                                   max_val <- avg_8hr
                                   maxhr   <- hr
                                 }
                               } 
                               
                               world8hr_day[i, j] <- max_val
                               worldhr_day[i, j]  <- maxhr
                             } 
                           } 
                           
                           # Return output for the parallel backend
                           list(world8hr_day = world8hr_day, worldhr_day = worldhr_day)
                         } 
  
  close(pb)
  
  # --- Reconstruct the 3D Array from Parallel Workers ---
  world8hr <- array(0, dim = c(nlon, nlat, ndays))
  worldhr  <- array(0, dim = c(nlon, nlat, ndays))
  
  for (d in 1:ndays) {
    world8hr[,,d] <- day_results[[d]]$world8hr_day
    worldhr[,,d]  <- day_results[[d]]$worldhr_day
  }
  
  # ==============================================================================
  # 6. MONTHLY AGGREGATION & EXPORT
  # ==============================================================================
  cat("\n  -> Aggregating daily MDA8 into monthly averages...\n")
  
  # Calculate monthly means accounting for leap years
  if (year %% 4 == 0) {  
    su8hr01 <- apply(world8hr[,,1:31],   c(1,2), mean)
    su8hr02 <- apply(world8hr[,,32:60],  c(1,2), mean)
    su8hr03 <- apply(world8hr[,,61:91],  c(1,2), mean)
    su8hr04 <- apply(world8hr[,,92:121], c(1,2), mean)
    su8hr05 <- apply(world8hr[,,122:152],c(1,2), mean)
    su8hr06 <- apply(world8hr[,,153:182],c(1,2), mean)
    su8hr07 <- apply(world8hr[,,183:213],c(1,2), mean)
    su8hr08 <- apply(world8hr[,,214:244],c(1,2), mean)
    su8hr09 <- apply(world8hr[,,245:274],c(1,2), mean)
    su8hr10 <- apply(world8hr[,,275:305],c(1,2), mean)
    su8hr11 <- apply(world8hr[,,306:336],c(1,2), mean)
    su8hr12 <- apply(world8hr[,,337:366],c(1,2), mean)
  } else {  
    su8hr01 <- apply(world8hr[,,1:31],   c(1,2), mean)
    su8hr02 <- apply(world8hr[,,32:59],  c(1,2), mean)
    su8hr03 <- apply(world8hr[,,60:90],  c(1,2), mean)
    su8hr04 <- apply(world8hr[,,91:120], c(1,2), mean)
    su8hr05 <- apply(world8hr[,,121:151],c(1,2), mean)
    su8hr06 <- apply(world8hr[,,152:181],c(1,2), mean)
    su8hr07 <- apply(world8hr[,,182:212],c(1,2), mean)
    su8hr08 <- apply(world8hr[,,213:243],c(1,2), mean)
    su8hr09 <- apply(world8hr[,,244:273],c(1,2), mean)
    su8hr10 <- apply(world8hr[,,274:304],c(1,2), mean)
    su8hr11 <- apply(world8hr[,,305:335],c(1,2), mean)
    su8hr12 <- apply(world8hr[,,336:365],c(1,2), mean)
  }
  
  # Format output and convert units (mol/mol to ppbv)
  loc <- make.surface.grid(list(x = lon, y = lat))
  
  df <- data.frame(
    lon = loc[,1], lat = loc[,2], 
    DMA8_1  = as.vector(su8hr01) * 10^9, 
    DMA8_2  = as.vector(su8hr02) * 10^9, 
    DMA8_3  = as.vector(su8hr03) * 10^9,
    DMA8_4  = as.vector(su8hr04) * 10^9, 
    DMA8_5  = as.vector(su8hr05) * 10^9, 
    DMA8_6  = as.vector(su8hr06) * 10^9,
    DMA8_7  = as.vector(su8hr07) * 10^9, 
    DMA8_8  = as.vector(su8hr08) * 10^9, 
    DMA8_9  = as.vector(su8hr09) * 10^9,
    DMA8_10 = as.vector(su8hr10) * 10^9, 
    DMA8_11 = as.vector(su8hr11) * 10^9, 
    DMA8_12 = as.vector(su8hr12) * 10^9
  )
  
  write.csv(df, output_filename, row.names = FALSE)
  
  cat(sprintf("  -> Export successful: %s\n", basename(output_filename)))
  
  rm(day_results, world8hr, worldhr, su8hr01, su8hr02, su8hr03, su8hr04, 
     su8hr05, su8hr06, su8hr07, su8hr08, su8hr09, su8hr10, su8hr11, su8hr12, df)
  gc()
}

stopCluster(cl)