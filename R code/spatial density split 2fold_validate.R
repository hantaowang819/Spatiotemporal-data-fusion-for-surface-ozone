# ==============================================================================
# TITLE: Spatial 2fold-validation: cluster vs. sparse
# DESCRIPTION: Implements a density-based spatial cross-validation framework. 
#              Splits ground observation sites into spatially dense (train) and 
#              sparse (validation) sets using k-NN density to rigorously evaluate 
#              the spatial generalization of data fusion models
# ==============================================================================

library(data.table)
library(ggplot2)
library(tidyverse)

# Source the fRAMP functions file
source("fRAMP_with_additive.R")

# ==============================================================================
# 1. CONFIGURATION & PATHS
# ==============================================================================
# Set working paths for collocation and grid files
PATH_OBS_SAT  <- "./data/ready_to_ramp/toar_omi_tcol_2005_2022.csv"
PATH_GRID_SAT <- "./data/ready_to_ramp/omi_tcol_2005_2022.csv"
PATH_OBS_M3   <- "./data/ready_to_ramp/toar_M3fusion_1990_2022.csv"
PATH_GRID_M3  <- "./data/ready_to_ramp/M3fusion_1990_2022.csv"

# Define output directory
DIR_OUT <- "./data/kfold_validation_files/cluster_vs_sparse/OMI_M3"
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# Clustering parameters
THRESHOLD_DEG <- 0.5   # Search radius in degrees (~55km)
MIN_NEIGHBORS <- 3     # Density threshold (>=3: Train/Cluster, <3: Val/Sparse)
TARGET_YEARS  <- 2005:2022

# ==============================================================================
# 2. LOAD & SYNC DATASETS
# ==============================================================================
message(">>> Loading and synchronizing datasets...")

obs_sat <- fread(PATH_OBS_SAT)[year %in% TARGET_YEARS]
obs_m3  <- fread(PATH_OBS_M3)[year %in% TARGET_YEARS]

# Standardize grid columns
cols_keep <- c("lat", "lon", paste0("Model_", TARGET_YEARS))

grid_sat <- fread(PATH_GRID_SAT)[, ..cols_keep]
fwrite(grid_sat, file.path(DIR_OUT, "Sat_Grid_Filtered.csv"))

# Generate unique spatiotemporal keys for exact matching
obs_sat[, match_key := paste(lat_toar, lon_toar, year, sep = "_")]
obs_m3[, match_key := paste(lat_toar, lon_toar, year, sep = "_")]

# Remove internal duplicates
obs_sat <- unique(obs_sat, by = "match_key")
obs_m3  <- unique(obs_m3, by = "match_key")

# Enforce strict intersection to ensure models are evaluated on identical observations
common_keys <- intersect(obs_sat$match_key, obs_m3$match_key)
obs_sat <- obs_sat[match_key %in% common_keys]
obs_m3  <- obs_m3[match_key %in% common_keys]

cat(sprintf("    Sync complete | Sat Rows: %d | M3 Rows: %d | Common Intersect: %d\n", 
            nrow(obs_sat), nrow(obs_m3), length(common_keys)))

# ==============================================================================
# 3. SPATIAL CLUSTERING LOGIC (Density-Based Splitting)
# ==============================================================================
assignment_file <- file.path(DIR_OUT, "cluster_sparse_assignments.csv")

if (file.exists(assignment_file)) {
  message(">>> Loading existing Cluster/Sparse assignments...")
  unique_locs <- fread(assignment_file)
} else {
  message(">>> Calculating spatial density for Train/Val split...")
  
  # Extract unique coordinates
  unique_locs <- unique(obs_sat[, .(lat_toar, lon_toar)])
  
  # Compute Euclidean distance matrix
  coords_mat <- as.matrix(unique_locs[, .(lat_toar, lon_toar)])
  dist_mat   <- as.matrix(dist(coords_mat)) 
  
  # Count neighbors within search radius (subtract 1 to exclude self)
  neighbor_counts <- colSums(dist_mat < THRESHOLD_DEG) - 1
  
  # Apply density mask
  unique_locs[, is_clustered := neighbor_counts >= MIN_NEIGHBORS]
  unique_locs[, neighbor_count := neighbor_counts]
  
  fwrite(unique_locs, assignment_file)
}

# Report spatial split statistics
n_train <- sum(unique_locs$is_clustered)
n_val   <- sum(!unique_locs$is_clustered)

cat(sprintf("    Split Summary | Total Stations: %d | TRAIN (Clustered): %d | VAL (Sparse): %d\n", 
            nrow(unique_locs), n_train, n_val))

# Merge cluster assignments back to main datasets
obs_sat <- merge(obs_sat, unique_locs, by = c("lat_toar", "lon_toar"), all.x = TRUE)

key_cluster_map <- obs_sat[, .(match_key, is_clustered)]
obs_m3 <- merge(obs_m3, key_cluster_map, by = "match_key")

# ==============================================================================
# 4. MODEL PIPELINE EXECUTION & EVALUATION
# ==============================================================================

# Helper function to extract, merge, and evaluate predictions
evaluate_model <- function(model_name, obs_data, grid_path, output_dir) {
  
  message(sprintf("\n>>> Executing %s Pipeline (Cluster Train -> Sparse Val)...", toupper(model_name)))
  
  path_train <- file.path(output_dir, sprintf("train_%s_cluster.csv", model_name))
  path_val   <- file.path(output_dir, sprintf("verify_%s_sparse.csv", model_name))
  
  # Apply spatial split
  df_train <- obs_data[is_clustered == TRUE]
  df_val   <- obs_data[is_clustered == FALSE]
  
  fwrite(df_train, path_train)
  fwrite(df_val, path_val)
  
  # Execute RAMP correction (silence verbose outputs)
  capture.output({
    get_ramp_correction(
      collocated_df_path = path_train,
      model_grid_df_path = grid_path,
      results_dir = output_dir, 
      plotting = FALSE 
    )
  })
  
  # Handle output artifacts
  default_out <- file.path(output_dir, "ramp_yearly_combined_results.csv")
  final_out   <- file.path(output_dir, sprintf("result_%s_sparse_val.csv", model_name))
  if (file.exists(default_out)) file.rename(default_out, final_out)
  
  # Process and reshape model predictions
  message("    Matching predictions with observations...")
  grid_res <- fread(final_out)
  model_cols <- grep("^Model_[0-9]{4}$", names(grid_res), value = TRUE)
  years_found <- sort(as.integer(gsub("Model_", "", model_cols)))
  
  target_l1 <- paste0("L1_Model_", years_found)
  target_l2 <- paste0("L2_Model_", years_found)
  valid_idx <- which(target_l2 %in% names(grid_res))
  
  res_selected <- grid_res %>%
    dplyr::select(lon, lat, all_of(target_l1[valid_idx]), all_of(target_l2[valid_idx])) %>%
    pivot_longer(cols = -c(lon, lat), names_to = c(".value", "year"), names_pattern = "(.*)_(\\d{4})") %>%
    rename(predicted_L1 = L1_Model, predicted_L2 = L2_Model) %>%
    mutate(year = as.integer(year)) %>%
    as.data.table()
  
  # Coordinate rounding to ensure exact floating-point matching
  df_val[, c("lon", "lat") := .(round(lon, 4), round(lat, 4))]
  res_selected[, c("lon", "lat") := .(round(lon, 4), round(lat, 4))]
  
  eval_df <- merge(df_val, res_selected, by = c("lon", "lat", "year"), all.x = FALSE)
  eval_df <- eval_df[!is.na(predicted_L1) & !is.na(predicted_L2)]
  
  # Compute error metrics
  eval_df[, deviation := observed_ozone - predicted_L1]
  eval_df[, z_score   := deviation / sqrt(predicted_L2)]
  
  lm_mod <- lm(observed_ozone ~ predicted_L1, data = eval_df)
  r2_val <- summary(lm_mod)$r.squared
  
  summary_stats <- data.table(
    n_obs  = nrow(eval_df),
    z_mean = mean(eval_df$z_score, na.rm = TRUE),
    z_sd   = sd(eval_df$z_score, na.rm = TRUE),
    rmse   = sqrt(mean(eval_df$deviation^2, na.rm = TRUE)),
    r2     = r2_val  
  )
  
  # Export results
  fwrite(eval_df, file.path(output_dir, sprintf("final_eval_%s_sparse.csv", model_name)))
  fwrite(summary_stats, file.path(output_dir, sprintf("summary_stats_%s.csv", model_name)))
  
  cat(sprintf("    Pipeline complete. RMSE: %.3f | R2: %.3f\n", summary_stats$rmse, summary_stats$r2))
}

# --- Execute Part A: Satellite Model ---
evaluate_model("sat", obs_sat, file.path(DIR_OUT, "Sat_Grid_Filtered.csv"), DIR_OUT)

# --- Execute Part B: M3Fusion Model ---
# Filter and export M3 grid once
grid_m3 <- fread(PATH_GRID_M3)[, ..cols_keep]
path_grid_m3_filtered <- file.path(DIR_OUT, "M3_Grid_Filtered.csv")
fwrite(grid_m3, path_grid_m3_filtered)

evaluate_model("m3", obs_m3, path_grid_m3_filtered, DIR_OUT)

message("\n>>> All spatial cross-validation tasks completed successfully.")