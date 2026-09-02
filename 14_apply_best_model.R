################################################################################
# 14_apply_best_model.R
#
# Applies each county's SELECTED best model (from best_per_county_q1,
# quality_ratio < 1 constraint) to every year 2009-2020, producing the
# final corrected classification raster per county-year.
#
# The winning model per county is fixed across all years (it was selected
# using only the 2012 evaluation in 100_Reclass_0.R).
#
# Inputs:  outputs/<VERSION>/classified/<STATEFP>/Classified_<YYYY>_<GEOID>.tif
#          data/county_lookup.csv
#          data/SF/counties_2016.rds
#          data/metadata/<VERSION>/reclass_skipped_counties.csv
#          data/NationalCSB_2013-2020_rev23/CSB1320.gdb (for CSB counties)
#          outputs/<VERSION>/post_processing/best_per_county_q1
#
# Outputs: outputs/<VERSION>/final_corrected/<STATEFP>/Final_<YYYY>_<GEOID>.tif
################################################################################

rm(list = ls())
setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(terra)
library(sf)

VERSION <- "v2"
TARGET_YEARS <- 2009:2020
PIXEL_ACRES  <- 0.222395
CROPLAND_CODES <- c(1L, 2L, 3L)

################################################################################
# PATHS
################################################################################
CLASSIFIED_DIR <- file.path("outputs", VERSION, "classified")
CLASSIFIED_PATH <- function(year, geoid, statefp) {
  file.path(CLASSIFIED_DIR, statefp, paste0("Classified_", year, "_", geoid, ".tif"))
}

FINAL_DIR <- file.path("outputs", VERSION, "final_corrected")
FINAL_PATH <- function(year, geoid, statefp) {
  file.path(FINAL_DIR, statefp, paste0("Final_", year, "_", geoid, ".tif"))
}

COUNTY_LOOKUP_PATH <- "data/county_lookup.csv"
SKIPPED_PATH <- file.path("data", "metadata", VERSION, "reclass_skipped_counties.csv")

CSB_GDB   <- "data/NationalCSB_2013-2020_rev23/CSB1320.gdb"
CSB_LAYER <- "national1320"

BEST_MODEL_ASSIGNMENT_PATH <- file.path("outputs", VERSION, "post_processing",
                                        "best_per_county_q1.csv")

################################################################################
# LOAD WINNING MODEL PER COUNTY + PARSE MMU PARAMETERS
################################################################################

best_assignment <- read_csv(BEST_MODEL_ASSIGNMENT_PATH, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    model_type = case_when(
      best_model == "Baseline" ~ "Baseline",
      best_model == "CSB_field_modal" ~ "CSB",
      str_detect(best_model, "^MMU_") ~ "MMU",
      TRUE ~ NA_character_
    ),
    mmu_acres      = if_else(model_type == "MMU",
                             as.numeric(str_extract(best_model, "(?<=MMU_)[0-9.]+(?=ac)")), NA_real_),
    mmu_neighbors  = if_else(model_type == "MMU",
                             as.integer(str_extract(best_model, "(?<=nb)[0-9]+")), NA_integer_),
    mmu_window     = if_else(model_type == "MMU",
                             as.integer(str_extract(best_model, "(?<=w)[0-9]+")), NA_integer_),
    mmu_iter       = if_else(model_type == "MMU",
                             as.integer(str_extract(best_model, "(?<=iter)[0-9]+")), NA_integer_)
  )

cat("Model assignment breakdown:\n")
print(table(best_assignment$model_type, useNA = "ifany"))

################################################################################
# COUNTY LOOKUP, ORDERED BY GEOID (keeps states together for CSB batching)
################################################################################

county_lookup <- read.csv(COUNTY_LOOKUP_PATH, stringsAsFactors = FALSE) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)),
         STATEFP = sprintf("%02d", as.integer(STATEFP))) %>%
  arrange(GEOID) %>%
  left_join(best_assignment %>% select(GEOID, best_model, model_type,
                                       mmu_acres, mmu_neighbors, mmu_window, mmu_iter),
            by = "GEOID")
#exclude the skipped counties of county lookup
skipped_exclusion <- read_csv(SKIPPED_PATH, show_col_types = FALSE)
county_lookup <- county_lookup %>% anti_join(skipped_exclusion, by = "GEOID")

missing_assignment <- sum(is.na(county_lookup$model_type))
if (missing_assignment > 0) {
  cat("WARNING:", missing_assignment, "counties have no model assignment -- ",
      "these will be skipped.\n")
}

counties_sf <- readRDS("data/SF/counties_2016.rds")

for (statefp in unique(county_lookup$STATEFP)) {
  d <- file.path(FINAL_DIR, statefp)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

################################################################################
# CORE APPLICATION FUNCTIONS (adapted from 100_Reclass_0.R, single model only)
################################################################################

apply_mmu <- function(original_category, mmu_acres, neighbors, window, max_fill_iter,
                      cropland_codes = CROPLAND_CODES, pixel_acres = PIXEL_ACRES) {
  output <- terra::deepcopy(original_category)
  for (code in cropland_codes) {
    class_raster <- ifel(original_category == code, 1, NA)
    patch_raster <- terra::patches(class_raster, directions = neighbors, zeroAsNA = TRUE)
    freq_tbl  <- terra::freq(patch_raster)
    small_ids <- freq_tbl$value[freq_tbl$count * pixel_acres < mmu_acres]
    if (length(small_ids) == 0) next
    temp <- ifel(patch_raster %in% small_ids, NA, original_category)
    for (iter in 1:max_fill_iter) {
      holes   <- is.na(temp) & !is.na(original_category)
      n_holes <- global(holes, "sum", na.rm = TRUE)[1, 1]
      if (n_holes == 0) break
      temp <- terra::focal(temp, w = window, fun = "modal", na.rm = TRUE, na.policy = "only")
    }
    output <- ifel(patch_raster %in% small_ids, temp, output)
  }
  ifel(is.na(original_category), NA, output)
}

apply_csb <- function(original_category, county_poly, csb_state,
                      cropland_codes = CROPLAND_CODES) {
  county_poly_proj <- st_transform(county_poly, st_crs(csb_state))
  csb_county <- st_filter(csb_state, county_poly_proj)
  if (nrow(csb_county) == 0) return(NULL)
  
  csb_vect   <- vect(st_transform(csb_county, crs(original_category)))
  modal_cats <- terra::extract(original_category, csb_vect, fun = "modal", na.rm = TRUE)
  csb_county$modal_cat <- as.integer(modal_cats[[names(original_category)]])
  csb_county$modal_cat[is.na(csb_county$modal_cat)] <- 99L
  
  csb_vect_fb  <- vect(st_transform(csb_county, crs(original_category)))
  model_raster <- rasterize(csb_vect_fb, original_category, field = "modal_cat")
  model_raster <- mask(model_raster, original_category)
  ifel(!is.na(original_category) & is.na(model_raster), original_category, model_raster)
}

################################################################################
# PASS 1: BASELINE + MMU counties -- parallelized at COUNTY-YEAR level
################################################################################

baseline_mmu_tasks <- county_lookup %>%
  filter(model_type %in% c("Baseline", "MMU")) %>%
  tidyr::crossing(year = TARGET_YEARS)

cat("\nPass 1 (Baseline + MMU): ", nrow(baseline_mmu_tasks), "county-year tasks\n")

process_baseline_or_mmu <- function(geoid, statefp, year, model_type,
                                    mmu_acres, neighbors, window, iter) {
  out_path <- FINAL_PATH(year, geoid, statefp)
  if (file.exists(out_path)) return(invisible(NULL))
  
  in_path <- CLASSIFIED_PATH(year, geoid, statefp)
  if (!file.exists(in_path)) {
    cat("  Missing classified raster, skipping:", in_path, "\n")
    return(invisible(NULL))
  }
  
  classified <- rast(in_path)
  
  if (model_type == "Baseline") {
    writeRaster(classified, out_path, overwrite = TRUE, datatype = "INT1U")
    return(invisible(NULL))
  }
  
  result <- tryCatch(
    apply_mmu(classified, mmu_acres, neighbors, window, iter),
    error = function(e) {
      cat("  MMU FAILED for", geoid, year, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(result)) {
    writeRaster(as.int(result), out_path, overwrite = TRUE, datatype = "INT1U")
  }
  invisible(NULL)
}

library(parallel)
source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c("baseline_mmu_tasks", "CLASSIFIED_PATH", "FINAL_PATH",
                    "apply_mmu", "process_baseline_or_mmu",
                    "CROPLAND_CODES", "PIXEL_ACRES"))

invisible(parLapplyLB(cl, seq_len(nrow(baseline_mmu_tasks)), function(i) {
  library(terra)
  row <- baseline_mmu_tasks[i, ]
  process_baseline_or_mmu(
    geoid = row$GEOID, statefp = row$STATEFP, year = row$year,
    model_type = row$model_type,
    mmu_acres = row$mmu_acres, neighbors = row$mmu_neighbors,
    window = row$mmu_window, iter = row$mmu_iter
  )
}))

stopCluster(cl)
cat("Pass 1 complete.\n")

################################################################################
# PASS 2: CSB counties -- parallelized at STATE level (load CSB once per state)
################################################################################

csb_counties <- county_lookup %>% filter(model_type == "CSB")
csb_states <- split(csb_counties$GEOID, csb_counties$STATEFP)

cat("\nPass 2 (CSB): ", length(csb_states), "states,",
    nrow(csb_counties), "counties\n")

process_state_csb_final <- function(statefp, geoids_in_state, counties_sf, years) {
  
  # skip the CSB_GDB read entirely if this state's work is already done
  remaining <- expand.grid(geoid = geoids_in_state, year = years, stringsAsFactors = FALSE)
  remaining$done <- file.exists(FINAL_PATH(remaining$year, remaining$geoid, statefp))
  remaining <- remaining[!remaining$done, ]
  if (nrow(remaining) == 0) {
    cat("STATEFP", statefp, "-- all done, skipping.\n")
    return(invisible(NULL))
  }
  
  cat("STATEFP", statefp, "--", nrow(remaining), "county-year tasks remaining\n")
  
  csb_state <- st_read(
    CSB_GDB, layer = CSB_LAYER,
    query = paste0("SELECT * FROM ", CSB_LAYER, " WHERE CSBID LIKE '", statefp, "%'"),
    quiet = TRUE
  )
  
  for (geoid in unique(remaining$geoid)) {
    county_poly <- counties_sf %>% filter(GEOID == geoid)
    for (year in years) {
      out_path <- FINAL_PATH(year, geoid, statefp)
      if (file.exists(out_path)) next
      
      in_path <- CLASSIFIED_PATH(year, geoid, statefp)
      if (!file.exists(in_path)) {
        cat("  Missing classified raster, skipping:", in_path, "\n")
        next
      }
      
      classified <- rast(in_path)
      result <- tryCatch(
        apply_csb(classified, county_poly, csb_state),
        error = function(e) {
          cat("  CSB FAILED for", geoid, year, ":", conditionMessage(e), "\n")
          NULL
        }
      )
      if (!is.null(result)) {
        writeRaster(as.int(result), out_path, overwrite = TRUE, datatype = "INT1U")
      }
    }
  }
  invisible(NULL)
}

cl <- createCluster()

clusterExport(cl, c("csb_states", "counties_sf", "TARGET_YEARS",
                    "CLASSIFIED_PATH", "FINAL_PATH", "apply_csb",
                    "process_state_csb_final", "CSB_GDB", "CSB_LAYER",
                    "CROPLAND_CODES"))

invisible(parLapplyLB(cl, names(csb_states), function(statefp) {
  library(terra)
  library(sf)
  library(dplyr)
  process_state_csb_final(
    statefp = statefp,
    geoids_in_state = csb_states[[statefp]],
    counties_sf = counties_sf,
    years = TARGET_YEARS
  )
}))

stopCluster(cl)
cat("Pass 2 complete.\n")

cat("\nALL DONE: final corrected rasters in", FINAL_DIR, "\n")