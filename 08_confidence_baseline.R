################################################################################
# 08_confidence_baseline.R
#
#   Inputs:  outputs/<VERSION>/class_and_conf/<STATEFP>/Conf_stacked_<YYYY>_<GEOID>.tif
#            outputs/<VERSION>/baseline_bias.csv
#            data/metadata/<VERSION>/missing_geoid_census.csv
#            data/county_lookup.csv
#
#   Outputs: outputs/<VERSION>/quality_baseline.csv
#            outputs/<VERSION>/baseline_measures.csv
################################################################################
rm(list = ls())
setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(terra)

VERSION <- "v2"
################################################################################
# PATHS BUILDER
################################################################################
COUNTY_LOOKUP_PATH   <- "data/county_lookup.csv"                                      
MISSING_GEOID_PATH   <- file.path("data/metadata", VERSION, "missing_geoid_census.csv")

QUALITY_BASELINE_PATH  <- file.path("outputs", VERSION, "quality_baseline.csv")
BASELINE_BIAS_PATH     <- file.path("outputs", VERSION, "baseline_bias.csv")
BASELINE_MEASURES_PATH <- file.path("outputs", VERSION, "baseline_measures.csv")

CONF_STACKED_DIR     <- file.path("outputs", VERSION, "class_and_conf")
CONF_STACKED_PATH <- function(geoid, statefp){
  file.path(
    CONF_STACKED_DIR, statefp, paste0(
      "Conf_stacked_", TARGET_YEAR, "_", geoid, ".tif"))
}
################################################################################

#Configuration

TARGET_YEAR <- 2012
PIXEL_ACRES <- 0.222395

# Category codes for cropland (excludes 0=NonCrop, 99=Unclassified)
CROPLAND_CODES <- c(1L, 2L, 3L)
CODE_TO_NAME <- c("1" = "GM", "2" = "Tolerant", "3" = "Vulnerable")


#Counties to process

missing_geoid <- read_csv(
  MISSING_GEOID_PATH
) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  pull(GEOID) %>%
  unique()

counties <- read.csv(COUNTY_LOOKUP_PATH) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    STATEFP = sprintf("%02d", as.integer(STATEFP))
  ) %>%
  filter(!GEOID %in% missing_geoid) %>%
  select(GEOID, STATEFP)


#Main loop

results <- list()

for (i in seq_len(nrow(counties))) {
  geoid <- counties$GEOID[i]
  statefp <- counties$STATEFP[i]
  path <- CONF_STACKED_PATH(geoid, statefp)
  
  if (!file.exists(path)) {
    cat("WARNING: missing stacked raster for GEOID", geoid, "— skipping.\n")
    next
  }
  
  cat("Processing GEOID", geoid, "...\n")
  stacked <- rast(path)
  cat_vals <- values(stacked[[1]]) # category band
  conf_vals <- values(stacked[[2]]) # confidence band
  
  # Q_i: county-level quality index
  # Mean confidence of all cropland pixels in the county
  is_cropland <- !is.na(cat_vals) & cat_vals %in% CROPLAND_CODES
  mean_conf_county <- mean(conf_vals[is_cropland], na.rm = TRUE)
  
  # Q_i,c: per-category quality index
  # Mean confidence of cropland pixels per category
  cat_results <- lapply(CROPLAND_CODES, function(code) {
    is_cat <- !is.na(cat_vals) & cat_vals == code
    data.frame(
      GEOID = geoid,
      Category_Code = code,
      Category = CODE_TO_NAME[as.character(code)],
      n_pixels = sum(is_cat),
      mean_conf_category = mean(conf_vals[is_cat], na.rm = TRUE)
    )
  })
  
  cat_df <- bind_rows(cat_results) %>%
    mutate(mean_conf_county = mean_conf_county)
  
  results[[geoid]] <- cat_df
}

# bind output

quality_baseline <- bind_rows(results)
write_csv(quality_baseline, QUALITY_BASELINE_PATH)

# Q_i,c = mean_conf_category / mean_conf_county  (per category per county)
# Q_i   = 1 by construction at baseline (all cropland / all cropland)
# These denominators are needed when we later compute quality for a model:
# Q_model_i,c = mean_conf(reclassified pixels in category c) / mean_conf_category
# Q_model_i   = mean_conf(all reclassified pixels) / mean_conf_county

cat("\nBaseline confidence by county and category:\n")
print(quality_baseline)

#
quality_baseline <- read_csv(QUALITY_BASELINE_PATH)


bias <- read_csv(BASELINE_BIAS_PATH)
baseline <- bias %>%
  select(
    GEOID,
    Category,
    n_suppressed,
    acres_cdl,
    acres_census,
    acres_census_total,
    delta_acres,
    delta_pct_census,
    delta_pct_CDL,
    delta_pct_total,
    cdl_share,
    accuracy_baseline,
    weighted_accuracy_baseline
  ) %>%
  left_join(
    quality_baseline %>%
      select(GEOID, Category, mean_conf_category, mean_conf_county),
    by = c("GEOID", "Category")
  )

write_csv(baseline, BASELINE_MEASURES_PATH)

####
baseline <- read_csv(BASELINE_MEASURES_PATH)
min(baseline$mean_conf_category, na.rm = T)
max(baseline$mean_conf_category, na.rm = T)

min(baseline$mean_conf_county, na.rm = T)
max(baseline$mean_conf_county, na.rm = T)






