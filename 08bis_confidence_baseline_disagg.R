################################################################################
# 08bis_confidence_baseline_disagg.R
# Mean confidence per county x crop-type (raw CDL code), rolled up to
# county x commodity level, analogous to 08_confidence_baseline.R's
# county x category level.
#
#   Inputs:  outputs/<VERSION>/class_and_conf_disagg/<STATEFP>/Conf_stacked_disagg_<YYYY>_<GEOID>.tif
#            data/metadata/<VERSION>/missing_geoid_census.csv
#            data/county_lookup.csv
#            data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")
#            outputs/<VERSION>/baseline_bias_disagg_commodity.csv  (from 07bis)

#   Outputs: outputs/<VERSION>/quality_baseline_disagg.csv
#            outputs/<VERSION>/quality_baseline_disagg_commodity.csv
#            outputs/<VERSION>/baseline_measures_disagg_commodity.csv
################################################################################

rm(list = ls())
setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(terra)

TARGET_YEAR <- 2012
VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
COUNTY_LOOKUP_PATH <- "data/county_lookup.csv"
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"

CONF_STACKED_DISAGG_DIR <- file.path("outputs", VERSION, "class_and_conf_disagg")
CONF_STACKED_DISAGG_PATH <- function(geoid, statefp){
  file.path(
    CONF_STACKED_DISAGG_DIR, statefp, paste0(
      "Conf_stacked_disagg_", TARGET_YEAR, "_", geoid, ".tif"))
}

CENSUS_MISSING_PATH <- file.path("data/metadata", VERSION, "missing_geoid_census.csv")

BASELINE_BIAS_DISAGG_PATH         <- file.path("outputs", VERSION, "baseline_bias_disagg_commodity.csv")
QUALITY_BASELINE_DISAG_PATH       <- file.path("outputs", VERSION, "quality_baseline_disagg.csv")
QUALITY_BASELINE_COMMO_DISAG_PATH <- file.path("outputs", VERSION, "quality_baseline_disagg_commodity.csv")
BASELINE_MEASURES_DISAGG_PATH     <- file.path("outputs", VERSION, "baseline_measures_disagg_commodity.csv")

################################################################################
missing_geoid <- read_csv(CENSUS_MISSING_PATH) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  pull(GEOID) %>% unique()

counties <- read.csv(COUNTY_LOOKUP_PATH) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    STATEFP = sprintf("%02d", as.integer(STATEFP))
  ) %>%
  filter(!GEOID %in% missing_geoid) %>%
  select(GEOID, STATEFP)

results <- list()

for (i in seq_len(nrow(counties))) {
  geoid   <- counties$GEOID[i]
  statefp <- counties$STATEFP[i]
  path    <- CONF_STACKED_DISAGG_PATH(geoid, statefp)
  
  if (!file.exists(path)) {
    cat("WARNING: missing disaggregated stacked raster for GEOID", geoid, "-- skipping.\n")
    next
  }
  
  cat("Processing GEOID", geoid, "...\n")
  stacked   <- rast(path)
  code_vals <- values(stacked[[1]])  # raw CDL code band
  conf_vals <- values(stacked[[2]])  # confidence band
  
  is_valid <- !is.na(code_vals)
  mean_conf_county <- mean(conf_vals[is_valid], na.rm = TRUE)
  
  df <- tibble(cdl_code = code_vals[is_valid], confidence = conf_vals[is_valid]) %>%
    group_by(cdl_code) %>%
    summarise(
      n_pixels = n(),
      mean_conf_cdl_code = mean(confidence, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(GEOID = geoid, mean_conf_county = mean_conf_county)
  
  results[[geoid]] <- df
}

quality_baseline_disagg <- bind_rows(results)
write_csv(quality_baseline_disagg, QUALITY_BASELINE_DISAG_PATH)

cat("\nSaved", QUALITY_BASELINE_DISAG_PATH)

################################################################################
# Roll up to commodity level: pixel-weighted mean confidence per
# (GEOID, commodity), since a commodity can span multiple CDL codes.
################################################################################
library(readxl)



code_to_commodity <- read_excel(METADATA_XLSX, sheet = VERSION) %>%
  filter(`Has Census` == 1) %>%
  mutate(`Census Commodity` = trimws(`Census Commodity`)) %>%
  select(`CDL Code`, `Census Commodity`) %>%
  distinct()

quality_baseline_commodity <- quality_baseline_disagg %>%
  left_join(code_to_commodity, by = c("cdl_code" = "CDL Code")) %>%
  filter(!is.na(`Census Commodity`)) %>%
  group_by(GEOID, commodity = `Census Commodity`) %>%
  summarise(
    mean_conf_commodity = weighted.mean(mean_conf_cdl_code, w = n_pixels),
    n_pixels = sum(n_pixels),
    .groups = "drop"
  ) %>%
  left_join(
    quality_baseline_disagg %>% distinct(GEOID, mean_conf_county),
    by = "GEOID"
  )

write_csv(quality_baseline_commodity, QUALITY_BASELINE_COMMO_DISAG_PATH)
quality_baseline_commodity <- read_csv(QUALITY_BASELINE_COMMO_DISAG_PATH)

################################################################################
# Merge into baseline_bias_disagg_commodity.csv from 07bis
################################################################################
bias_disagg <- read_csv(BASELINE_BIAS_DISAGG_PATH)

baseline_disagg <- bias_disagg %>%
  left_join(
    quality_baseline_commodity %>% select(GEOID, commodity, mean_conf_commodity, mean_conf_county),
    by = c("GEOID", "commodity")
  )

write_csv(baseline_disagg, BASELINE_MEASURES_DISAGG_PATH)

cat("Saved", BASELINE_MEASURES_DISAGG_PATH)
