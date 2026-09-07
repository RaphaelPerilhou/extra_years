################################################################################
# 19_classification_summary_final.R
#
# Builds a per-county-year classification summary directly from the FINAL
# (best-model-corrected) rasters -- analog of 01bis_classification_dataframes.R
# but for the corrected output. Since Final_*.tif is already in aggregated
# category space, this is a direct freq() count, no crosswalk needed.
#
# Inputs:  outputs/<VERSION>/final_corrected/<statefp>/Final_<year>_<geoid>.tif
#
# Outputs: outputs/<VERSION>/classification_summary_final_parts/<geoid>.csv
#          outputs/<VERSION>/classification_summary_final.csv (long)
#          outputs/<VERSION>/classification_summary_final_wide.csv (wide)
################################################################################

rm(list = ls())
setwd("/users/rperilhou/extra_years")

library(terra)
library(dplyr)
library(tidyr)
library(readr)

VERSION <- "v2"
TARGET_YEARS <- c(2009:2020)

################################################################################
# PATHS
################################################################################
FINAL_DIR <- file.path("outputs", VERSION, "final_corrected")
FINAL_PATH <- function(year, geoid, statefp) {
  file.path(FINAL_DIR, statefp, paste0("Final_", year, "_", geoid, ".tif"))
}

SUMMARY_PARTS_DIR <- file.path("outputs", VERSION, "classification_summary_final_parts")
summary_part_path <- function(geoid) file.path(SUMMARY_PARTS_DIR, paste0(geoid, ".csv"))

SUMMARY_LONG_PATH <- file.path("outputs", VERSION, "classification_summary_final.csv")
SUMMARY_WIDE_PATH <- file.path("outputs", VERSION, "classification_summary_final_wide.csv")

if (!dir.exists(SUMMARY_PARTS_DIR)) dir.create(SUMMARY_PARTS_DIR, recursive = TRUE)

category_labels <- c("0"  = "NonCrop",
                     "1"  = "GM",
                     "2"  = "Tolerant",
                     "3"  = "Vulnerable",
                     "99" = "Unclassified")

################################################################################
# SELECT COUNTIES -- from whatever has final_corrected output
################################################################################

county_lookup <- read.csv("data/county_lookup.csv", stringsAsFactors = FALSE) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)),
         STATEFP = sprintf("%02d", as.integer(STATEFP)))

tasks <- county_lookup %>%
  filter(dir.exists(file.path(FINAL_DIR, STATEFP)))

cat("Counties to process:", nrow(tasks), "\n")

################################################################################
# PER-COUNTY FUNCTION: pixel counts by category for every year
################################################################################

summarize_county_final <- function(geoid, statefp, years) {
  
  out_path <- summary_part_path(geoid)
  if (file.exists(out_path)) {
    cat("  GEOID", geoid, "already done, skipping.\n")
    return(invisible(NULL))
  }
  
  all_rows <- list()
  for (year in years) {
    in_path <- FINAL_PATH(year, geoid, statefp)
    if (!file.exists(in_path)) {
      cat("  Missing final raster, skipping:", in_path, "\n")
      next
    }
    r <- rast(in_path)
    freq_tbl <- terra::freq(r)
    
    all_rows[[as.character(year)]] <- data.frame(
      statefp  = statefp,
      geoid    = geoid,
      year     = year,
      category = category_labels[as.character(freq_tbl$value)],
      n_pixels = freq_tbl$count
    )
  }
  
  if (length(all_rows) == 0) return(invisible(NULL))
  
  county_df <- bind_rows(all_rows)
  write.csv(county_df, out_path, row.names = FALSE)
  cat("  GEOID", geoid, "done, saved to", out_path, "\n")
  invisible(NULL)
}

################################################################################
# RUN (ANUBIS) -- parallelized at county level
################################################################################

library(parallel)
source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c("tasks", "TARGET_YEARS", "summarize_county_final",
                    "FINAL_PATH", "summary_part_path", "category_labels"))

parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  library(terra)
  summarize_county_final(tasks$GEOID[i], tasks$STATEFP[i], TARGET_YEARS)
})

stopCluster(cl)
cat("Per-county summaries complete.\n")

################################################################################
# COMBINE PARTS -> LONG AND WIDE FORMAT (mirrors 01bis's combine step)
################################################################################

part_files <- list.files(SUMMARY_PARTS_DIR, full.names = TRUE)
long <- read_csv(part_files, show_col_types = FALSE)
write_csv(long, SUMMARY_LONG_PATH)

wide <- long %>%
  pivot_wider(
    id_cols = c(statefp, geoid, year),
    names_from = category,
    values_from = n_pixels,
    values_fill = 0
  ) %>%
  mutate(total = NonCrop + GM + Tolerant + Vulnerable + Unclassified) %>%
  select(statefp, geoid, year, NonCrop, GM, Tolerant, Vulnerable, Unclassified, total)

write_csv(wide, SUMMARY_WIDE_PATH)
cat("ALL DONE: wide summary saved to", SUMMARY_WIDE_PATH, "\n")