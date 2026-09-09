################################################################################
# 16_TM_final.R
#
# Computes category-level transition matrices between consecutive years
# on the FINAL (best-model-corrected) rasters, using terra::crosstab().
# Unlike 02_TM_disag.R, these rasters are already in aggregated category
# space {0,1,2,3,99} -- no disagg->collapse step needed.
#
# Both years are masked to the county's UNION MASK before crosstab (the
# union mask is identical for every year of a given county, since
# 01_mask_disag.R builds it once across all study years). Off-mask
# pixels are dropped entirely before tabulation -- never counted at all,
# not even as a "NoData" category. This means any NA remaining after
# masking can ONLY be genuine model-induced data loss (MMU's incomplete
# hole-filling after max_fill_iter iterations; Baseline and CSB never
# introduce new in-mask NA), never off-mask background. Without this
# masking step, off-mask pixels would flood every county's "NoData"
# category and swamp any genuine model-induced signal, since off-mask
# area is typically much larger than the agricultural mask itself.
#
# Inputs:  outputs/<VERSION>/final_corrected/<statefp>/Final_<year>_<geoid>.tif
#          outputs/<VERSION>/classified/<statefp>/Classified_<year>_<geoid>.tif
#          (used only as the union-mask reference; any one year works
#          since the mask is identical across years for a given county)
#
# Outputs: outputs/<VERSION>/transition_final/<statefp>/TM_<year_from><year_to>_<geoid>.csv
################################################################################

rm(list = ls())
setwd("/users/rperilhou/extra_years")

library(terra)
library(dplyr)

VERSION <- "v2"
TARGET_YEARS <- c(2009:2020)

################################################################################
# PATHS
################################################################################
FINAL_DIR <- file.path("outputs", VERSION, "final_corrected")
FINAL_PATH <- function(year, geoid, statefp) {
  file.path(FINAL_DIR, statefp, paste0("Final_", year, "_", geoid, ".tif"))
}

CLASSIFIED_DIR <- file.path("outputs", VERSION, "classified")
CLASSIFIED_PATH <- function(year, geoid, statefp) {
  file.path(CLASSIFIED_DIR, statefp, paste0("Classified_", year, "_", geoid, ".tif"))
}

TRANSITION_FINAL_DIR <- file.path("outputs", VERSION, "transition_final")
TRANSITION_FINAL_PATH <- function(year_from, year_to, geoid, statefp) {
  file.path(TRANSITION_FINAL_DIR, statefp,
            paste0("TM_", year_from, year_to, "_", geoid, ".csv"))
}

category_labels <- c("0"  = "NonCrop",
                     "1"  = "GM",
                     "2"  = "Tolerant",
                     "3"  = "Vulnerable",
                     "99" = "Unclassified")

# Helper: map raw dimnames (character codes, or NA/"NA" for no-data) to
# readable labels. After masking to the union mask, "NoData" can only
# mean genuine model-induced NA (never off-mask background).
relabel_dimnames <- function(x) {
  ifelse(
    is.na(x) | x == "NA",
    "NoData",
    category_labels[as.character(x)]
  )
}

################################################################################
# SELECT COUNTIES -- from whatever counties actually have final_corrected output
# (naturally excludes skipped counties, since those never got a Final_*.tif)
################################################################################

county_lookup <- read.csv("data/county_lookup.csv", stringsAsFactors = FALSE) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)),
         STATEFP = sprintf("%02d", as.integer(STATEFP)))

tasks <- county_lookup %>%
  filter(dir.exists(file.path(FINAL_DIR, STATEFP)))

cat("Counties to process:", nrow(tasks), "\n")

for (statefp in unique(tasks$STATEFP)) {
  d <- file.path(TRANSITION_FINAL_DIR, statefp)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

################################################################################
# TRANSITION FUNCTION: compute 5x5(+NoData) category TM for one year pair,
# restricted to the county's union mask extent
################################################################################

compute_transition_final <- function(year_from, year_to, geoid, statefp, force = FALSE) {
  
  out_path <- TRANSITION_FINAL_PATH(year_from, year_to, geoid, statefp)
  
  if (file.exists(out_path) && !force) {
    cat("  Already exists, skipping:", out_path, "\n")
    return(as.matrix(read.csv(out_path, row.names = 1, check.names = FALSE)))
  }
  
  path_from <- FINAL_PATH(year_from, geoid, statefp)
  path_to   <- FINAL_PATH(year_to,   geoid, statefp)
  
  if (!file.exists(path_from)) {
    cat("  Missing final raster, skipping:", path_from, "\n")
    return(NULL)
  }
  if (!file.exists(path_to)) {
    cat("  Missing final raster, skipping:", path_to, "\n")
    return(NULL)
  }
  
  # Union mask reference: identical every year for this county, so any
  # single year's Classified raster carries the correct extent.
  ref_path <- CLASSIFIED_PATH(year_from, geoid, statefp)
  if (!file.exists(ref_path)) {
    cat("  Missing classified reference raster, skipping:", ref_path, "\n")
    return(NULL)
  }
  union_mask_ref <- rast(ref_path)
  
  r_from <- mask(rast(path_from), union_mask_ref)
  r_to   <- mask(rast(path_to),   union_mask_ref)
  
  if (!compareGeom(r_from, r_to, stopOnError = FALSE)) {
    cat("  Rasters not aligned for GEOID", geoid, "years", year_from, "-", year_to, "\n")
    return(NULL)
  }
  
  stacked <- c(r_from, r_to)
  tm      <- crosstab(stacked, useNA = TRUE)
  
  # Relabel raw category codes (0,1,2,3,99) AND the NA/"NA" no-data
  # category to readable labels -- nothing gets dropped. NoData here is
  # guaranteed genuine (MMU fill-failure), never off-mask background.
  rownames(tm) <- relabel_dimnames(rownames(tm))
  colnames(tm) <- relabel_dimnames(colnames(tm))
  
  write.csv(tm, out_path)
  cat("  Saved:", out_path, "\n")
  
  return(tm)
}

################################################################################
# RUN (ANUBIS)
################################################################################

library(parallel)
source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c("tasks", "TARGET_YEARS", "compute_transition_final",
                    "FINAL_PATH", "CLASSIFIED_PATH", "TRANSITION_FINAL_PATH",
                    "category_labels", "relabel_dimnames"))

parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  library(terra)
  geoid   <- tasks$GEOID[i]
  statefp <- tasks$STATEFP[i]
  for (j in seq_len(length(TARGET_YEARS) - 1)) {
    compute_transition_final(TARGET_YEARS[j], TARGET_YEARS[j + 1], geoid, statefp)
  }
})

stopCluster(cl)
cat("ALL DONE: Final-model transition matrices saved in", TRANSITION_FINAL_DIR, "\n")