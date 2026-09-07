################################################################################
# 18_TM_final.R
#
# Computes category-level transition matrices between consecutive years
# on the FINAL (best-model-corrected) rasters, using terra::crosstab().
# Unlike 02_TM_disag.R, these rasters are already in aggregated category
# space {0,1,2,3,99} -- no disagg->collapse step needed.
#
# Inputs:  outputs/<VERSION>/final_corrected/<statefp>/Final_<year>_<geoid>.tif
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
# TRANSITION FUNCTION: compute 5x5 category TM for one year pair
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
  
  r_from <- rast(path_from)
  r_to   <- rast(path_to)
  
  if (!compareGeom(r_from, r_to, stopOnError = FALSE)) {
    cat("  Rasters not aligned for GEOID", geoid, "years", year_from, "-", year_to, "\n")
    return(NULL)
  }
  
  stacked <- c(r_from, r_to)
  tm      <- crosstab(stacked)
  
  # Relabel raw category codes (0,1,2,3,99) to their names, matching the
  # aggregated output convention from 02_TM_disag.R's collapse step
  rownames(tm) <- category_labels[as.character(rownames(tm))]
  colnames(tm) <- category_labels[as.character(colnames(tm))]
  
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
                    "FINAL_PATH", "TRANSITION_FINAL_PATH", "category_labels"))

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