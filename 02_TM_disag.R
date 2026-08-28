################################################################################
# 02_TM_disag.R
# Disaggregated variant of 02_transition_matrix.R
#
# Computes pixel-level transition matrices between consecutive years
# at the RAW CDL CODE level using terra::crosstab(). Then, we derive the
# aggregated 5x5 transition matrix from it by collapsing codes to
# categories. 
#
# Inputs:  outputs/<VERSION>/classified_disagg/<statefp>/Disagg_<year>_<geoid>.tif
#          data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")

# Outputs: outputs/<VERSION>/transition_disagg/<statefp>/TMD_<year_from><year_to>_<geoid>.csv
#          outputs/<VERSION>/transitions/<statefp>/TM_<year_from><year_to>_<geoid>.csv
#
# Remark: the disaggregated matrices dimensions vary by
# county and year pair: rows = CDL codes present in year t, columns =
# CDL codes present in year t+1. Row/column names are the CDL codes
# themselves.
################################################################################
################################################################################
# CONFIGURATION
################################################################################
rm(list = ls())
setwd("/users/rperilhou/extra_years")

library(terra)
library(dplyr)
library(readxl)

VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"
STATES_SF_PATH <- "data/SF/states_2016.rds"
COUNTIES_SF_PATH <- "data/SF/counties_2016.rds"

DISAGG_DIR <- file.path("outputs", VERSION, "classified_disagg")
DISAGG_PATH <- function(year, geoid, statefp) {
  file.path(DISAGG_DIR, statefp,
            paste0("Disagg_", year, "_", geoid, ".tif"))
}

TRANSITION_DISAGG_DIR <- file.path("outputs", VERSION, "transition_disagg")
TRANSITION_DISAGG_PATH <- function(year_from, year_to, geoid, statefp) {
  file.path(TRANSITION_DISAGG_DIR, statefp,
            paste0("TMD_", year_from, year_to, "_", geoid, ".csv"))
}

# Path to the aggregated 5x5 matrices
TRANSITION_AGG_DIR <- file.path("outputs", VERSION, "transition")
TRANSITION_AGG_PATH <- function(year_from, year_to, geoid, statefp) {
  file.path(TRANSITION_AGG_DIR, statefp,
            paste0("TM_", year_from, year_to, "_", geoid, ".csv"))
}

################################################################################
#select counties and years
states_sf   <- readRDS(STATES_SF_PATH)
counties_sf <- readRDS(COUNTIES_SF_PATH)

# All 48 contiguous states (STATEFP codes)
contiguous_statefps <- states_sf %>%
  sf::st_drop_geometry() %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  pull(STATEFP)

cat("Contiguous states available:", length(contiguous_statefps), "\n")

TARGET_STATEFPS <- contiguous_statefps
TARGET_YEARS    <- c(2009:2020)


# Create output directories
for (statefp in TARGET_STATEFPS) {
  d_disagg <- file.path(TRANSITION_DISAGG_DIR, statefp)
  d_agg    <- file.path(TRANSITION_AGG_DIR, statefp)
  if (!dir.exists(d_disagg)) dir.create(d_disagg, recursive = TRUE)
  if (!dir.exists(d_agg))    dir.create(d_agg,    recursive = TRUE)
}

################################################################################
# CODE -> CATEGORY MAPPING
# Identical to reclass_table in the aggregated scripts. Any code
# absent from this map corresponds to 99 (Unclassified).
################################################################################

cdl_map <- read_excel(METADATA_XLSX, sheet = VERSION) %>%
  filter(Keep == 1)

code_to_category <- setNames(cdl_map$`Category Code`, as.character(cdl_map$`CDL Code`))


category_labels <- c("0"  = "NonCrop",
                     "1"  = "GM",
                     "2"  = "Tolerant",
                     "3"  = "Vulnerable",
                     "99" = "Unclassified")

################################################################################
# COLLAPSE DISAGGREGATED TM TO 5x5 AGGREGATED TM
# Replaces the old 01_mask_and_classify.R + 02_transition_matrix.R output.
################################################################################
collapse_to_aggregated <- function(tmd, year_from, year_to, geoid, statefp) {
  
  out_path <- TRANSITION_AGG_PATH(year_from, year_to, geoid, statefp)
  
  if (file.exists(out_path)) {
    cat("  Already exists, skipping:", out_path, "\n")
    return(as.matrix(read.csv(out_path, row.names = 1, check.names = FALSE)))
  }
  
  # Map each CDL code to its category; unmapped codes -> 99
  map_codes <- function(codes) {
    cats <- code_to_category[codes]
    cats[is.na(cats)] <- 99
    category_labels[as.character(cats)]
  }
  
  row_cats <- map_codes(rownames(tmd))
  col_cats <- map_codes(colnames(tmd))
  
  # Collapse rows then columns by category
  collapsed <- rowsum(tmd, group = row_cats)
  collapsed <- t(rowsum(t(collapsed), group = col_cats))
  
  write.csv(collapsed, out_path)
  cat("  Saved aggregated TM:", out_path, "\n")
  
  return(collapsed)
}

################################################################################
# TRANSITION FUNCTION: Compute disaggregated TM for one year pair
################################################################################
compute_transition_disagg <- function(year_from, year_to, geoid, statefp, force = FALSE) {
  
  out_path <- TRANSITION_DISAGG_PATH(year_from, year_to, geoid, statefp)
  
  if (file.exists(out_path) && !force) {
    cat("  Already exists, skipping:", out_path, "\n")
    return(as.matrix(read.csv(out_path, row.names = 1, check.names = FALSE)))
  }
  
  # Check inputs exist
  path_from <- DISAGG_PATH(year_from, geoid, statefp)
  path_to   <- DISAGG_PATH(year_to,   geoid, statefp)
  

  
  if (!file.exists(path_from)) {
    stop("Missing disaggregated raster: ", path_from,
         "\nRun 01_mask_and_keep_disagg.R first.")
  }
  if (!file.exists(path_to)) {
    stop("Missing disaggregated raster: ", path_to,
         "\nRun 01_mask_and_keep_disagg.R first.")
  }
  
  # Load both rasters
  r_from <- rast(path_from)
  r_to   <- rast(path_to)
  
  # Verify alignment (same extent, resolution and CRS)
  if (!compareGeom(r_from, r_to, stopOnError = FALSE)) {
    stop("Rasters are not aligned for GEOID ", geoid,
         " years ", year_from, "-", year_to)
  }
  
  # Stack and crosstab: counts every (code in t, code in t+1) combination.
  # NAs (outside union mask) are ignored. In-mask pixels always carry a
  # raw CDL code, so the total equals the union mask size.
  cat("  Running crosstab", year_from, "->", year_to, "...\n")
  stacked <- c(r_from, r_to)
  tm      <- crosstab(stacked)
  
  # Row/column names are the raw CDL codes (characters). Dimensions vary
  # by county and year pair — only codes actually present appear.
  total <- sum(tm)
  cat("  Dimensions:", nrow(tm), "x", ncol(tm),
      "| Total transition pixels:", total, "\n")
  
  # Save as named matrix CSV (codes as row/col names)
  write.csv(tm, out_path)
  cat("  Saved:", out_path, "\n")
  
  return(tm)
}

################################################################################
# RUN locally
################################################################################
#cat("Starting disaggregated transition matrix pipeline...\n")
#cat("Years: ", paste(TARGET_YEARS, collapse = ", "), "\n\n")

#tasks <- counties_sf %>%
#  sf::st_drop_geometry() %>%
#  filter(STATEFP %in% TARGET_STATEFPS) %>%
#  select(GEOID, STATEFP)

#for (i in seq_len(nrow(tasks))) {
#  geoid   <- tasks$GEOID[i]
#  statefp <- tasks$STATEFP[i]
#  for (j in seq_len(length(TARGET_YEARS) - 1)) {
#    tmd <- compute_transition_disagg(TARGET_YEARS[j], TARGET_YEARS[j+1], geoid, statefp)
#    collapse_to_aggregated(tmd, TARGET_YEARS[j], TARGET_YEARS[j+1], geoid, statefp)
#  }
#}

#cat("ALL DONE: Disaggregated + aggregated transition matrices saved.\n")

################################################################################
# RUN (ANUBIS)
################################################################################
cat("Starting disaggregated transition matrix pipeline...\n")
cat("Years: ", paste(TARGET_YEARS, collapse = ", "), "\n\n")

library(parallel)

tasks <- counties_sf %>%
   sf::st_drop_geometry() %>%
   filter(STATEFP %in% TARGET_STATEFPS) %>%
   select(GEOID, STATEFP)

source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c("counties_sf", "TARGET_YEARS",
                     "compute_transition_disagg", "collapse_to_aggregated",
                     "DISAGG_PATH", "TRANSITION_DISAGG_PATH",
                     "TRANSITION_AGG_PATH", "code_to_category",
                     "category_labels", "tasks"))

parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
   library(terra)
   library(dplyr)
   geoid   <- tasks$GEOID[i]
   statefp <- tasks$STATEFP[i]
   for (j in seq_len(length(TARGET_YEARS) - 1)) {
     tmd <- compute_transition_disagg(TARGET_YEARS[j], TARGET_YEARS[j+1], geoid, statefp)
     collapse_to_aggregated(tmd, TARGET_YEARS[j], TARGET_YEARS[j+1], geoid, statefp)
   }
})

stopCluster(cl)
cat("ALL DONE: Disaggregated + aggregated transition matrices saved.\n")