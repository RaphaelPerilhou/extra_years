################################################################################
# 01_mask_and_keep_disagg.R
#
# Disaggregated variant of 01_mask_and_classify.R 
#
# Step 1: Build agricultural mask per year
# Step 2: Build union mask across all study years 
# Step 3: Apply union mask to CDL layers
# Step 4: KEEP RAW CDL CODES (no reclassification).
#         - Pixels inside union mask keep their CDL code that year
#           (including non-agricultural codes, e.g. 63 Forest, which
#           correspond to 99 in the aggregated pipeline)
#         - Pixels outside union mask stay NA (excluded entirely)
#
# Inputs:  data/clipped/<statefp>/CDL_<year>_<geoid>.tif
# Also need: data/SF/states_2016.rds
#            data/SF/counties_2016.rds
#            data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")
#
# Outputs: outputs/<VERSION>/classified_disagg/<statefp>/Disagg_<year>_<geoid>.tif
#          outputs/<VERSION>/classification_summary_disagg_parts/<geoid>.csv
#
# The union mask is built from the SAME ag_codes as the aggregated pipeline,
# so the pixel universe is identical and the aggregated categories can be
# recovered exactly by mapping codes through reclass_table (unmapped -> 99).
################################################################################
# when run locally: mem.maxVSize(vsize = 64000)  # In MB, so 64000 = 64GB
rm(list = ls())

setwd("/users/rperilhou/extra_years")

library(terra)
library(tidyverse)
library(readxl)
################################################################################
# CONFIGURATION
################################################################################

VERSION <- "v2"   # or "v2" (=choose the sheet of CDL_CENSUS_MAP_meta.xlsx)

################################################################################
# PATHS: every directory/file this script reads or writes.
################################################################################
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"

STATES_SF_PATH <- "data/SF/states_2016.rds"
COUNTIES_SF_PATH <- "data/SF/counties_2016.rds"

CLIPPED_DIR <- "data/clipped" 
DISAGG_DIR  <- file.path("outputs", VERSION, "classified_disagg")
SUMMARY_PARTS_DIR <- file.path("outputs", VERSION, "classification_summary_disagg_parts")

################################################################################
# PATH BUILDERS
################################################################################
clipped_path <- function(year, geoid, statefp) {
  file.path(CLIPPED_DIR, statefp, paste0("CDL_", year, "_", geoid, ".tif"))
}

disagg_path <- function(year, geoid, statefp) {
  file.path(DISAGG_DIR, statefp, paste0("Disagg_", year, "_", geoid, ".tif"))
}

summary_part_path <- function(geoid) {
  file.path(SUMMARY_PARTS_DIR, paste0(geoid, ".csv"))
}

################################################################################
#Select counties and years
states_sf <- readRDS(STATES_SF_PATH)
counties_sf <- readRDS(COUNTIES_SF_PATH)

# All 48 contiguous states (STATEFP codes)
contiguous_statefps <- states_sf %>%
  sf::st_drop_geometry() %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  pull(STATEFP)

cat("Contiguous states available:", length(contiguous_statefps), "\n")

TARGET_YEARS <- c(2009:2020)
TARGET_STATEFPS <- contiguous_statefps

# Create output directories
for (statefp in TARGET_STATEFPS) {
  d <- file.path(DISAGG_DIR, statefp)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
if (!dir.exists(SUMMARY_PARTS_DIR)) dir.create(SUMMARY_PARTS_DIR, recursive = TRUE)

################################################################################
# 1/ AGRICULTURAL MASK CODES
################################################################################
#keep only the  cdl codes we defined in the metadata as in the agricultural mask. 

ag_codes <- read_excel(METADATA_XLSX,
                       sheet = VERSION) %>%
  filter(Keep == 1) %>%
  pull(`CDL Code`) %>%
  sort()
################################################################################
# 2/ CORE FUNCTIONS
################################################################################
# 2.1) Build agricultural mask for one year 
make_mask <- function(year, geoid, statefp) {
  cdl <- rast(clipped_path(year, geoid, statefp))
  ifel(cdl %in% ag_codes, 1, NA)
}

# 2.2) Build union mask across all years for one county
# Pixel = 1 if EVER agricultural in any year
make_union_mask <- function(years, geoid, statefp) {
  cat("  Building union mask across", length(years), "years\n")
  masks <- lapply(years, function(y) make_mask(y, geoid, statefp))
  
  union <- masks[[1]]
  for (i in seq_along(masks)[-1]) {
    union <- ifel(!is.na(union) | !is.na(masks[[i]]), 1, NA)
  }
  
  n_ag_pixels <- sum(values(union) == 1, na.rm = TRUE)
  cat("  Union mask:", n_ag_pixels, "agricultural pixels retained\n")
  return(union)
}

# 2.3) Keep raw CDL codes for one year within the union mask
keep_year <- function(year, geoid, statefp, union_mask) {
  
  out_path <- disagg_path(year, geoid, statefp)
  
  if (file.exists(out_path)) {
    cat("  Already exists, skipping:", out_path, "\n")
    return(out_path)
  }
  
  # Load clipped CDL and apply union mask.
  # No reclassification: every in-mask pixel keeps its raw CDL code
  # Out-of-mask pixels stay NA.
  cdl    <- rast(clipped_path(year, geoid, statefp))
  masked <- mask(cdl, union_mask)
  masked <- as.int(masked)
  
  # Per-code pixel counts (long format)
  code_freq <- freq(masked) 
  
  total <- sum(code_freq$count)
  cat("  Distinct CDL codes in mask:", nrow(code_freq),
      "| Total pixels:", total,
      "(union mask has", sum(values(union_mask) == 1, na.rm = TRUE), ")\n")
  
  # Per-county summary file, avoids concurrent writes across parallel workers
  # We will bind each file after.
  summary_path <- summary_part_path(geoid)
  
  counts_rows <- data.frame(
    statefp  = statefp,
    geoid    = geoid,
    year     = year,
    cdl_code = code_freq$value,
    n_pixels = code_freq$count
  )
  
  if (file.exists(summary_path)) {
    write.table(counts_rows, summary_path, append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)
  } else {
    write.csv(counts_rows, summary_path, row.names = FALSE)
  }
  
  writeRaster(masked, out_path, overwrite = TRUE, datatype = "INT1U")
  cat("  Saved:", out_path, "\n")
  
  return(out_path)
}

# Full pipeline for one county
run_county <- function(geoid, statefp, years) {
  cat("####################################\n")
  cat("STATEFP:", statefp, "| GEOID:", geoid, "| Years:", min(years), "-", max(years), "\n")
  cat("####################################\n")
  
  # Check all clipped files exist
  missing <- years[!file.exists(sapply(years, clipped_path, geoid = geoid, statefp = statefp))]
  if (length(missing) > 0) {
    stop("Missing clipped files for GEOID ", geoid, " years: ",
         paste(missing, collapse = ", "),
         "\nRun 00_setup_and_clip.R first.")
  }
  
  # Skip entire county if all disaggregated files already exist
  # (avoids recomputing the union mask for completed counties on restart)
  all_done <- all(file.exists(sapply(years, disagg_path, geoid = geoid, statefp = statefp)))
  if (all_done) {
    cat("  All years already done, skipping county.\n")
    return(invisible(NULL))
  }
  
  # Build union mask once for all years
  union_mask <- make_union_mask(years, geoid, statefp)
  
  # Mask each year, keeping raw codes
  paths <- lapply(years, function(y) {
    cat("    - Year:", y, "\n")
    keep_year(y, geoid, statefp, union_mask)
  })
  
  cat("County", geoid, "complete.\n")
  return(paths)
}
################################################################################
# 3/ RUN locally
################################################################################
#cat("Years: ", paste(TARGET_YEARS, collapse = ", "), "\n")

#tasks <- counties_sf %>%
#  sf::st_drop_geometry() %>%
#  filter(STATEFP %in% TARGET_STATEFPS) %>%
#  select(GEOID, STATEFP)

#for (i in seq_len(nrow(tasks))) {
#  geoid   <- tasks$GEOID[i]
#  statefp <- tasks$STATEFP[i]
#  run_county(geoid, statefp, TARGET_YEARS)
#}

#cat("ALL DONE: Disaggregated files saved in outputs/<VERSION>/classified_disagg/\n")


################################################################################
# 3/ RUN (cluster)
################################################################################
cat("Years: ", paste(TARGET_YEARS, collapse = ", "), "\n")

library(parallel)

tasks <- counties_sf %>%
  sf::st_drop_geometry() %>%
  filter(STATEFP %in% TARGET_STATEFPS) %>%
  select(GEOID, STATEFP)

source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c("counties_sf", "TARGET_YEARS",
                    "run_county", "make_union_mask", "make_mask",
                    "keep_year", "clipped_path", "disagg_path", "summary_part_path",
                    "ag_codes", "tasks"))

parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  library(terra)
  library(dplyr)
  geoid   <- tasks$GEOID[i]
  statefp <- tasks$STATEFP[i]
  run_county(geoid, statefp, TARGET_YEARS)
})

stopCluster(cl)
cat("ALL DONE: Disaggregated files saved in outputs/", VERSION, "/classified_disagg/\n", sep = "")
