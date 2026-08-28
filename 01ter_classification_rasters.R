#01ter_classification_rasters.R

#   inputs:  outputs/<VERSION>/classified_disagg/<STATEFP>/Disagg_<YYYY>_<GEOID>.tif
#            data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")

#   outputs: outputs/<VERSION>/classified/<STATEFP>/Classified_<YYYY>_<GEOID>.tif

################################################################################
# CONFIGURATION
################################################################################
rm(list=ls())
library(tidyverse)
library(terra)
library(readxl)

setwd("/users/rperilhou/extra_years")

VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"
STATES_SF_PATH <- "data/SF/states_2016.rds"
COUNTIES_SF_PATH <- "data/SF/counties_2016.rds"

DISAGG_DIR <- file.path("outputs", VERSION, "classified_disagg")
DISAGG_PATH <- function(year, geoid, statefp){
  file.path(
    DISAGG_DIR, statefp, paste0("Disagg_", year, "_", geoid, ".tif")
  )
}


CLASSIFIED_DIR <- file.path("outputs", VERSION, "classified")
CLASSIFIED_PATH <- function(year, geoid, statefp){
  file.path(
    CLASSIFIED_DIR, statefp, paste0("Classified_", year, "_", geoid,".tif")
  )
}
################################################################################
#select counties

states_sf   <- readRDS(STATES_SF_PATH)
counties_sf <- readRDS(COUNTIES_SF_PATH)

contiguous_statefps <- states_sf %>%
  sf::st_drop_geometry() %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  pull(STATEFP)

TARGET_STATEFPS <- contiguous_statefps
TARGET_YEARS    <- c(2009:2020)

for (statefp in TARGET_STATEFPS) {
  d <- file.path(CLASSIFIED_DIR, statefp)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

#map crop types to categories

cdl_map <- read_excel(METADATA_XLSX, sheet = VERSION) %>%
  filter(Keep == 1)

reclass_table <- as.matrix(cdl_map[, c("CDL Code", "Category Code")])
dimnames(reclass_table) <- NULL

category_labels <- c("0"  = "NonCrop",
                     "1"  = "GM",
                     "2"  = "Tolerant",
                     "3"  = "Vulnerable",
                     "99" = "Unclassified")

#########
# Apply reclass_table to the disaggregated raster to produce the
# aggregated {0,1,2,3,99} category raster.
#########

classify_from_disagg <- function(year, geoid, statefp) {
  
  out_path <- CLASSIFIED_PATH(year, geoid, statefp)
  
  if (file.exists(out_path)) {
    cat("  Already exists, skipping:", out_path, "\n")
    return(out_path)
  }
  
  in_path <- DISAGG_PATH(year, geoid, statefp)
  
  if (!file.exists(in_path)) {
    stop("Missing disaggregated raster: ", in_path,
         "\nRun 01_mask_and_keep_disagg.R first.")
  }
  
  disagg <- rast(in_path)
  
  # Classify raw CDL codes into categories. Unmapped codes -> NA first.
  classified <- classify(disagg, reclass_table, others = NA)
  
  # Pixels that were in-mask (not NA in disagg) but unmapped -> Unclassified (99).
  # Pixels that were NA in disagg (outside union mask) stay NA.
  classified <- ifel(
    !is.na(disagg) & is.na(classified), 99L,
    classified
  )
  
  classified <- as.int(classified)
  
  writeRaster(classified, out_path, overwrite = TRUE, datatype = "INT1U")
  cat("  Saved:", out_path, "\n")
  
  return(out_path)
}



#########
# Run (locally)
#########

#tasks <- counties_sf %>%
#  sf::st_drop_geometry() %>%
#  filter(STATEFP %in% TARGET_STATEFPS) %>%
#  select(GEOID, STATEFP)

#for (i in seq_len(nrow(tasks))) {
#  geoid   <- tasks$GEOID[i]
#  statefp <- tasks$STATEFP[i]
#  cat("County:", geoid, "| State:", statefp, "\n")
#  for (year in TARGET_YEARS) {
#    classify_from_disagg(year, geoid, statefp)
#  }
#}

#cat("ALL DONE: Classified rasters saved in outputs/", VERSION, "classified/\n")

#########
# Run (ANUBIS)
#########

cat("Starting raster reclassification pipeline...\n")
cat("Years: ", paste(TARGET_YEARS, collapse = ", "), "\n")

library(parallel)

tasks <- counties_sf %>%
  sf::st_drop_geometry() %>%
  filter(STATEFP %in% TARGET_STATEFPS) %>%
  select(GEOID, STATEFP)

source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c("tasks", "TARGET_YEARS",
                    "classify_from_disagg", "CLASSIFIED_PATH", "DISAGG_PATH",
                    "reclass_table", "category_labels"))

parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  library(terra)
  geoid   <- tasks$GEOID[i]
  statefp <- tasks$STATEFP[i]
  for (year in TARGET_YEARS) {
    classify_from_disagg(year, geoid, statefp)
  }
})

stopCluster(cl)
cat("ALL DONE: Classified rasters saved in outputs/", VERSION, "classified/")


