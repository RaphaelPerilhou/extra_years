########
# I) Clip confidence raster to county boundaries and stack with classified raster.

# Confidence is clipped from the national CDL confidence raster using the
# same county polygon + touches=FALSE as in 00_setup_and_clip.R,
# guaranteeing pixel-perfect alignment with the classified raster.
#
#   Inputs:  outputs/<VERSION>/classified/<STATEFP>/Classified_<YYYY>_<GEOID>.tif
#            data/confidence_NAT/CDL_conf_2012_national_aligned.tif

#   Outputs: outputs/<VERSION>/class_and_conf/<STATEFP>/Conf_stacked_<YYYY>_<GEOID>.tif
#########
setwd("/users/rperilhou/extra_years")
rm(list = ls())

library(terra)
library(dplyr)
library(sf)

VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
COUNTIES_SF_PATH <- "data/SF/counties_2016.rds"
COUNTY_LOOKUP_PATH <- "data/county_lookup.csv"

CONF_NAT_PATH   <- "data/confidence_NAT/CDL_conf_2012_national_aligned.tif"

CONF_OUTPUT_DIR <- file.path("outputs", VERSION, "class_and_conf")
CONF_OUTPUT_PATH <- function(year, geoid, statefp) {
  file.path(CONF_OUTPUT_DIR, statefp, paste0("Conf_stacked_", year, "_", geoid, ".tif"))
}

CLASSIFIED_DIR <- file.path("outputs", VERSION, "classified")
CLASSIFIED_PATH <- function(year, geoid, statefp) {
  file.path(CLASSIFIED_DIR, statefp, paste0("Classified_", year, "_", geoid, ".tif"))
}

################################################################################
# Parameters
TARGET_YEAR <- 2012
# Load data

counties_sf   <- readRDS(COUNTIES_SF_PATH)
county_lookup <- read.csv(COUNTY_LOOKUP_PATH, colClasses = "character")

cat("Loading national confidence raster...\n")
conf_nat <- rast(CONF_NAT_PATH)
cat("  Extent:", as.vector(ext(conf_nat)), "\n")
cat("  CRS:   ", crs(conf_nat, describe = TRUE)$name, "\n")

# Pre-project all counties once, instead of re-projecting per county
cat("Pre-projecting county boundaries...\n")
counties_proj <- project(vect(counties_sf), crs(conf_nat))

# Loop over each county
tasks <- county_lookup %>% select(GEOID, STATEFP)
cat("Counties:", nrow(tasks), "\n")

for (i in seq_len(nrow(tasks))) {
  
  geoid   <- tasks$GEOID[i]
  statefp <- tasks$STATEFP[i]
  year    <- TARGET_YEAR
  
  # skip if output already exists
  out_path <- CONF_OUTPUT_PATH(year, geoid, statefp)
  if (file.exists(out_path)) next
  
  # skip if classified raster is missing
  classified_file <- CLASSIFIED_PATH(year, geoid, statefp)
  if (!file.exists(classified_file)) next
  
  # create output directory if it doesn't exist
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # load classified raster and clip confidence raster to county boundary
  classified  <- rast(classified_file)
  county_proj <- counties_proj[counties_proj$GEOID == geoid, ]
  confidence  <- crop(conf_nat, county_proj) %>% mask(county_proj, touches = FALSE)
  
  if (!compareGeom(classified, confidence, stopOnError = FALSE)) {
    warning("Extent mismatch GEOID ", geoid, " -- skipping.")
    next
  }
  
  stacked        <- c(classified, confidence)
  names(stacked) <- c("category", "confidence")
  writeRaster(stacked, out_path, overwrite = TRUE, datatype = "INT1U")
  cat("Saved:", out_path, "\n")
}

cat("Done.\n")