########
# I) Clip confidence raster to county boundaries and stack with disaggregated raster.
# Confidence is clipped from the national CDL confidence raster using the
# same county polygon + touches=FALSE as in 00_setup_and_clip.R,
# guaranteeing pixel-perfect alignment with the disaggregated raster.
#
#   Inputs:  outputs/<VERSION>/classified_disagg/<STATEFP>/Disagg_<YYYY>_<GEOID>.tif
#            data/confidence_NAT/CDL_conf_2012_national_aligned.tif

#   Outputs: outputs/<VERSION>/class_and_conf_disagg/<STATEFP>/Conf_stacked_disagg_<YYYY>_<GEOID>.tif
#########
rm(list = ls())
setwd("/users/rperilhou/extra_years")

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

CONF_OUTPUT_DIR <- file.path("outputs", VERSION, "class_and_conf_disagg")

DISAGG_DIR <- file.path("outputs", VERSION, "classified_disagg")
DISAGG_PATH <- function(year, geoid, statefp) {
  file.path(DISAGG_DIR, statefp, paste0("Disagg_", year, "_", geoid, ".tif"))
}

CONF_OUTPUT_PATH <- function(year, geoid, statefp) {
  file.path(CONF_OUTPUT_DIR, statefp, paste0("Conf_stacked_disagg_", year, "_", geoid, ".tif"))
}

################################################################################
# Load data
TARGET_YEAR     <- 2012

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
  
  # skip if disaggregated raster is missing
  disagg_file <- DISAGG_PATH(year, geoid, statefp)
  if (!file.exists(disagg_file)) next
  
  # create output directory if it doesn't exist
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # load disaggregated raster and clip confidence raster to county boundary
  disagg      <- rast(disagg_file)
  county_proj <- counties_proj[counties_proj$GEOID == geoid, ]
  confidence  <- crop(conf_nat, county_proj) %>% mask(county_proj, touches = FALSE)
  
  if (!compareGeom(disagg, confidence, stopOnError = FALSE)) {
    warning("Extent mismatch GEOID ", geoid, " -- skipping.")
    next
  }
  
  stacked        <- c(disagg, confidence)
  names(stacked) <- c("cdl_code", "confidence")
  writeRaster(stacked, out_path, overwrite = TRUE, datatype = "INT1U")
  cat("Saved:", out_path, "\n")
}

cat("Done.\n")