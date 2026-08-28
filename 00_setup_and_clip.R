##############################################################################
# 00_setup_and_clip.R
# Step 1: Export county lookup table
# Step 2: Create directory structure
# Step 3: Verify CDL files exist + Projection is correct.
# Step 4: Clip national CDL to each county
#
# Run this script first before anything else.
# Later on TSE server: change TARGET_YEARS only.
##############################################################################
setwd("/users/rperilhou/extra_years/")
rm(list = ls())

library(terra)
library(dplyr)

# 1/ LOAD OFFICIAL STATE BOUNDARIES FROM CENSUS

cat("Loading state boundaries (from Census TIGER)\n")
states_sf   <- readRDS("data/SF/states_2016.rds")
counties_sf <- readRDS("data/SF/counties_2016.rds")

# All 48 contiguous states (STATEFP codes)
contiguous_statefps <- states_sf %>%
  sf::st_drop_geometry() %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  pull(STATEFP)

cat("Contiguous states:", length(contiguous_statefps), "\n")

TARGET_STATEFPS <- contiguous_statefps
TARGET_YEARS    <- c(2009:2020)

cat("Running for years:", paste(TARGET_YEARS, collapse = ", "), "\n")

# 2/ EXPORT COUNTY LOOKUP TABLE
# One-time reference for GEOID -> NAME mapping (for human inspection only;
# the pipeline uses GEOID and STATEFP exclusively).

counties_sf %>%
  sf::st_drop_geometry() %>%
  filter(STATEFP %in% TARGET_STATEFPS) %>%
  select(GEOID, STATEFP, COUNTYFP, NAME, LSAD) %>%
  write.csv("data/county_lookup.csv", row.names = FALSE)

cat("County lookup table saved to data/county_lookup.csv\n")

# 3/ CREATE DIRECTORY STRUCTURE

cat("Creating directory structure\n")

created <- 0
for (statefp in TARGET_STATEFPS) {
  dirs <- c(
    file.path("data/clipped",        statefp),
    file.path("outputs/classified",  statefp),
    file.path("outputs/transitions", statefp)
  )
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      created <- created + 1
    }
  }
}

cat("Done:", created, "directories created.\n")
# Does not overwrite so if there is any modification on the directory structure,
# need to delete the folders and re-run.

# 4/ VERIFY CDL FILES EXIST
# Expected: data/<year>_30m_cdls/<year>_30m_cdls.tif

cat("Checking CDL files:\n")

get_cdl_path <- function(year) {
  paste0("data/", year, "_30m_cdls/", year, "_30m_cdls.tif")
}

missing_files <- c()
for (year in TARGET_YEARS) {
  path <- get_cdl_path(year)
  if (file.exists(path)) {
    size_mb <- round(file.info(path)$size / 1e6, 1)
    cat(" CORRECT:", path, "(", size_mb, "MB)\n")
  } else {
    cat(" MISSING:", path, "\n")
    missing_files <- c(missing_files, path)
  }
}

if (length(missing_files) > 0) {
  stop("Missing CDL files. Please download the following before continuing:",
       paste(missing_files, collapse = "\n"))
} else {
  cat("All CDL files found.\n")
}

# 5/ VERIFY the projection used matches the official CDL CRS
# Can be found here: https://www.nass.usda.gov/Research_and_Science/Cropland/sarsfaqs2.php#common.2

rasters <- lapply(TARGET_YEARS, function(y) rast(get_cdl_path(y)))

cat("CRS:       ", crs(rasters[[1]], describe = TRUE)$name, "\n")
cat("Resolution:", res(rasters[[1]]), "metres\n\n")

for (i in seq_along(rasters)[-1]) {
  cat("vs", TARGET_YEARS[i],
      "— CRS:", same.crs(rasters[[1]], rasters[[i]]),
      "| res:", all(res(rasters[[1]]) == res(rasters[[i]])),
      "| ext:", ext(rasters[[1]]) == ext(rasters[[i]]), "\n")
}

# 6/ CLIP FUNCTION

clipped_path <- function(year, geoid, statefp) {
  file.path("data/clipped", statefp, paste0("CDL_", year, "_", geoid, ".tif"))
}

clip_county <- function(year, geoid, statefp, counties_sf) {

  out_path <- clipped_path(year, geoid, statefp)

  if (file.exists(out_path)) {
    cat("  Skipping (exists):", out_path, "\n")
    return(out_path)
  }

  cdl         <- rast(get_cdl_path(year))
  county_vect <- counties_sf %>% filter(GEOID == geoid) %>% vect()
  county_proj <- project(county_vect, crs(cdl))
  clipped     <- crop(cdl, county_proj) %>% mask(county_proj, touches = FALSE)

  writeRaster(clipped, out_path, overwrite = TRUE, datatype = "INT1U")
  cat("  Saved:", out_path, "\n")
  return(out_path)
}


# 7/ RUN CLIPPING (locally)

#tasks <- counties_sf %>%
#  sf::st_drop_geometry() %>%
#  filter(STATEFP %in% TARGET_STATEFPS) %>%
#  select(GEOID, STATEFP)

#for (i in seq_len(nrow(tasks))) {
#  geoid   <- tasks$GEOID[i]
#  statefp <- tasks$STATEFP[i]
#  for (year in TARGET_YEARS) {
#    clip_county(year, geoid, statefp, counties_sf)
#  }
#}

#cat("ALL DONE: Clipped files saved in data/clipped/\n")


# 7/ RUN CLIPPING (cluster)
library(parallel)

tasks <- counties_sf %>%
  sf::st_drop_geometry() %>%
  filter(STATEFP %in% TARGET_STATEFPS) %>%
  select(GEOID, STATEFP)

# ANUBIS cluster setup
source("/softs/R/createCluster.R")
cl <- createCluster()

# Export necessary objects to all workers
clusterExport(cl, c("counties_sf", "TARGET_YEARS",
                    "clip_county", "clipped_path", "get_cdl_path", "tasks"))

# Run clipping in parallel across all counties
parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  library(terra)
  library(dplyr)
  geoid   <- tasks$GEOID[i]
  statefp <- tasks$STATEFP[i]
  for (year in TARGET_YEARS) {
    clip_county(year, geoid, statefp, counties_sf)
  }
})

stopCluster(cl)
cat("ALL DONE: Clipped files saved in data/clipped/\n")

# At this stage each folder data/clipped/<STATEFP>/ should contain a .tif per county per year,
# named like: CDL_<year>_<GEOID>.tif.
