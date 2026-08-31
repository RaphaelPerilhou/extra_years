################################################################################
# 05bis_census_area.R
#
# Builds census_shortdesc_map and census_area_2012 directly from the
# CDL_CENSUS_MAP_meta.xlsx "Census Short_desc" column.
#
#   Inputs:  data/qs.census2012.txt
#            data/county_lookup.csv
#            data/SF/states_2016.rds
#            data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")
#
#   Outputs: data/metadata/<VERSION>/census_shortdesc_map.csv
#            data/CENSUS/<VERSION>/census_area_2012.csv
#            data/metadata/<VERSION>/missing_geoid_census.csv
################################################################################
rm(list = ls())

setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(readxl)
library(data.table)

VERSION <- "v2"   # or "v1"

################################################################################
# PATHS BUILDER
################################################################################
CENSUS_PATH        <- "data/qs.census2012.txt"                                
COUNTY_LOOKUP_PATH <- "data/county_lookup.csv"                                 
STATES_SF_PATH      <- "data/SF/states_2016.rds"                               
METADATA_XLSX       <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"                

V_META_DIR   <- file.path("data/metadata", VERSION)
V_CENSUS_DIR <- file.path("data/CENSUS", VERSION)

SHORTDESC_MAP_PATH  <- file.path(V_META_DIR, "census_shortdesc_map.csv")
CENSUS_AREA_PATH    <- file.path(V_CENSUS_DIR, "census_area_2012.csv")
CENSUS_MISSING_PATH <- file.path(V_META_DIR, "missing_geoid_census.csv")

dir.create(V_META_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(V_CENSUS_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1/ LOAD CENSUS + COUNTY REFERENCE (unchanged from 05)
################################################################################
census <- fread(CENSUS_PATH, sep = "\t", quote = "")

availability <- census[
  SOURCE_DESC == "CENSUS" &
    AGG_LEVEL_DESC == "COUNTY" &
    YEAR == 2012 &
    DOMAIN_DESC == "TOTAL"
][, .(
  GEOID = sprintf("%02d%03d", STATE_FIPS_CODE, COUNTY_CODE),
  county_name = COUNTY_NAME,
  state_fips = STATE_FIPS_CODE,
  commodity = COMMODITY_DESC,
  statisticcat = STATISTICCAT_DESC,
  short_desc = SHORT_DESC,
  value = VALUE
)]

county_lookup <- read.csv(COUNTY_LOOKUP_PATH, stringsAsFactors = FALSE)
county_lookup$GEOID_padded <- sprintf("%05d", as.integer(county_lookup$GEOID))

states_sf <- readRDS(STATES_SF_PATH)
contiguous_statefps <- states_sf %>%
  sf::st_drop_geometry() %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  pull(STATEFP)

availability_clean <- availability[GEOID %in% county_lookup$GEOID_padded]

keep_stats <- c(
  "AREA HARVESTED", "AREA BEARING", "AREA NON-BEARING",
  "AREA BEARING & NON-BEARING", "AREA IN PRODUCTION",
  "AREA GROWN", "AREA NOT HARVESTED", "AREA"
)
availability_clean <- availability_clean[statisticcat %in% keep_stats]

missing_geoids <- setdiff(county_lookup$GEOID_padded, unique(availability$GEOID))

################################################################################
# 2/ BUILD census_shortdesc_map DIRECTLY FROM THE EXCEL (VERSION-dependent)
################################################################################
cdl_census_map <- read_excel(METADATA_XLSX, sheet = VERSION)

census_shortdesc_map <- cdl_census_map %>%
  filter(`Has Census` == 1) %>%
  mutate(`Census Commodity` = trimws(`Census Commodity`)) %>%
  select(commodity = `Census Commodity`, short_desc = `Census Short_desc`) %>%
  distinct() %>%
  separate_longer_delim(short_desc, delim = " + ") %>%
  mutate(short_desc = trimws(short_desc))

write.csv(census_shortdesc_map, SHORTDESC_MAP_PATH, row.names = FALSE)

################################################################################
# 3/ PULL THE SELECTED ROWS AND SAVE census_area_2012
################################################################################
area_stats <- availability_clean[
  statisticcat %in%
    c("AREA HARVESTED", "AREA BEARING & NON-BEARING", "AREA IN PRODUCTION",
      "AREA GROWN") &
    grepl(
      "ACRES HARVESTED$|ACRES BEARING & NON-BEARING$|ACRES IN PRODUCTION$|ACRES GROWN$",
      short_desc
    ) &
    !grepl("IRRIGATED", short_desc)
]

selected_descs <- census_shortdesc_map$short_desc

census_area <- area_stats[short_desc %in% selected_descs] %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)))

#Counties with no commodities of interest in the CENSUS
missing_county_shortdesc <- setdiff(
  unique(county_lookup$GEOID_padded),
  unique(census_area$GEOID)
) 

write.csv(census_area, CENSUS_AREA_PATH, row.names = FALSE)

################################################################################
# 4/ MISSING-GEOID BOOKKEEPING
################################################################################
census_missing <- county_lookup[
  county_lookup$GEOID_padded %in% missing_geoids,
] %>%
  mutate(reason = "Not in Census")

write.csv(census_missing, CENSUS_MISSING_PATH, row.names = FALSE)

cat("\nDone. VERSION =", VERSION, "\n")

