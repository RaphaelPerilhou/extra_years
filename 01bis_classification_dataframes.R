#01bis_classification_dataframes.R
#
#   input:  outputs/<VERSION>/classification_summary_disagg_parts/<GEOID>.csv
#           data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")

#   output: outputs/<VERSION>/classification_summary_disagg.csv (long format)
#           outputs/<VERSION>/classification_summary_disagg_wide.csv (wide format)
#           outputs/<VERSION>/classification_summary.csv (recovered aggregated counts)

# Then we use the wide dataframe to recover the aggregated dataframe using the
# crop_type -> category {0,1,2,3,99} mapping (reclass_table)

rm(list=ls())

setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(terra)
library(readxl)

################################################################################
# CONFIGURATION
################################################################################

VERSION <- "v2"   # or "v2" (=choose the sheet of CDL_CENSUS_MAP_meta.xlsx)

################################################################################
# PATHS: every directory/file this script reads or writes.
################################################################################
STATES_SF_PATH <- "data/SF/states_2016.rds"
COUNTIES_SF_PATH <- "data/SF/counties_2016.rds"
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"

SUMMARY_DISAG_PARTS_PATH <- file.path(
  "outputs",VERSION, "classification_summary_disagg_parts")

SUMMARY_DISAG_PATH <- file.path(
  "outputs", VERSION, "classification_summary_disagg.csv")

WIDE_SUMMARY_DISAG_PATH <- file.path(
  "outputs", VERSION, "classification_summary_disagg_wide.csv")

AGGREGATED_WIDE_PATH <- file.path(
  "outputs", VERSION, "classification_summary.csv")



states_sf <- readRDS(STATES_SF_PATH)
counties_sf <- readRDS(COUNTIES_SF_PATH)

#Select counties
# All 48 contiguous states (STATEFP codes)
contiguous_statefps <- states_sf %>%
  sf::st_drop_geometry() %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  pull(STATEFP)

cat("Contiguous states available:", length(contiguous_statefps), "\n")

TARGET_STATEFPS <- contiguous_statefps
TARGET_YEARS    <- c(2009:2020)

################################################################################
# DISAGGREGATED SUMMARY -> WIDE DATAFRAME
################################################################################

part_files <- list.files(
  SUMMARY_DISAG_PARTS_PATH,
  full.names = TRUE)

long <- read_csv(part_files, show_col_types = F)

write_csv(long, SUMMARY_DISAG_PATH)

wide <- long %>% 
  pivot_wider(
    id_cols = c(statefp, geoid, year),
    names_from = cdl_code,
    values_from = n_pixels,
    values_fill = 0
    )

write_csv(wide, WIDE_SUMMARY_DISAG_PATH)
cat("Wide-format disaggregated classification is complete:", 
    "one row per county-year pair.")


################################################################################
# DISAGGREGATED -> AGGREGATED
# Recover the aggregated category counts directly from the long table,
# by mapping cdl_code -> category and summing pixel counts.
# No raster reclassification needed since we already have per-code counts.
################################################################################

# Map crop types to categories
cdl_map <- read_excel(METADATA_XLSX, sheet = VERSION) %>%
  filter(Keep == 1)

code_to_category <- setNames(cdl_map$`Category Code`, as.character(cdl_map$`CDL Code`))

category_labels <- c("0"  = "NonCrop",
                     "1"  = "GM",
                     "2"  = "Tolerant",
                     "3"  = "Vulnerable",
                     "99" = "Unclassified")

aggregated <- long %>%
  mutate(
    category = code_to_category[as.character(cdl_code)],
    category = ifelse(is.na(category), 99, category)  # unmapped codes -> Unclassified
  ) %>%
  group_by(statefp, geoid, year, category) %>%
  summarise(n_pixels = sum(n_pixels), .groups = "drop")

aggregated_wide <- aggregated %>%
  mutate(category_label = category_labels[as.character(category)]) %>%
  select(-category) %>%
  pivot_wider(
    id_cols = c(statefp, geoid, year),
    names_from = category_label,
    values_from = n_pixels,
    values_fill = 0
  ) %>%
  mutate(total = NonCrop + GM + Tolerant + Vulnerable + Unclassified) %>%
  select(statefp, geoid, year, NonCrop, GM, Tolerant, Vulnerable, Unclassified, total)

write_csv(aggregated_wide, AGGREGATED_WIDE_PATH)
