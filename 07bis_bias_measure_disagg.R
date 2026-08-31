################################################################################
# 07bis_bias_measure_disagg.R
# Bias measures at the county x commodity level (disaggregated), analogous
# to 07_Bias_measure.R's county x category level.
#
#   Inputs:  outputs/<VERSION>/classification_summary_disagg.csv
#            outputs/<VERSION>/commodity_census_acres_2012.csv
#            outputs/<VERSION>/category_census_acres_2012.csv
#            data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")
#            data/metadata/<VERSION>/missing_geoid_census.csv

#   Outputs: outputs/<VERSION>/baseline_bias_disagg_commodity.csv
################################################################################
rm(list = ls())

setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(readxl)
VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"

SUMMARY_DISAG_PATH <- file.path(
  "outputs", VERSION, "classification_summary_disagg.csv")

COMMODITY_CENSUS_PATH <- file.path("outputs", VERSION, "commodity_census_acres_2012.csv")
CENSUS_MISSING_PATH <- file.path("data/metadata", VERSION, "missing_geoid_census.csv")
CENSUS_CAT_PATH <- file.path("outputs", VERSION, "category_census_acres_2012.csv")

BASELINE_BIAS_DISAGG_PATH <- file.path(
  "outputs", VERSION, "baseline_bias_disagg_commodity.csv")

################################################################################

# CDL code -> Census commodity mapping
cdl_census_map <- read_excel(METADATA_XLSX, sheet = VERSION) %>%
  filter(`Has Census` == 1) %>%
  mutate(`Census Commodity` = trimws(`Census Commodity`))

code_to_commodity <- cdl_census_map %>%
  select(`CDL Code`, `Census Commodity`, Category) %>%
  distinct()

cdl_commodity <- read_csv(SUMMARY_DISAG_PATH) %>%
  filter(year == 2012) %>%
  mutate(GEOID = sprintf("%05d", as.integer(geoid))) %>%
  left_join(code_to_commodity, by = c("cdl_code" = "CDL Code")) %>%
  filter(!is.na(`Census Commodity`)) %>%
  group_by(GEOID, commodity = `Census Commodity`, Category) %>%
  summarise(pixels_cdl = sum(n_pixels), .groups = "drop") %>%
  mutate(acres_cdl = pixels_cdl * 0.222395)

commodity_census <- read_csv(COMMODITY_CENSUS_PATH) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  select(GEOID, commodity, acres_census = acres, n_suppressed, n_total, fully_suppressed)

missing_geoid <- read_csv(CENSUS_MISSING_PATH) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  pull(GEOID) %>% unique()

fulldf <- cdl_commodity %>%
  filter(!GEOID %in% missing_geoid) %>%
  full_join(commodity_census, by = c("GEOID", "commodity")) %>%
  mutate(
    acres_cdl    = replace_na(acres_cdl, 0),
    pixels_cdl   = replace_na(pixels_cdl, 0),
    acres_census = replace_na(acres_census, 0),
    n_suppressed = replace_na(n_suppressed, 0),
    n_total      = replace_na(n_total, 0),
    fully_suppressed = replace_na(fully_suppressed, FALSE)
  )

commodity_to_category <- code_to_commodity %>%
  select(commodity = `Census Commodity`, Category) %>%
  distinct()

fulldf <- fulldf %>%
  select(-Category) %>%
  left_join(commodity_to_category, by = "commodity")

category_census_totals <- read_csv(CENSUS_CAT_PATH) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  group_by(GEOID) %>%
  summarise(acres_census_total = sum(acres_census, na.rm = TRUE), .groups = "drop")

bias_disagg <- fulldf %>%
  left_join(category_census_totals, by = "GEOID") %>%
  group_by(GEOID, Category) %>%
  mutate(acres_cdl_category_total = sum(acres_cdl, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    delta_acres = acres_cdl - acres_census,
    delta_pct_census = case_when(
      acres_census == 0 & acres_cdl == 0 ~ 0,
      acres_census == 0 & acres_cdl != 0 ~ NA_real_,
      TRUE ~ (acres_cdl - acres_census) / acres_census
    ),
    delta_pct_CDL = case_when(
      acres_census == 0 & acres_cdl == 0 ~ 0,
      acres_census != 0 & acres_cdl == 0 ~ NA_real_,
      TRUE ~ (acres_cdl - acres_census) / acres_cdl
    ),
    delta_pct_total = ifelse(acres_census_total > 0,
                             (acres_cdl - acres_census) / acres_census_total,
                             NA_real_),
    cdl_share_within_category = ifelse(acres_cdl_category_total > 0,
                                       acres_cdl / acres_cdl_category_total,
                                       NA_real_),
    accuracy_baseline = ifelse(
      acres_cdl == 0 & acres_census == 0,
      1,
      pmin(acres_cdl, acres_census) / pmax(acres_cdl, acres_census)
    )
  )

write_csv(bias_disagg, BASELINE_BIAS_DISAGG_PATH)

cat("Commodity-level rows:", nrow(bias_disagg), "| Counties:", n_distinct(bias_disagg$GEOID),
    "| Commodities:", n_distinct(bias_disagg$commodity), "\n")

