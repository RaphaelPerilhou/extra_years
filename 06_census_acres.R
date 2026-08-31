################################################################################
# 06_Census_Acres.R
# Aggregates raw Census area records to commodity- and category-level acres,
# distinguishing genuine disclosure suppression (D) from near-zero
# rounding (Z).
#
#   Inputs:  data/CENSUS/<VERSION>/census_area_2012.csv
#            data/county_lookup.csv
#            data/metadata/<VERSION>/census_shortdesc_map.csv
#            data/metadata/CDL_CENSUS_MAP_meta.xlsx (sheet "<VERSION>")
#
#   Outputs: outputs/<VERSION>/commodity_census_acres_2012.csv
#            outputs/<VERSION>/category_census_acres_2012.csv
#            outputs/<VERSION>/census_suppressed_log.csv
################################################################################
rm(list = ls())

setwd("/users/rperilhou/extra_years")

library(tidyverse)
library(readxl)
VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
COUNTY_LOOKUP_PATH <- "data/county_lookup.csv"
METADATA_XLSX <- "data/metadata/CDL_CENSUS_MAP_meta.xlsx"

CENSUS_AREA_PATH <- file.path("data/CENSUS", VERSION, "census_area_2012.csv")
SHORTDESC_MAP_PATH <- file.path("data/metadata", VERSION, "census_shortdesc_map.csv")

AGG_COMMO_PATH <- file.path("outputs", VERSION, "commodity_census_acres_2012.csv")
CENSUS_CAT_PATH <- file.path("outputs", VERSION, "category_census_acres_2012.csv")
SUPPRESSED_PATH <- file.path("outputs", VERSION, "census_suppressed_log.csv")


# Load data
census_area <- read.csv(CENSUS_AREA_PATH) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)))

county_lookup <- read.csv(COUNTY_LOOKUP_PATH, stringsAsFactors = FALSE) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)))

shortdesc_map <- read.csv(SHORTDESC_MAP_PATH)
cdl_census_map <- read_excel(METADATA_XLSX, sheet = VERSION)

length(setdiff(county_lookup$GEOID, census_area$GEOID))

# Step 1: Flag and clean suppressed values
# (D) = withheld to avoid disclosing data for individual operations
# (Z) = less than half the unit shown
# Others can be found here: https://quickstats.nass.usda.gov/src/glossary.pdf

census_area <- census_area %>%
  mutate(
    value = gsub(",", "", value),
    is_Z = grepl("\\(Z\\)", value),          # treat as 0, not suppressed
    is_D = grepl("\\(D\\)", value),          # true suppression
    suppressed = is_D,                        # only (D) counts as suppressed going forward
    value_clean = case_when(
      is_Z ~ 0,                               # (Z) becomes a real, known 0
      is_D ~ 0,                               # (D) stays zero-filled, but flagged as suppressed
      grepl("^[0-9]+$", value) ~ as.numeric(value),  # genuine numeric values
      TRUE ~ NA_real_                         # anything else unexpected NA, worth investigating
    )
  )

sum(is.na(census_area$value_clean))

suppressed_log <- census_area %>%
  filter(suppressed) %>%
  select(GEOID, commodity, short_desc, value)

cat("Suppressed non_numeric value:", nrow(suppressed_log), "\n")

# Step 2: Aggregate subcategories of commodity (e.g bell pepper + chile pepper)
agg_commodity <- census_area %>%
  group_by(GEOID, commodity) %>%
  summarise(
    acres = sum(value_clean, na.rm = TRUE),
    n_suppressed = sum(suppressed),
    n_total = n(),
    fully_suppressed = n_suppressed == n_total,
    .groups = "drop"
  )

# Step 3: Join category from CDL map
commodity_category <- cdl_census_map %>%
  filter(`Has Census` == 1) %>%
  mutate(`Census Commodity` = trimws(`Census Commodity`)) %>%
  select(`Census Commodity`, Category, `Category Code`) %>%
  distinct()

commodity_category <- agg_commodity %>%
  left_join(commodity_category, by = c("commodity" = "Census Commodity"))

# Step 4: Aggregate to category level
census_category <- commodity_category %>%
  group_by(GEOID, `Category Code`, Category) %>%
  summarise(
    acres_census = sum(acres, na.rm = TRUE),
    n_suppressed = sum(n_suppressed),
    .groups = "drop"
  ) %>%
  arrange(`Category Code`)

length(unique(census_category$GEOID)) #Not 3 rows per county because some county don't have the 3 categories

# Output
# wanna check for some counties:
#if want every GEOID:
# county_list <- unique(census_category$GEOID)
county_list <- c("01033")

cat("\n── Commodity-level acres ──\n")
for (i in county_list) {
  name <- county_lookup$NAME[county_lookup$GEOID == i][1]
  cat("\n", i, "—", name, "\n")
  agg_commodity %>%
    filter(GEOID == i) %>%
    select(commodity, acres, n_suppressed) %>%
    print()
}

cat("\n── Category-level acres ──\n")
for (i in county_list) {
  name <- county_lookup$NAME[county_lookup$GEOID == i][1]
  cat("\n", i, "—", name, "\n")
  census_category %>%
    filter(GEOID == i) %>%
    select(Category, acres_census, n_suppressed) %>%
    print()
}

cat("\n── Suppressed short_desc ──\n")
for (i in county_list) {
  name <- county_lookup$NAME[county_lookup$GEOID == i][1]
  cat("\n", i, "-", name, "\n")
  suppressed_log %>%
    filter(GEOID == i) %>%
    select(commodity, value) %>%
    print()
}

#Save outputs (for all GEOID)

write_csv(
  agg_commodity %>% arrange(GEOID, commodity),
  AGG_COMMO_PATH
)

write_csv(
  census_category %>% arrange(GEOID, `Category Code`),
  CENSUS_CAT_PATH
)

write_csv(
  suppressed_log %>% arrange(GEOID, commodity),
  SUPPRESSED_PATH
)

#
n_total <- nrow(census_area)
n_suppressed <- sum(census_area$suppressed)
pct_suppressed <- round(100 * n_suppressed / n_total, 2)

cat("Total observations:", n_total, "\n")
cat("Suppressed observations:", n_suppressed, "\n")
cat("Percent suppressed:", pct_suppressed, "%\n")

n_commodity_obs <- nrow(agg_commodity)
n_fully_suppressed <- sum(agg_commodity$fully_suppressed)
n_partially_suppressed <- sum(agg_commodity$n_suppressed > 0 & !agg_commodity$fully_suppressed)

pct_fully_suppressed <- round(100 * n_fully_suppressed / n_commodity_obs, 2)
pct_partially_suppressed <- round(100 * n_partially_suppressed / n_commodity_obs, 2)

cat("Total county-commodity pairs:", n_commodity_obs, "\n")
cat("Fully suppressed (acres = 0, true unknown):", n_fully_suppressed, "(", pct_fully_suppressed, "%)\n")
cat("Partially suppressed (real acreage retained):", n_partially_suppressed, "(", pct_partially_suppressed, "%)\n")

