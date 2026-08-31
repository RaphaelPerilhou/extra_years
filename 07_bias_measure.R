################################################################################
# 07_bias_measure.R
# County x category bias table: CDL acres vs Census acres.
#
#   Inputs: outputs/<VERSION>/category_census_acres_2012.csv
#           outputs/<VERSION>/classification_summary.csv
#           data/metadata/<VERSION>/missing_geoid_census.csv
#
#   Outputs: outputs/<VERSION>/baseline_bias.csv
#            data/metadata/<VERSION>/zero_census_acres_geoids.csv
################################################################################

rm(list = ls())
setwd("/users/rperilhou/extra_years")
library(tidyverse)
library(DescTools)
VERSION <- "v2"

################################################################################
# PATHS BUILDER
################################################################################
CENSUS_ACRES_PATH <- file.path(
  "outputs", VERSION, "category_census_acres_2012.csv")

CLASSIFICATION_PATH <- file.path(
  "outputs", VERSION, "classification_summary.csv")

CENSUS_MISSING_PATH <- file.path(
  "data/metadata", VERSION, "missing_geoid_census.csv")

ZERO_CENSUS_PATH <- file.path(
  "data/metadata", VERSION, "zero_census_acres_geoids.csv")

BASELINE_BIAS_PATH <- file.path(
  "outputs", VERSION, "baseline_bias.csv")

################################################################################


census_acres <- read_csv(CENSUS_ACRES_PATH) %>% 
  select(GEOID, Category, acres_census, n_suppressed)

sum(is.na(census_acres$acres_census))

classification_summary <- read_csv(CLASSIFICATION_PATH) %>%
  filter(year == 2012) %>%
  mutate(GEOID = sprintf("%05d", as.integer(geoid))) %>%
  select(GEOID, GM, Tolerant, Vulnerable) %>%
  pivot_longer(
    cols = c(GM, Tolerant, Vulnerable),
    names_to = "Category",
    values_to = "pixels_cdl"
  ) %>%
  mutate(
    acres_cdl = pixels_cdl * 0.222395)

## Check for missing rows (example; a county with only 2 categories has 2 rows instead of 3).
# We want to replace the missing row by 0 so it gets a value. 

classification_summary %>% 
  count(GEOID, name = "n_row_counties") %>% 
  count(n_row_counties, name = "n_counties") #No missing row in CDL

census_acres %>%
  count(GEOID, name = "n_row_counties") %>%
  count(n_row_counties, name = "n_counties") #Missing rows in Census, we replace by 0. 


#create the bias variable and keep only counties available in census.
missing_geoid <- read_csv(CENSUS_MISSING_PATH) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)))
missing_geoid <- unique(missing_geoid$GEOID)

fulldf <- classification_summary %>% 
  select(GEOID, Category, pixels_cdl, acres_cdl) %>% 
  filter(!GEOID %in% missing_geoid) %>% 
  left_join(census_acres, by = c("GEOID", "Category"))

fulldf %>% 
  count(GEOID, name = "n_row_counties") %>% 
  count(n_row_counties, name = "n_counties") #Now, we have 3 rows per county,
#but we have NAs. 

fulldf <- fulldf %>% 
  mutate(acres_census = replace_na(acres_census, 0),
         n_suppressed = replace_na(n_suppressed, 0)
  )

bias <- fulldf %>% 
  select(GEOID, Category, n_suppressed, pixels_cdl, acres_cdl, acres_census) %>%
  group_by(GEOID) %>% 
  mutate(
    # raw delta
    delta_acres = acres_cdl - acres_census,
    
    # delta pct census
    delta_pct_census = case_when(
      acres_census == 0 & acres_cdl == 0 ~ 0,
      acres_census == 0 & acres_cdl != 0 ~ NA_real_,
      TRUE ~ (acres_cdl - acres_census) / acres_census
    ),
    
    # delta pct cdl
    delta_pct_CDL = case_when(
      acres_census == 0 & acres_cdl == 0 ~ 0,
      acres_census != 0 & acres_cdl == 0 ~ NA_real_,
      TRUE ~ (acres_cdl - acres_census) / acres_cdl
    ),
    
    # delta pct total
    acres_census_total = sum(acres_census),
    delta_pct_total = (acres_cdl - acres_census) / acres_census_total,
    
    # accuracy and weighted accuracy
    acres_cdl_total = sum(acres_cdl),
    cdl_share = acres_cdl / acres_cdl_total,
    accuracy_baseline = ifelse(
      acres_cdl == 0 & acres_census == 0,
      1,
      pmin(acres_cdl, acres_census) / pmax(acres_cdl, acres_census)
    ),
    weighted_accuracy_baseline = accuracy_baseline * cdl_share
  ) %>% 
  ungroup()

#Some checks
sum(is.na(bias$acres_census)) #No 0 because we replace missing rows (no values reported) by 0.
sum(is.na(bias$acres_cdl))

sum(is.na(bias$delta_pct_census)) #We check if it's only when cdl != 0 and census == 0
sum(bias$acres_census == 0 & bias$acres_cdl != 0) #good

sum(is.na(bias$delta_pct_CDL)) #We check if it's only when cdl == 0 and census != 0
sum(bias$acres_census != 0 & bias$acres_cdl == 0) #good

sum(is.na(bias$accuracy_baseline))
mean(bias$accuracy_baseline)

sum(is.na(bias$weighted_accuracy_baseline))
sum(is.na(bias$cdl_share)) #The only possibility is that cdl_share is NA because 
# acres_cdl == 0 for all 3 categories. 
bias %>% filter(is.na(cdl_share))  #it is only 1 county, seemingly irrelevant for our analysis (5 acres census total)
#County concerned: GEOID = 27031

county_level <- bias %>%
  filter(GEOID != "27031") %>% 
  group_by(GEOID) %>%
  summarise(weighted_accuracy_baseline = sum(weighted_accuracy_baseline, na.rm = TRUE),
            accuracy_baseline = mean(accuracy_baseline, na.rm = T), .groups = "drop")

mean(county_level$weighted_accuracy_baseline, na.rm = TRUE)
mean(county_level$accuracy_baseline, na.rm = TRUE)

#That was just to get a magnitude intuition.

# Flag counties where Census reports zero across all 3 categories
# Do NOT remove as zero Census may reflect suppression or land use/cover mismatch
# rather than true absence of crops. Models may still improve CDL for these counties.
bias <- bias %>%
  group_by(GEOID) %>%
  mutate(zero_census_total = sum(acres_census, na.rm = TRUE) == 0) %>%
  ungroup()

zero_counties_length <- bias %>% 
  filter(zero_census_total) %>% 
  pull(GEOID) %>% 
  n_distinct()
#We identify 20 counties with zero acres (Census) accross the 3 categories. 

# Save list of flagged counties for reference in downstream scripts
zero_census_counties <- bias %>%
  filter(zero_census_total) %>%
  distinct(GEOID)

write_csv(zero_census_counties, ZERO_CENSUS_PATH)
write_csv(bias, BASELINE_BIAS_PATH)

