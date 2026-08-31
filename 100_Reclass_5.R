################################################################################
# 11_post_processing_models.R  (RESUMABLE VERSION)
#
# Two post-classification correction models, scored per county vs. Census:
#   - MMU sweep: patch removal + modal fill over a grid of params (mmu_model_grid)
#   - CSB field-modal: reclassify to modal CDL category within CSB polygons
#
# Also computes quality_ratio per (GEOID, model): confidence of reclassified
# pixels vs. all cropland pixels, to flag coincidental "improvements".
#
# best_model_per_county selection happens in 12_merge_post_processing_parts.R,
# after all parts are combined -- not here.
#
# WHICHPART: set manually below. 0 = test (N_TEST_COUNTIES counties).
# 1..N_PARTS = that slice of county_lookup (~3000 counties / N_PARTS each).
#
# ---------------------------------------------------------------------------
# RESUMABILITY -- what changed vs. the original version
# ---------------------------------------------------------------------------
# Previously, all 311+ counties' results were held in memory across the
# whole run and only written to disk in one big CSV at the very end. If the
# cluster crashed (e.g. the "unserialize: error reading from connection"
# socket error) partway through, EVERYTHING for that part was lost and the
# whole part had to be rerun from county #1.
#
# Now, each county's MMU results and each (state, county) pair's CSB result
# is written to its OWN small .rds file the moment it finishes:
#   outputs/<VERSION>/parts/part_<N>/mmu_by_county/<GEOID>.rds
#   outputs/<VERSION>/parts/part_<N>/csb_by_county/<GEOID>.rds
#
# Both process_county_mmu() and the per-county CSB step check FIRST whether
# that GEOID's .rds already exists -- if so, they skip it immediately
# (no raster read, no computation) rather than redoing work. This means:
#   - If the job crashes, just resubmit the SAME script (same WHICHPART) --
#     it will skip every county already done and pick up where it left off.
#   - The parLapplyLB workers now return NULL instead of a full results
#     tibble, so far less data crosses the socket per county -- this should
#     also make the "unserialize" crash itself less likely to happen again.
#   - The final combine step at the bottom no longer relies on what THIS
#     run's parLapplyLB call returned -- it just reads every .rds file
#     sitting in mmu_by_county/ and csb_by_county/, so it's correct even if
#     those files were written across several different (crashed/resumed)
#     runs.
#
# Inputs:  outputs/<VERSION>/class_and_conf/<STATEFP>/Conf_stacked_<YEAR>_<GEOID>.tif
#          outputs/<VERSION>/baseline_measures.csv
#          data/county_lookup.csv, data/SF/counties_2016.rds, CSB_GDB
#
# Outputs: outputs/<VERSION>/parts/part_<WHICHPART>/mmu_by_county/<GEOID>.rds
#          outputs/<VERSION>/parts/part_<WHICHPART>/csb_by_county/<GEOID>.rds
#          outputs/<VERSION>/parts/part_<WHICHPART>/models_processing_summary.csv
#          outputs/<VERSION>/parts/part_<WHICHPART>/models_processing_details.csv
#          figures/<VERSION>/parts/part_<WHICHPART>/post_processing/*.png (if MAKE_PLOTS)
################################################################################
rm(list = ls())
setwd("/users/rperilhou/extra_years")
library(tidyverse)
library(terra)
#library(landscapemetrics)
library(sf)

################################################################################
# Configuration
################################################################################
PIXEL_ACRES    <- 0.222395
TARGET_YEAR    <- 2012
CROPLAND_CODES <- c(1L, 2L, 3L)
CODE_TO_NAME   <- c("1" = "GM", "2" = "Tolerant", "3" = "Vulnerable")

VERSION <- "v2"

################################################################################
# WHICHPART / chunking configuration
################################################################################
N_PARTS         <- 10   # number of full-scale chunks (~300 counties each on ~3000 counties)
N_TEST_COUNTIES <- 2    # how many counties WHICHPART = 0 processes

WHICHPART <- 5
PART_TAG <- paste0("part_", WHICHPART)

MAKE_PLOTS <- TRUE
################################################################################
# PATHS BUILDER -- part-specific output locations
################################################################################
CONF_STACKED_DIR <- file.path("outputs", VERSION, "class_and_conf")
CONF_STACKED_PATH <- function(geoid, statefp) {
  file.path(CONF_STACKED_DIR, statefp,
            paste0("Conf_stacked_", TARGET_YEAR, "_", geoid, ".tif"))
}
BASELINE_MEASURES_PATH <- file.path("outputs", VERSION, "baseline_measures.csv")

# figures + summary/detail outputs all nest under a part_<WHICHPART> folder
FIGURES_DIR <- file.path("figures", VERSION, "parts", PART_TAG, "post_processing")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

PART_OUTPUT_DIR <- file.path("outputs", VERSION, "parts", PART_TAG)
dir.create(PART_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# per-county checkpoint directories -- one small .rds per GEOID, written as
# soon as that county finishes, checked for existence before (re)computing
MMU_COUNTY_DIR <- file.path(PART_OUTPUT_DIR, "mmu_by_county")
CSB_COUNTY_DIR <- file.path(PART_OUTPUT_DIR, "csb_by_county")
dir.create(MMU_COUNTY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CSB_COUNTY_DIR, recursive = TRUE, showWarnings = FALSE)

mmu_county_path <- function(geoid) file.path(MMU_COUNTY_DIR, paste0(geoid, ".rds"))
csb_county_path <- function(geoid) file.path(CSB_COUNTY_DIR, paste0(geoid, ".rds"))

COUNTY_LOOKUP_PATH <- "data/county_lookup.csv"

CSB_GDB   <- "data/NationalCSB_2013-2020_rev23/CSB1320.gdb"
CSB_LAYER <- "national1320"

MODELS_SUMMARY_PATH <- file.path(PART_OUTPUT_DIR, "models_processing_summary.csv")
MODELS_DETAILS_PATH <- file.path(PART_OUTPUT_DIR, "models_processing_details.csv")

cat_labels <- c("0" = "NonCrop", "1" = "GM", "2" = "Tolerant",
                "3" = "Vulnerable", "99" = "Unclassified")
cat_colors <- c("0" = "grey85", "1" = "#e41a1c", "2" = "#ff7f00",
                "3" = "#4daf4a", "99" = "#377eb8")
breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 99.5)

################################################################################
# MODEL GRID: MMU parameter combinations
################################################################################
mmu_model_grid <- tribble(
  ~mmu_acres, ~neighbors, ~window, ~max_fill_iter,
  35,          8,          3,       3,
  35,          8,          9,       3,
  35,          4,          3,       3,
  35,          4,          9,       3,
  25,          8,          3,       3,
  25,          8,          9,       3,
  25,          4,          3,       3,
  25,          4,          9,       3,
  15,          8,          3,       3,
  15,          8,          9,       3,
  15,          4,          3,       3,
  15,          4,          9,       3,
  5,           8,          3,       3,
  5,           8,          9,       3,
  5,           4,          3,       3,
  2,           8,          3,       3,
  1,           8,          3,       3,
  1,           8,          3,       1
) %>%
  mutate(model_name = paste0("MMU_", mmu_acres, "ac_nb", neighbors, "_w", window, "_iter", max_fill_iter))

################################################################################
# run one MMU model over a category raster
################################################################################
run_mmu_model <- function(original_category, mmu_acres, neighbors, window, max_fill_iter,
                          model_name, cropland_codes = CROPLAND_CODES,
                          pixel_acres = PIXEL_ACRES) {
  
  output <- terra::deepcopy(original_category)
  
  for (code in cropland_codes) {
    class_raster <- ifel(original_category == code, 1, NA)
    patch_raster <- terra::patches(class_raster, directions = neighbors, zeroAsNA = TRUE)
    
    freq_tbl  <- terra::freq(patch_raster)
    small_ids <- freq_tbl$value[freq_tbl$count * pixel_acres < mmu_acres]
    if (length(small_ids) == 0) next
    
    temp <- ifel(patch_raster %in% small_ids, NA, original_category)
    
    for (iter in 1:max_fill_iter) {
      holes   <- is.na(temp) & !is.na(original_category)
      n_holes <- global(holes, "sum", na.rm = TRUE)[1, 1]
      if (n_holes == 0) break
      temp <- terra::focal(temp, w = window, fun = "modal", na.rm = TRUE, na.policy = "only")
    }
    fill <- temp
    
    output <- ifel(patch_raster %in% small_ids, fill, output)
  }
  
  model_raster <- ifel(is.na(original_category), NA, output)
  
  freq_tbl <- terra::freq(model_raster)
  freq_tbl <- freq_tbl[freq_tbl$value %in% cropland_codes, ]
  model_acres <- tibble(
    Category        = unname(CODE_TO_NAME[as.character(freq_tbl$value)]),
    acres_cdl_model = freq_tbl$count * pixel_acres
  )
  
  list(raster = model_raster, acres = model_acres, model_name = model_name)
}

################################################################################
# run the CSB field-modal model over a category raster
################################################################################
run_csb_model <- function(original_category, county_poly, csb_state,
                          cropland_codes = CROPLAND_CODES, pixel_acres = PIXEL_ACRES) {
  
  county_poly_proj <- st_transform(county_poly, st_crs(csb_state))
  csb_county <- st_filter(csb_state, county_poly_proj)
  
  if (nrow(csb_county) == 0) {
    return(NULL)
  }
  
  csb_vect   <- vect(st_transform(csb_county, crs(original_category)))
  modal_cats <- terra::extract(original_category, csb_vect, fun = "modal", na.rm = TRUE)
  csb_county$modal_cat <- as.integer(modal_cats[[names(original_category)]])
  csb_county$modal_cat[is.na(csb_county$modal_cat)] <- 99L
  
  csb_vect_fb <- vect(st_transform(csb_county, crs(original_category)))
  model_raster <- rasterize(csb_vect_fb, original_category, field = "modal_cat")
  model_raster <- mask(model_raster, original_category)
  model_raster <- ifel(!is.na(original_category) & is.na(model_raster), original_category, model_raster)
  
  freq_tbl <- terra::freq(model_raster)
  freq_tbl <- freq_tbl[freq_tbl$value %in% cropland_codes, ]
  model_acres <- tibble(
    Category        = unname(CODE_TO_NAME[as.character(freq_tbl$value)]),
    acres_cdl_model = freq_tbl$count * pixel_acres
  )
  
  list(raster = model_raster, acres = model_acres, model_name = "CSB_field_modal")
}

################################################################################
# FUNCTION: save a diagnostic plot for one model's output raster
# path: figures/<VERSION>/parts/<PART_TAG>/post_processing/<STATEFP>_<GEOID>_<model_name>.png
################################################################################
plot_model <- function(model_raster, geoid, statefp, model_name, out_dir = FIGURES_DIR) {
  png(file.path(out_dir, paste0(statefp, "_", geoid, "_", model_name, ".png")),
      width = 1400, height = 700, res = 150, type = "cairo")
  par(mfrow = c(1, 1))
  plot(model_raster, breaks = breaks, col = unname(cat_colors), legend = FALSE,
       main = model_name)
  legend("bottomleft", legend = unname(cat_labels), fill = unname(cat_colors),
         cex = 0.7, bg = "white")
  dev.off()
}

################################################################################
# FUNCTION: confidence of reclassified pixels vs. all cropland (shared by
# every model type -- MMU or CSB, since both just produce a category raster
# to compare against the original). One row per (GEOID, model) -- county-
# wide, NOT split by Category.
################################################################################
compute_quality <- function(stacked_raster, model_category_raster, geoid, model_name) {
  orig_cat <- stacked_raster[["category"]]
  conf     <- stacked_raster[["confidence"]]
  
  reclass_ind <- ifel(!is.na(orig_cat) & !is.na(model_category_raster) &
                        orig_cat != model_category_raster, 1L, 0L)
  crop_ind    <- ifel(orig_cat %in% CROPLAND_CODES, 1L, 0L)
  
  z_re <- zonal(conf, reclass_ind, fun = "mean", na.rm = TRUE)
  z_cr <- zonal(conf, crop_ind,    fun = "mean", na.rm = TRUE)
  n_re <- terra::freq(reclass_ind)
  n_reclassified <- sum(n_re$count[n_re$value == 1])
  
  mean_conf_reclassified <- if (1 %in% z_re[[1]]) z_re[z_re[[1]] == 1, 2] else NA_real_
  mean_conf_all_cropland <- if (1 %in% z_cr[[1]]) z_cr[z_cr[[1]] == 1, 2] else NA_real_
  
  tibble(
    GEOID = geoid,
    model = model_name,
    n_reclassified = n_reclassified,
    mean_conf_reclassified = mean_conf_reclassified,
    mean_conf_all_cropland = mean_conf_all_cropland,
    quality_ratio = ifelse(
      !is.na(mean_conf_reclassified) & !is.na(mean_conf_all_cropland) &
        mean_conf_all_cropland != 0,
      mean_conf_reclassified / mean_conf_all_cropland, NA_real_
    )
  )
}

################################################################################
# FUNCTION: bundle a model's raster + acres into the same results/quality
# shape, regardless of which model produced it (MMU or CSB)
################################################################################
score_model <- function(model_output, stacked, baseline, geoid, statefp, make_plots) {
  if (is.null(model_output)) return(list(results = NULL, quality = NULL))
  
  model_name_used <- model_output$model_name
  
  if (make_plots) plot_model(model_output$raster, geoid, statefp, model_name_used)
  
  results <- baseline %>%
    left_join(model_output$acres, by = "Category") %>%
    mutate(
      acres_cdl_model   = replace_na(acres_cdl_model, 0),
      delta_acres_model = acres_cdl_model - acres_census,
      improvement_acres = abs(delta_acres_baseline) - abs(delta_acres_model),
      improvement_pct   = 100 * (abs(delta_acres_baseline) - abs(delta_acres_model)) / abs(delta_acres_baseline),
      model = model_name_used
    )
  
  quality <- compute_quality(stacked, model_output$raster, geoid, model_name_used)
  
  list(results = results, quality = quality)
}

################################################################################
# FUNCTION: MMU sweep for ONE county -- parallelizable at county grain,
# since it needs nothing shared across counties (no state-level read).
# RESUMABLE: checks mmu_county_path(geoid) FIRST -- if it already exists,
# skips all computation for this county entirely. On success, writes the
# county's results/quality to that path before returning.
################################################################################
process_county_mmu <- function(geoid, statefp, mmu_model_grid, baseline_all,
                               make_plots = FALSE) {
  
  out_path <- mmu_county_path(geoid)
  if (file.exists(out_path)) {
    cat("  GEOID", geoid, "MMU already done, skipping.\n")
    return(invisible(NULL))
  }
  
  cat("Processing GEOID (MMU):", geoid, "\n")
  
  raster_path <- CONF_STACKED_PATH(geoid, statefp)
  if (!file.exists(raster_path)) {
    cat("  Missing stacked raster, skipping:", raster_path, "\n")
    return(invisible(NULL))
  }
  
  baseline <- baseline_all %>%
    filter(GEOID == geoid) %>%
    select(GEOID, Category, acres_cdl_baseline, acres_census, delta_acres_baseline)
  
  if (nrow(baseline) == 0) {
    cat("  No baseline row for this GEOID, skipping.\n")
    return(invisible(NULL))
  }
  
  stacked            <- rast(raster_path)
  original_category  <- stacked[["category"]]
  
  county_results  <- list()
  quality_results <- list()
  
  mmu_failed <- FALSE
  for (i in seq_len(nrow(mmu_model_grid))) {
    row <- mmu_model_grid[i, ]
    model_output <- tryCatch(
      run_mmu_model(
        original_category = original_category,
        mmu_acres  = row$mmu_acres, neighbors = row$neighbors,
        window = row$window, max_fill_iter = row$max_fill_iter,
        model_name = row$model_name
      ),
      error = function(e) {
        cat("  GEOID", geoid, "MMU model", row$model_name,
            "FAILED:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(model_output)) {
      mmu_failed <- TRUE
      next  # skip this grid row, keep trying the rest of the grid
    }
    scored <- score_model(model_output, stacked, baseline, geoid, statefp, make_plots)
    county_results[[row$model_name]]  <- scored$results
    quality_results[[row$model_name]] <- scored$quality
  }
  
  if (mmu_failed) {
    # ANY grid row failed -- do NOT write a checkpoint file. This means the
    # NEXT run will see no file for this county and reprocess it entirely
    # (all 18 grid rows again, not just the ones that failed) -- slightly
    # wasteful for the rows that already succeeded, but guarantees a county
    # is never silently left "done" with some models missing. Check
    # mmu_failures.csv to see exactly which (GEOID, model) pairs failed
    # without having to grep the full job log.
    cat("  GEOID", geoid, "MMU: at least one grid row failed, NOT checkpointing ",
        "-- full county will retry on next run.\n")
    return(invisible(NULL))
  }
  
  # write this county's checkpoint file -- only reached when ALL grid rows
  # succeeded. Do this BEFORE returning, so a crash on the NEXT county
  # still leaves this one safely on disk.
  saveRDS(
    list(results = bind_rows(county_results), quality = bind_rows(quality_results)),
    out_path
  )
  
  cat("  GEOID", geoid, "MMU sweep complete, saved to", out_path, "\n")
  # return NULL (not the data) -- the data already lives on disk; this
  # keeps what crosses the parLapplyLB socket per county tiny, which should
  # also reduce the odds of a socket/unserialize crash like the one before
  invisible(NULL)
}

################################################################################
# FUNCTION: CSB field-modal model for an entire STATE -- reads the state's
# CSB polygons ONCE, then loops every county in that state. This is the
# unit meant to be handed to one core (parLapply/future_map over states).
# RESUMABLE: each county's result is checkpointed to csb_county_path(geoid)
# and skipped on subsequent runs if already present -- same logic as the
# MMU function above, just applied per-county inside the per-state loop.
# NOTE: since counties are split across parts by WHICHPART, a given state's
# counties may be spread across multiple parts -- each part still reads and
# filters the CSB polygons for that state independently, scoped to this
# part's counties. Some redundant I/O across parts, but no coordination
# needed between nodes.
################################################################################
process_state_csb <- function(statefp, geoids_in_state, baseline_all, counties_sf,
                              make_plots = FALSE) {
  
  # skip the (possibly expensive) CSB_GDB read entirely if every county in
  # this state (within this part) is already checkpointed
  remaining_geoids <- geoids_in_state[!file.exists(csb_county_path(geoids_in_state))]
  if (length(remaining_geoids) == 0) {
    cat("STATEFP (CSB):", statefp, "-- all", length(geoids_in_state),
        "counties already done, skipping state entirely.\n")
    return(invisible(NULL))
  }
  
  cat("####################################\n")
  cat("STATEFP (CSB):", statefp, "|", length(remaining_geoids), "of",
      length(geoids_in_state), "counties remaining\n")
  cat("####################################\n")
  
  csb_state <- st_read(
    CSB_GDB, layer = CSB_LAYER,
    query = paste0("SELECT * FROM ", CSB_LAYER, " WHERE CSBID LIKE '", statefp, "%'"),
    quiet = TRUE
  )
  
  invisible(lapply(remaining_geoids, function(geoid) {
    out_path <- csb_county_path(geoid)
    # re-check in case of a race with another concurrent job; cheap safety net
    if (file.exists(out_path)) {
      cat("  GEOID", geoid, "CSB already done, skipping.\n")
      return(invisible(NULL))
    }
    
    cat("  Processing GEOID (CSB):", geoid, "\n")
    
    raster_path <- CONF_STACKED_PATH(geoid, statefp)
    if (!file.exists(raster_path)) {
      cat("    Missing stacked raster, skipping:", raster_path, "\n")
      return(invisible(NULL))
    }
    
    baseline <- baseline_all %>%
      filter(GEOID == geoid) %>%
      select(GEOID, Category, acres_cdl_baseline, acres_census, delta_acres_baseline)
    
    if (nrow(baseline) == 0) {
      cat("    No baseline row for this GEOID, skipping.\n")
      return(invisible(NULL))
    }
    
    stacked           <- rast(raster_path)
    original_category <- stacked[["category"]]
    county_poly       <- counties_sf %>% filter(GEOID == geoid)
    
    csb_output <- tryCatch(
      run_csb_model(original_category, county_poly, csb_state),
      error = function(e) {
        cat("    CSB model failed for this county:", conditionMessage(e), "\n")
        NULL
      }
    )
    
    scored <- score_model(csb_output, stacked, baseline, geoid, statefp, make_plots)
    if (is.null(scored$results) && is.null(scored$quality)) return(invisible(NULL))
    
    saveRDS(scored, out_path)
    cat("    GEOID", geoid, "CSB complete, saved to", out_path, "\n")
    invisible(NULL)
  }))
  
  invisible(NULL)
}

################################################################################
# RUN (ANUBIS)
################################################################################

library(parallel)

county_lookup <- read.csv(COUNTY_LOOKUP_PATH, stringsAsFactors = FALSE) %>%
  mutate(
    GEOID   = sprintf("%05d", as.integer(GEOID)),
    STATEFP = sprintf("%02d", as.integer(STATEFP))
  )

baseline_all <- read_csv(BASELINE_MEASURES_PATH, show_col_types = FALSE) %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  rename(acres_cdl_baseline = acres_cdl, delta_acres_baseline = delta_acres)

counties_sf <- readRDS("data/SF/counties_2016.rds")

################################################################################
# Slice county_lookup according to WHICHPART
################################################################################
if (WHICHPART == 0) {
  tasks <- county_lookup %>% slice_head(n = N_TEST_COUNTIES)
  cat("TEST MODE: running on", nrow(tasks), "counties ->",
      paste(tasks$GEOID, collapse = ", "), "\n")
} else {
  part_index <- cut(
    seq_len(nrow(county_lookup)),
    breaks = N_PARTS,
    labels = FALSE
  )
  tasks <- county_lookup[part_index == WHICHPART, ]
  cat("WHICHPART", WHICHPART, "of", N_PARTS, "-> ", nrow(tasks),
      "counties (GEOID range:", first(tasks$GEOID), "-", last(tasks$GEOID), ")\n")
}

if (nrow(tasks) == 0) {
  stop("No counties assigned to WHICHPART = ", WHICHPART, ". Check N_PARTS / county_lookup.")
}

n_mmu_done <- sum(file.exists(mmu_county_path(tasks$GEOID)))
n_csb_done <- sum(file.exists(csb_county_path(tasks$GEOID)))
cat("Resume check: MMU already done for", n_mmu_done, "/", nrow(tasks),
    "counties; CSB already done for", n_csb_done, "/", nrow(tasks), "counties.\n")

################################################################################
# PASS 1: MMU sweep, parallelized one-core-per-GEOID
################################################################################
cat("Starting MMU pass:", nrow(tasks), "counties\n")

source("/softs/R/createCluster.R")
cl <- createCluster()

clusterExport(cl, c(
  "PIXEL_ACRES", "TARGET_YEAR", "CROPLAND_CODES", "CODE_TO_NAME",
  "CONF_STACKED_DIR", "CONF_STACKED_PATH", "FIGURES_DIR",
  "mmu_county_path", "csb_county_path",
  "cat_labels", "cat_colors", "breaks",
  "mmu_model_grid", "run_mmu_model", "plot_model", "compute_quality",
  "score_model", "process_county_mmu",
  "tasks", "baseline_all", "MAKE_PLOTS"
))

invisible(parLapplyLB(cl, seq_len(nrow(tasks)), function(i) {
  library(tidyverse)
  library(terra)
  process_county_mmu(
    geoid          = tasks$GEOID[i],
    statefp        = tasks$STATEFP[i],
    mmu_model_grid = mmu_model_grid,
    baseline_all   = baseline_all,
    make_plots     = MAKE_PLOTS
  )
}))

stopCluster(cl)
cat("MMU pass complete.\n")

################################################################################
# PASS 2: CSB model, parallelized one-core-per-STATEFP (scoped to this part's
# counties only). Resumable per-county inside process_state_csb.
################################################################################
tasks_by_state <- split(tasks$GEOID, tasks$STATEFP)
cat("Starting CSB pass:", length(tasks_by_state), "states (within this part)\n")

cl <- createCluster()

clusterExport(cl, c(
  "PIXEL_ACRES", "TARGET_YEAR", "CROPLAND_CODES", "CODE_TO_NAME",
  "CONF_STACKED_DIR", "CONF_STACKED_PATH", "FIGURES_DIR",
  "CSB_GDB", "CSB_LAYER",
  "mmu_county_path", "csb_county_path",
  "cat_labels", "cat_colors", "breaks",
  "run_csb_model", "plot_model", "compute_quality",
  "score_model", "process_state_csb",
  "tasks_by_state", "baseline_all", "counties_sf", "MAKE_PLOTS"
))

invisible(parLapplyLB(cl, names(tasks_by_state), function(statefp) {
  library(tidyverse)
  library(terra)
  library(sf)
  process_state_csb(
    statefp          = statefp,
    geoids_in_state  = tasks_by_state[[statefp]],
    baseline_all     = baseline_all,
    counties_sf      = counties_sf,
    make_plots       = MAKE_PLOTS
  )
}))

stopCluster(cl)
cat("CSB pass complete.\n")

################################################################################
# COMBINE + WRITE this part's outputs -- reads EVERY .rds checkpoint file
# currently on disk for this part, regardless of which run (this one or a
# previous crashed one) produced it. This is what makes the whole thing
# safe to resume: even if THIS run only processed a handful of remaining
# counties, the combine step below still picks up all of them together.
################################################################################
read_checkpoints <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) return(list(results = tibble(), quality = tibble()))
  all_data <- lapply(files, readRDS)
  list(
    results = bind_rows(lapply(all_data, `[[`, "results")),
    quality = bind_rows(lapply(all_data, `[[`, "quality"))
  )
}

mmu_combined <- read_checkpoints(MMU_COUNTY_DIR)
csb_combined <- read_checkpoints(CSB_COUNTY_DIR)

n_mmu_files <- length(list.files(MMU_COUNTY_DIR, pattern = "\\.rds$"))
n_csb_files <- length(list.files(CSB_COUNTY_DIR, pattern = "\\.rds$"))
cat("Combining checkpoints: ", n_mmu_files, "MMU county files,",
    n_csb_files, "CSB county files.\n")

all_results <- bind_rows(mmu_combined$results, csb_combined$results)
all_quality <- bind_rows(mmu_combined$quality, csb_combined$quality)

models_processing_details <- all_results %>%
  left_join(all_quality, by = c("GEOID", "model"))

models_processing_summary <- models_processing_details %>%
  group_by(GEOID, model) %>%
  summarise(
    bias_base      = sum(abs(delta_acres_baseline), na.rm = TRUE),
    bias_corrected = sum(abs(delta_acres_model), na.rm = TRUE),
    imp_share_pct  = 100 * sum(improvement_acres, na.rm = TRUE) / sum(abs(delta_acres_baseline), na.rm = TRUE),
    quality_ratio  = first(quality_ratio),
    .groups = "drop"
  )

write_csv(models_processing_details, MODELS_DETAILS_PATH)
write_csv(models_processing_summary, MODELS_SUMMARY_PATH)

################################################################################
# INCOMPLETE-COUNTY SUMMARY -- tells you at a glance, without opening any
# log file, whether you need to resubmit this part. A county only ever gets
# a checkpoint file when it fully succeeds (all 18 MMU grid rows, or the
# CSB model), so "missing from disk" == "still needs to (re)run".
################################################################################
missing_mmu <- setdiff(tasks$GEOID, sub("\\.rds$", "", basename(
  list.files(MMU_COUNTY_DIR, pattern = "\\.rds$"))))
missing_csb <- setdiff(tasks$GEOID, sub("\\.rds$", "", basename(
  list.files(CSB_COUNTY_DIR, pattern = "\\.rds$"))))

if (length(missing_mmu) == 0 && length(missing_csb) == 0) {
  cat("\nALL COUNTIES COMPLETE for part", WHICHPART, "-- no need to rerun.\n")
} else {
  cat("\nINCOMPLETE -- rerun this same script (WHICHPART =", WHICHPART,
      ") to retry the counties below (already-done counties will be skipped):\n")
  if (length(missing_mmu) > 0) {
    cat("  MMU missing (", length(missing_mmu), "counties):",
        paste(missing_mmu, collapse = ", "), "\n")
  }
  if (length(missing_csb) > 0) {
    cat("  CSB missing (", length(missing_csb), "counties):",
        paste(missing_csb, collapse = ", "), "\n")
  }
}

cat("\nPart", WHICHPART, "combine step done.\n")
cat("  Details written to:", MODELS_DETAILS_PATH, "\n")
cat("  Summary written to:", MODELS_SUMMARY_PATH, "\n")
cat("  Counties in tasks:", nrow(tasks),
    "| MMU checkpoints on disk:", n_mmu_files,
    "| CSB checkpoints on disk:", n_csb_files, "\n")