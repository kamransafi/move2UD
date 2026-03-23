# Helper script to generate reference test fixtures from the move package.
# Run this ONCE on a machine with both `move` and `move2` installed:
#
#   setwd("/home/kami/Documents/Research/Projects/move2UD")
#   source("tests/testthat/helper-generate-fixtures.R")
#   generate_fixtures()
#
# This creates .rds files in tests/testthat/fixtures/ that the tests
# compare against.

generate_fixtures <- function() {
  if (!requireNamespace("move", quietly = TRUE)) {
    stop("The 'move' package must be installed to generate fixtures.")
  }
  if (!requireNamespace("move2", quietly = TRUE)) {
    stop("The 'move2' package must be installed to generate fixtures.")
  }

  library(move)
  library(move2)
  library(sf)

  fixture_dir <- file.path("tests", "testthat", "fixtures")
  if (!dir.exists(fixture_dir)) dir.create(fixture_dir, recursive = TRUE)

  # --- Load Leroy via move (the old way) ---
  leroy_move <- move::move(system.file("extdata", "leroy.csv.gz", package = "move"))
  leroy_proj_move <- spTransform(leroy_move, center = TRUE)

  # --- Create equivalent move2 object ---
  # Read same CSV, build move2 object manually to ensure identical data
  leroy_df <- read.csv(system.file("extdata", "leroy.csv.gz", package = "move"),
                        as.is = TRUE)
  leroy_df <- leroy_df[!is.na(leroy_df$location.long) &
                          !is.na(leroy_df$location.lat), ]
  leroy_df$timestamp <- as.POSIXct(leroy_df$timestamp,
                                     format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  leroy_sf <- st_as_sf(leroy_df, coords = c("location.long", "location.lat"),
                         crs = 4326)
  leroy_move2 <- mt_as_move2(leroy_sf,
                               time_column = "timestamp",
                               track_id_column = "individual.local.identifier")

  # Project using AEQD centred on track (same as spTransform(center=T))
  bb <- st_bbox(leroy_move2)
  centre_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  centre_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  aeqd_crs <- st_crs(paste0("+proj=aeqd +lon_0=", centre_lon,
                              " +lat_0=", centre_lat, " +units=m"))
  leroy_move2_proj <- st_transform(leroy_move2, aeqd_crs)

  # --- dBBMM variance ---
  cat("Computing dBBMM variance with move...\n")
  move_var <- brownian.motion.variance.dyn(leroy_proj_move,
                                            location.error = 25,
                                            window.size = 31,
                                            margin = 11)

  ref_dbbmm_var <- list(
    variance = move::getMotionVariance(move_var),
    window_size = 31,
    margin = 11,
    location_error = 25
  )
  saveRDS(ref_dbbmm_var, file.path(fixture_dir, "ref_dbbmm_var.rds"))
  cat("  Saved ref_dbbmm_var.rds\n")

  # --- dBGB variance ---
  cat("Computing dBGB variance with move...\n")
  suppressWarnings(
    move_bgb_var <- move::dynBGBvariance(leroy_proj_move, locErr = 25,
                                          windowSize = 31, margin = 15)
  )

  ref_dbgb_var <- list(
    para_var = move::getMotionVariance(move_bgb_var)[, "para"],
    orth_var = move::getMotionVariance(move_bgb_var)[, "orth"],
    window_size = 31,
    margin = 15,
    location_error = 25
  )
  saveRDS(ref_dbgb_var, file.path(fixture_dir, "ref_dbgb_var.rds"))
  cat("  Saved ref_dbgb_var.rds\n")

  # --- Save the projected move2 object ---
  saveRDS(leroy_move2_proj, file.path(fixture_dir, "leroy_projected.rds"))
  cat("  Saved leroy_projected.rds\n")

  cat("Done. Fixtures saved to", fixture_dir, "\n")
}

# Only run if sourced directly
if (sys.nframe() == 0) {
  generate_fixtures()
}
