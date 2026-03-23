# Helper script to generate reference test fixtures from the move package.
# Run this ONCE on a machine with both `move` and `move2` installed:
#
#   source("tests/testthat/helper-generate-fixtures.R")
#
# This creates .rds files in tests/testthat/fixtures/ that the tests
# compare against. These fixtures freeze the "correct" output from the
# original move package so we can verify numerical equivalence.

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

  # --- Load Leroy data via move ---
  leroy_move <- move(system.file("extdata", "leroy.csv.gz", package = "move"))
  leroy_proj <- spTransform(leroy_move, center = TRUE)

  # --- Load same data via move2 ---
  leroy_move2 <- mt_read(mt_example())
  # Filter to Leroy only
  leroy_move2 <- leroy_move2[mt_track_id(leroy_move2) == "Leroy", ]
  # Remove empty geometries
  leroy_move2 <- leroy_move2[!st_is_empty(leroy_move2), ]
  # Project using AEQD centred on the track
  bb <- st_bbox(leroy_move2)
  centre_lon <- (bb["xmin"] + bb["xmax"]) / 2
  centre_lat <- (bb["ymin"] + bb["ymax"]) / 2
  aeqd_crs <- st_crs(paste0("+proj=aeqd +lon_0=", centre_lon,
                              " +lat_0=", centre_lat, " +units=m"))
  leroy_move2_proj <- st_transform(leroy_move2, aeqd_crs)

  # --- dBBMM variance ---
  cat("Computing dBBMM variance with move...\n")
  move_var <- brownian.motion.variance.dyn(leroy_proj,
                                            location.error = 25,
                                            window.size = 31,
                                            margin = 11)

  ref_dbbmm_var <- list(
    variance = getMotionVariance(move_var),
    window_size = 31,
    margin = 11,
    location_error = 25
  )
  saveRDS(ref_dbbmm_var, file.path(fixture_dir, "ref_dbbmm_var.rds"))

  # --- dBGB variance ---
  cat("Computing dBGB variance with move...\n")
  suppressWarnings(
    move_bgb_var <- dynBGBvariance(leroy_proj, locErr = 25,
                                    windowSize = 31, margin = 15)
  )

  ref_dbgb_var <- list(
    para_var = getMotionVariance(move_bgb_var)[, "para"],
    orth_var = getMotionVariance(move_bgb_var)[, "orth"],
    window_size = 31,
    margin = 15,
    location_error = 25
  )
  saveRDS(ref_dbgb_var, file.path(fixture_dir, "ref_dbgb_var.rds"))

  # --- Save the projected move2 object for test input ---
  saveRDS(leroy_move2_proj, file.path(fixture_dir, "leroy_projected.rds"))

  # --- Save coordinates and timestamps for minimal reproducibility ---
  coords <- sp::coordinates(leroy_proj)
  ref_track <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    time_mins = as.numeric(move::timestamps(leroy_proj)) / 60
  )
  saveRDS(ref_track, file.path(fixture_dir, "ref_track_data.rds"))

  cat("Fixtures saved to", fixture_dir, "\n")
}

# Only run if sourced directly, not when loaded by testthat
if (sys.nframe() == 0) {
  generate_fixtures()
}
