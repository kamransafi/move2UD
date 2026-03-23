test_that("dbgb_variance_dyn returns expected structure", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "Leroy", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  result <- dbgb_variance_dyn(leroy_proj, loc_err = 25,
                               margin = 15, window_size = 31)

  expect_s3_class(result, "dbgb_var")
  expect_length(result$para_sd, nrow(leroy_proj))
  expect_length(result$orth_sd, nrow(leroy_proj))
  expect_true(any(!is.na(result$para_sd)))
  expect_true(any(!is.na(result$orth_sd)))
  expect_true(all(result$para_sd[!is.na(result$para_sd)] >= 0))
  expect_true(all(result$orth_sd[!is.na(result$orth_sd)] >= 0))
})

test_that("get_motion_variance works on dbgb_var", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "Leroy", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  var_obj <- dbgb_variance_dyn(leroy_proj, loc_err = 25,
                                margin = 15, window_size = 31)
  mv <- get_motion_variance(var_obj)
  expect_true(is.data.frame(mv))
  expect_named(mv, c("para", "orth"))
  expect_equal(mv$para, var_obj$para_sd^2)
  expect_equal(mv$orth, var_obj$orth_sd^2)
})

test_that("dbgb_ud returns a SpatRaster that sums to 1", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "Leroy", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  ud <- dbgb_ud(leroy_proj, loc_err = 25,
                 margin = 15, window_size = 31,
                 dim_size = 50)

  expect_s4_class(ud, "SpatRaster")
  expect_true(abs(sum(terra::values(ud), na.rm = TRUE) - 1) < 0.01)
})

test_that("dbgb_variance_dyn matches move package output", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")

  fixture_file <- test_path("fixtures", "ref_dbgb_var.rds")
  skip_if_not(file.exists(fixture_file),
              "Reference fixtures not generated. Run helper-generate-fixtures.R first.")

  library(move2)
  library(sf)

  ref <- readRDS(fixture_file)
  leroy_proj <- readRDS(test_path("fixtures", "leroy_projected.rds"))

  result <- dbgb_variance_dyn(leroy_proj,
                               loc_err = ref$location_error,
                               margin = ref$margin,
                               window_size = ref$window_size)

  mv <- get_motion_variance(result)
  valid_para <- !is.na(ref$para_var) & !is.na(mv$para)
  valid_orth <- !is.na(ref$orth_var) & !is.na(mv$orth)

  expect_equal(mv$para[valid_para], ref$para_var[valid_para],
               tolerance = 1e-4)
  expect_equal(mv$orth[valid_orth], ref$orth_var[valid_orth],
               tolerance = 1e-4)
})
