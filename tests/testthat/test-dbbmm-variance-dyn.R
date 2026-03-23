test_that("dbbmm_variance_dyn rejects lon/lat input", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  fishers <- fishers[mt_track_id(fishers) == "F1", ]

  expect_error(
    dbbmm_variance_dyn(fishers, location_error = 25,
                        window_size = 31, margin = 11),
    "projected CRS"
  )
})

test_that("dbbmm_variance_dyn works on projected data", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]

  # Project to AEQD centred on track
  bb <- st_bbox(leroy)
  crs_aeqd <- st_crs(paste0("+proj=aeqd +lon_0=",
                              (bb["xmin"] + bb["xmax"]) / 2,
                              " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2,
                              " +units=m"))
  leroy_proj <- st_transform(leroy, crs_aeqd)

  result <- dbbmm_variance_dyn(leroy_proj, location_error = 25,
                                window_size = 31, margin = 11)

  expect_s3_class(result, "dbbmm_var")
  expect_length(result$variance, nrow(leroy_proj))
  expect_true(any(!is.na(result$variance)))
  expect_true(all(result$variance[!is.na(result$variance)] >= 0))
  expect_true(any(result$interest))
  expect_equal(result$window_size, 31)
  expect_equal(result$margin, 11)
})

test_that("dbbmm_variance_dyn errors on even window_size", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    dbbmm_variance_dyn(leroy_proj, location_error = 25,
                        window_size = 30, margin = 11),
    "must both be odd"
  )
})

test_that("dbbmm_variance_dyn errors when window larger than track", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  # Take only a few points
  tiny <- fishers[1:10, ]
  bb <- st_bbox(tiny)
  tiny_proj <- st_transform(tiny, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    dbbmm_variance_dyn(tiny_proj, location_error = 25,
                        window_size = 31, margin = 11),
    "cannot be larger"
  )
})

test_that("get_motion_variance works on dbbmm_var", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  leroy <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(leroy)
  leroy_proj <- st_transform(leroy, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  var_obj <- dbbmm_variance_dyn(leroy_proj, location_error = 25,
                                 window_size = 31, margin = 11)
  mv <- get_motion_variance(var_obj)
  expect_identical(mv, var_obj$variance)
})

test_that("dbbmm_variance_dyn matches move package output", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")

  fixture_file <- test_path("fixtures", "ref_dbbmm_var.rds")
  skip_if_not(file.exists(fixture_file),
              "Reference fixtures not generated. Run helper-generate-fixtures.R first.")

  library(move2)
  library(sf)

  ref <- readRDS(fixture_file)
  leroy_proj <- readRDS(test_path("fixtures", "leroy_projected.rds"))

  result <- dbbmm_variance_dyn(leroy_proj,
                                location_error = ref$location_error,
                                window_size = ref$window_size,
                                margin = ref$margin)

  # Compare variance vectors — close but not identical due to
  # small differences in projection centering between sp::spTransform(center=T)
  # and sf::st_transform with manual AEQD
  valid <- !is.na(ref$variance) & !is.na(result$variance)
  # Use relative tolerance: values should agree within ~1%
  expect_equal(result$variance[valid], ref$variance[valid],
               tolerance = 0.01)
})
