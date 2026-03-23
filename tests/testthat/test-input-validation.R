test_that("multi-track input is rejected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  # This has 8 tracks
  bb <- st_bbox(fishers)
  fishers_proj <- st_transform(fishers, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    mt_dbbmm_variance(fishers_proj, location_error = 25,
                       window_size = 31, margin = 11),
    "single track"
  )
})

test_that("empty geometries are rejected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  # Keep empty geometries — F1 has some
  f1 <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(f1[!st_is_empty(f1), ])
  f1_proj <- st_transform(f1, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  expect_error(
    mt_dbbmm_variance(f1_proj, location_error = 25,
                       window_size = 31, margin = 11),
    "empty geometries"
  )
})

test_that("lon/lat input is rejected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  f1 <- fishers[mt_track_id(fishers) == "F1", ]

  expect_error(
    mt_dbbmm_variance(f1, location_error = 25,
                       window_size = 31, margin = 11),
    "projected CRS"
  )
})

test_that("window_size < 2*margin+1 is rejected", {
  skip_if_not_installed("move2")
  skip_if_not_installed("sf")
  library(move2)
  library(sf)

  fishers <- mt_read(mt_example())
  fishers <- fishers[!st_is_empty(fishers), ]
  f1 <- fishers[mt_track_id(fishers) == "F1", ]
  bb <- st_bbox(f1)
  f1_proj <- st_transform(f1, st_crs(paste0(
    "+proj=aeqd +lon_0=", (bb["xmin"] + bb["xmax"]) / 2,
    " +lat_0=", (bb["ymin"] + bb["ymax"]) / 2, " +units=m")))

  # margin=5, window=9 → 2*5+1=11 > 9
  expect_error(
    mt_dbbmm_variance(f1_proj, location_error = 25,
                       window_size = 9, margin = 5),
    "2 \\* margin"
  )
})

test_that("non-move2 input is rejected", {
  expect_error(
    mt_dbbmm_variance(data.frame(x = 1:10), location_error = 25,
                       window_size = 5, margin = 3),
    "move2 object"
  )
})
