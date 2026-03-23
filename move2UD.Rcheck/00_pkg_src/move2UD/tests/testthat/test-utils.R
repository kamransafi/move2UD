test_that(".expand_loc_error works with scalar", {
  result <- move2UD:::.expand_loc_error(25, 10)
  expect_length(result, 10)
  expect_true(all(result == 25))
})

test_that(".expand_loc_error works with vector of correct length", {
  err <- runif(10, 5, 50)
  result <- move2UD:::.expand_loc_error(err, 10)
  expect_identical(result, err)
})

test_that(".expand_loc_error errors on wrong length", {
  expect_error(
    move2UD:::.expand_loc_error(c(10, 20, 30), 10),
    "length 1 or equal"
  )
})

test_that(".expand_loc_error errors on NA", {
  expect_error(
    move2UD:::.expand_loc_error(c(10, NA, 30), 3),
    "must not contain NAs"
  )
})

test_that(".expand_loc_error errors on non-positive", {
  expect_error(
    move2UD:::.expand_loc_error(c(10, 0, 30), 3),
    "must be positive"
  )
  expect_error(
    move2UD:::.expand_loc_error(c(10, -5, 30), 3),
    "must be positive"
  )
})

test_that(".extcalc produces expanded bounding box", {
  skip_if_not_installed("sf")
  library(sf)
  pts <- st_as_sf(data.frame(x = c(0, 100), y = c(0, 100)),
                   coords = c("x", "y"), crs = 32632)
  bb <- move2UD:::.extcalc(pts, ext = 0.3)
  expect_true(bb["xmin"] < 0)
  expect_true(bb["xmax"] > 100)
  expect_true(bb["ymin"] < 0)
  expect_true(bb["ymax"] > 100)
})

test_that(".make_raster creates raster with correct CRS", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  library(sf)
  pts <- st_as_sf(data.frame(x = c(0, 1000), y = c(0, 1000)),
                   coords = c("x", "y"), crs = 32632)
  r <- move2UD:::.make_raster(pts, cell_size = 100, ext = 0.1)
  expect_s4_class(r, "SpatRaster")
  expect_true(grepl("32632", terra::crs(r, describe = TRUE)$code))
})
