test_that("delta_para_orth with movement along direction gives pure parallel", {
  # mu at origin, direction along x-axis, point along x-axis
  mu <- matrix(c(0, 0), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(5, 0), ncol = 2)

  result <- move2UD:::delta_para_orth(mu, direction, point)
  expect_equal(result[, "deltaOrth"], 0, tolerance = 1e-10)
  expect_equal(result[, "deltaPara"], 5, tolerance = 1e-10)
})

test_that("delta_para_orth with movement perpendicular to direction gives pure orthogonal", {
  # mu at origin, direction along x-axis, point along y-axis
  mu <- matrix(c(0, 0), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(0, 5), ncol = 2)

  result <- move2UD:::delta_para_orth(mu, direction, point)
  expect_equal(result[, "deltaPara"], 0, tolerance = 1e-10)
  expect_equal(result[, "deltaOrth"], 5, tolerance = 1e-10)
})

test_that("delta_para_orth at 45 degrees splits equally", {
  mu <- matrix(c(0, 0), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(5, 5), ncol = 2)

  result <- move2UD:::delta_para_orth(mu, direction, point)
  # At 45 degrees, para and orth should be equal
  expect_equal(result[, "deltaPara"], result[, "deltaOrth"],
               tolerance = 1e-10)
  # And the total distance should be sqrt(50) ≈ 7.07
  total <- sqrt(result[, "deltaPara"]^2 + result[, "deltaOrth"]^2)
  expect_equal(total, sqrt(50), tolerance = 1e-10)
})

test_that("delta_para_orth handles multiple points", {
  mu <- matrix(c(0, 0, 0, 0), ncol = 2)
  direction <- matrix(c(10, 0, 10, 0), ncol = 2)
  point <- matrix(c(5, 0, 0, 5), ncol = 2)

  result <- move2UD:::delta_para_orth(mu, direction, point)
  expect_equal(nrow(result), 2)
  # First point: pure parallel
  expect_equal(result[1, "deltaOrth"], 0, tolerance = 1e-10)
  # Second point: pure orthogonal
  expect_equal(result[2, "deltaPara"], 0, tolerance = 1e-10)
})

test_that("delta_para_orth handles same location (mu == point)", {
  mu <- matrix(c(5, 5), ncol = 2)
  direction <- matrix(c(10, 0), ncol = 2)
  point <- matrix(c(5, 5), ncol = 2)

  result <- move2UD:::delta_para_orth(mu, direction, point)
  expect_equal(result[, "deltaPara"], 0)
  expect_equal(result[, "deltaOrth"], 0)
})
