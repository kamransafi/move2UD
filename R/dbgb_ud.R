#' Dynamic Bivariate Gaussian Bridge — Utilisation Distribution
#'
#' Compute a utilisation distribution using the dynamic bivariate Gaussian
#' bridge model, which decomposes movement variance into parallel (along
#' the direction of travel) and orthogonal components.
#'
#' @param object A `move2` object (single track, projected CRS) or an
#'   `mt_dbgb_variance` object from [mt_dbgb_variance()].
#' @param raster A `terra::SpatRaster` defining the output grid, or a
#'   numeric cell size, or `NULL` to auto-compute.
#' @param location_error Numeric scalar or vector of location errors.
#' @param margin Integer (odd). Margin for variance window.
#' @param window_size Integer (odd). Window size for variance estimation.
#' @param ext Numeric. Extension factor for bounding box.
#' @param dim_size Integer. Cells along longest dimension (when raster is NULL).
#' @param time_step Numeric. Integration time step (minutes).
#'
#' @return A `terra::SpatRaster` containing the UD (values sum to 1).
#'
#' @references
#' Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian
#' bridges: directional factorization of diffusion in Brownian bridge
#' models. *Movement Ecology*, 2(1), 5.
#'
#' @export
mt_dbgb_ud <- function(object,
                        raster = NULL,
                        location_error,
                        margin = 15,
                        window_size = 31,
                        ext = 0.5,
                        dim_size = 100,
                        time_step = NULL) {
  UseMethod("mt_dbgb_ud")
}

#' @export
mt_dbgb_ud.move2 <- function(object,
                               raster = NULL,
                               location_error,
                               margin = 15,
                               window_size = 31,
                               ext = 0.5,
                               dim_size = 100,
                               time_step = NULL) {
  var_obj <- mt_dbgb_variance(object, location_error = location_error,
                               margin = margin, window_size = window_size)
  mt_dbgb_ud(var_obj, raster = raster, location_error = location_error,
              ext = ext, dim_size = dim_size, time_step = time_step)
}

#' @export
mt_dbgb_ud.mt_dbgb_variance <- function(object,
                                          raster = NULL,
                                          location_error,
                                          ext = 0.5,
                                          dim_size = 100,
                                          time_step = NULL,
                                          ...) {
  td <- object$track_data
  n <- td$n_locs
  location_error <- .expand_loc_error(location_error, n)

  points_interest <- object$seg_interest | rev(object$seg_interest)

  if (is.null(raster)) {
    pts <- sf::st_as_sf(
      data.frame(x = td$x, y = td$y),
      coords = c("x", "y"), crs = td$crs
    )
    raster <- .make_raster(pts, dim_size = dim_size, ext = ext)
  } else if (is.numeric(raster) && length(raster) == 1) {
    pts <- sf::st_as_sf(
      data.frame(x = td$x, y = td$y),
      coords = c("x", "y"), crs = td$crs
    )
    raster <- .make_raster(pts, cell_size = raster, ext = ext)
  }

  t_mins <- td$time_mins
  if (is.null(time_step)) {
    time_step <- min(diff(t_mins)) / 20.1
  }

  para_sd <- object$para_sd
  orth_sd <- object$orth_sd
  para_sd[is.na(para_sd)] <- 0
  orth_sd[is.na(orth_sd)] <- 0

  x_grid <- terra::xFromCol(raster, 1:terra::ncol(raster))
  y_grid <- sort(unique(terra::yFromRow(raster, 1:terra::nrow(raster))))

  ans <- .Call(
    "bgb_omp",
    td$x[points_interest], td$y[points_interest],
    para_sd[points_interest], orth_sd[points_interest],
    t_mins[points_interest],
    rep(location_error, length.out = n)[points_interest],
    x_grid, y_grid, time_step, 5
  )

  ans <- ans / sum(ans)
  terra::values(raster) <- ans

  raster
}
