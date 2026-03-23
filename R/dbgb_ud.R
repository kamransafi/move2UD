#' Dynamic Bivariate Gaussian Bridge — Utilisation Distribution
#'
#' Compute a utilisation distribution using the dynamic bivariate Gaussian
#' bridge model, which decomposes movement variance into parallel (along
#' the direction of travel) and orthogonal components.
#'
#' @param object A `move2` object (single track, projected CRS) or a
#'   `dbgb_var` object from [dbgb_variance_dyn()].
#' @param raster A `terra::SpatRaster` defining the output grid, or a
#'   numeric cell size, or `NULL` to auto-compute.
#' @param loc_err Numeric scalar or vector of location errors.
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
dbgb_ud <- function(object,
                    raster = NULL,
                    loc_err,
                    margin = 15,
                    window_size = 31,
                    ext = 0.5,
                    dim_size = 100,
                    time_step = NULL) {
  UseMethod("dbgb_ud")
}

#' @export
dbgb_ud.move2 <- function(object,
                           raster = NULL,
                           loc_err,
                           margin = 15,
                           window_size = 31,
                           ext = 0.5,
                           dim_size = 100,
                           time_step = NULL) {
  var_obj <- dbgb_variance_dyn(object, loc_err = loc_err,
                                margin = margin, window_size = window_size)
  dbgb_ud(var_obj, raster = raster, loc_err = loc_err,
          ext = ext, dim_size = dim_size, time_step = time_step)
}

#' @export
dbgb_ud.dbgb_var <- function(object,
                              raster = NULL,
                              loc_err,
                              ext = 0.5,
                              dim_size = 100,
                              time_step = NULL,
                              ...) {
  td <- object$track_data
  n <- td$n_locs
  loc_err <- .expand_loc_error(loc_err, n)

  # Points of interest: segments where both the segment and its reverse are valid
  points_interest <- object$seg_interest | rev(object$seg_interest)

  # Create raster if needed
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

  # Prepare variance vectors — replace NA with 0
  para_sd <- object$para_sd
  orth_sd <- object$orth_sd
  para_sd[is.na(para_sd)] <- 0
  orth_sd[is.na(orth_sd)] <- 0

  # Grid coordinates
  x_grid <- terra::xFromCol(raster, 1:terra::ncol(raster))
  y_grid <- sort(unique(terra::yFromRow(raster, 1:terra::nrow(raster))))

  # Call the C kernel
  ans <- .Call(
    "bgb",
    td$x[points_interest],
    td$y[points_interest],
    para_sd[points_interest],
    orth_sd[points_interest],
    t_mins[points_interest],
    rep(loc_err, length.out = n)[points_interest],
    x_grid,
    y_grid,
    time_step,
    5  # number of SDs for integration distance
  )

  # Normalise
  ans <- ans / sum(ans)
  terra::values(raster) <- ans

  raster
}
