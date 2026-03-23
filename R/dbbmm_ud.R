#' Dynamic Brownian Bridge Movement Model — Utilisation Distribution
#'
#' Compute a utilisation distribution (UD) from a movement track using
#' the dynamic Brownian bridge movement model.
#'
#' @param object A `move2` object (single or multi-track, projected CRS),
#'   an `mt_dbbmm_variance` object, or a named list of `mt_dbbmm_variance`
#'   objects (as returned by [mt_dbbmm_variance()] for multi-track input).
#' @param raster A `terra::SpatRaster` defining the output grid, or a
#'   numeric scalar giving the cell size in map units, or `NULL` to
#'   auto-compute.
#' @param location_error Numeric scalar or vector of location errors.
#' @param margin Integer (odd). Margin for variance estimation.
#' @param window_size Integer (odd). Window size for variance estimation.
#' @param ext Numeric. Extension factor for the bounding box.
#' @param dim_size Integer. Number of cells along the longest dimension.
#' @param time_step Numeric. Time step for integration (in minutes).
#' @param verbose Logical. Print computational size estimate.
#'
#' @return For single-track input: a `terra::SpatRaster` (values sum to 1).
#'   For multi-track input: a multi-layer `terra::SpatRaster` with one
#'   layer per track, all on a common grid. Each layer sums to 1.
#'
#' @export
mt_dbbmm_ud <- function(object,
                         raster = NULL,
                         location_error,
                         margin = 11,
                         window_size = 31,
                         ext = 0.5,
                         dim_size = 100,
                         time_step = NULL,
                         verbose = TRUE) {
  UseMethod("mt_dbbmm_ud")
}

#' @export
mt_dbbmm_ud.move2 <- function(object,
                                raster = NULL,
                                location_error,
                                margin = 11,
                                window_size = 31,
                                ext = 0.5,
                                dim_size = 100,
                                time_step = NULL,
                                verbose = TRUE) {
  var_obj <- mt_dbbmm_variance(object, location_error = location_error,
                                window_size = window_size, margin = margin)
  mt_dbbmm_ud(var_obj, raster = raster, location_error = location_error,
               ext = ext, dim_size = dim_size, time_step = time_step,
               verbose = verbose)
}

#' @export
mt_dbbmm_ud.list <- function(object,
                               raster = NULL,
                               location_error,
                               ext = 0.5,
                               dim_size = 100,
                               time_step = NULL,
                               verbose = TRUE,
                               ...) {
  # Validate: must be a named list of mt_dbbmm_variance objects
  if (!all(vapply(object, inherits, logical(1), "mt_dbbmm_variance"))) {
    stop("List must contain mt_dbbmm_variance objects.", call. = FALSE)
  }

  # Create a common grid from all tracks if not provided
  if (is.null(raster) || (is.numeric(raster) && length(raster) == 1)) {
    all_x <- unlist(lapply(object, function(v) v$track_data$x))
    all_y <- unlist(lapply(object, function(v) v$track_data$y))
    crs <- object[[1]]$track_data$crs
    pts <- sf::st_as_sf(data.frame(x = all_x, y = all_y),
                         coords = c("x", "y"), crs = crs)
    if (is.numeric(raster)) {
      common_raster <- .make_raster(pts, cell_size = raster, ext = ext)
    } else {
      common_raster <- .make_raster(pts, dim_size = dim_size, ext = ext)
    }
  } else {
    common_raster <- raster
  }

  # Compute UD for each track on the common grid
  layers <- lapply(names(object), function(nm) {
    if (verbose) message("Computing UD for track: ", nm)
    mt_dbbmm_ud(object[[nm]], raster = common_raster,
                 location_error = location_error,
                 ext = ext, dim_size = dim_size,
                 time_step = time_step, verbose = FALSE)
  })

  # Stack into multi-layer SpatRaster
  stk <- terra::rast(layers)
  names(stk) <- names(object)
  stk
}

#' @export
mt_dbbmm_ud.mt_dbbmm_variance <- function(object,
                                            raster = NULL,
                                            location_error,
                                            ext = 0.5,
                                            dim_size = 100,
                                            time_step = NULL,
                                            verbose = TRUE,
                                            ...) {
  td <- object$track_data
  n <- td$n_locs
  location_error <- .expand_loc_error(location_error, n)

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

  time_lag <- c(diff(td$time_mins), 0)

  if (is.null(time_step)) {
    time_step <- min(time_lag[-length(time_lag)]) / 15
  }

  if (verbose) {
    comp_size <- terra::ncell(raster) * (sum(time_lag[object$interest]) / time_step)
    message(sprintf("Computational size: %.1e", comp_size))
  }

  x_grid <- terra::xFromCol(raster, 1:terra::ncol(raster))
  y_grid <- terra::yFromRow(raster, terra::nrow(raster):1)

  sigma <- c(object$variance, 0)
  sigma[is.na(sigma)] <- 0

  ans <- .Call(
    "dbbmm2_omp",
    td$x, td$y, sigma,
    (td$time_mins - min(td$time_mins)),
    location_error, x_grid, y_grid,
    time_step, 4, object$interest
  )

  total <- sum(ans)
  if (total == 0 || !is.finite(total)) {
    stop("UD computation produced zero or non-finite values. ",
         "The raster extent may not overlap the track, or ext may be too large.",
         call. = FALSE)
  }
  ans <- ans / total
  terra::values(raster) <- ans

  raster
}
