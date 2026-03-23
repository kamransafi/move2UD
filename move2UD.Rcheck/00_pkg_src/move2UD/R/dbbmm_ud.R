#' Dynamic Brownian Bridge Movement Model — Utilisation Distribution
#'
#' Compute a utilisation distribution (UD) from a movement track using
#' the dynamic Brownian bridge movement model. Accepts a `move2` object
#' directly (estimates variance internally), a pre-computed variance
#' object, or a named list of variance objects for multi-track input.
#'
#' @param object One of:
#'   \itemize{
#'     \item A `move2` object (single or multi-track, projected CRS).
#'       Variance is estimated internally via [mt_dbbmm_variance()].
#'     \item An `mt_dbbmm_variance` object from [mt_dbbmm_variance()].
#'     \item A named list of `mt_dbbmm_variance` objects (multi-track).
#'   }
#' @param raster A `terra::SpatRaster` defining the output grid, or a
#'   numeric scalar giving the cell size in map units, or `NULL` to
#'   auto-compute from `dim_size` and `ext`. For multi-track input with
#'   `NULL`, a common grid is computed from the combined extent.
#' @param location_error Numeric scalar or vector of location errors in
#'   map units. Must be positive.
#' @param margin Integer (odd). Margin for variance estimation window.
#'   Only used when `object` is a `move2` object.
#' @param window_size Integer (odd). Window size for variance estimation.
#'   Only used when `object` is a `move2` object.
#' @param ext Numeric. Extension factor for the bounding box when
#'   auto-creating the raster. Increase if the C kernel reports that the
#'   grid is not large enough. Default 0.5.
#' @param dim_size Integer. Number of cells along the longest dimension
#'   of the auto-generated raster. Higher values give finer resolution
#'   but slower computation. Default 100.
#' @param time_step Numeric or `NULL`. Time step for the Brownian bridge
#'   integration, in minutes. Defaults to 1/15 of the minimum time lag.
#'   Smaller values give more precise UDs but are slower.
#' @param verbose Logical. If `TRUE`, print a computational size estimate.
#' @param ... Additional arguments passed to methods.
#'
#' @return For single-track input: a `terra::SpatRaster` where cell
#'   values sum to 1.0, representing the utilisation distribution.
#'
#'   For multi-track input: a multi-layer `terra::SpatRaster` with one
#'   named layer per track, all on a common grid. Each layer sums to 1.0.
#'
#' @details
#' The UD is computed by evaluating the Brownian bridge probability
#' density at regular time steps along each segment and accumulating
#' the density onto a raster grid. The computation is implemented in C
#' with OpenMP parallelisation of the inner grid loops for performance.
#'
#' For multi-track input, a common raster grid is computed from the
#' combined spatial extent of all tracks, ensuring that UDs are
#' directly comparable.
#'
#' @references
#' Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K.
#' (2012). A dynamic Brownian bridge movement model to estimate utilization
#' distributions for heterogeneous animal movement. *Journal of Animal
#' Ecology*, 81(4), 738-746. \doi{10.1111/j.1365-2656.2012.01955.x}
#'
#' @seealso [mt_dbbmm_variance()] to compute variance separately,
#'   [mt_dbgb_ud()] for the directional variant, [mt_motion_variance()]
#'   to extract variances.
#'
#' @examples
#' \donttest{
#' library(move2)
#' library(sf)
#'
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#'
#' # One-step
#' ud <- mt_dbbmm_ud(f1_proj, location_error = 25,
#'                    window_size = 31, margin = 11)
#' terra::plot(ud)
#'
#' # Two-step (re-use variance for different raster settings)
#' var_obj <- mt_dbbmm_variance(f1_proj, location_error = 25,
#'                               window_size = 31, margin = 11)
#' ud_fine <- mt_dbbmm_ud(var_obj, location_error = 25, dim_size = 500)
#' }
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
                         verbose = TRUE,
                         ...) {
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
                                verbose = TRUE,
                                ...) {
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
                               margin = 11,
                               window_size = 31,
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
                                            margin = 11,
                                            window_size = 31,
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
