#' Dynamic Bivariate Gaussian Bridge — Utilisation Distribution
#'
#' Compute a utilisation distribution using the dynamic bivariate Gaussian
#' bridge model, which decomposes movement variance into parallel and
#' orthogonal components relative to the direction of travel.
#'
#' @param object One of:
#'   \itemize{
#'     \item A `move2` object (single or multi-track, projected CRS).
#'     \item An `mt_dbgb_variance` object from [mt_dbgb_variance()].
#'     \item A named list of `mt_dbgb_variance` objects (multi-track).
#'   }
#' @param raster A `terra::SpatRaster` defining the output grid, a numeric
#'   cell size in map units, or `NULL` to auto-compute.
#' @param location_error Numeric scalar or vector of location errors in
#'   map units. Must be positive.
#' @param margin Integer (odd). Margin for variance estimation. Only used
#'   when `object` is a `move2` object.
#' @param window_size Integer (odd). Window size for variance estimation.
#'   Only used when `object` is a `move2` object.
#' @param ext Numeric. Extension factor for the bounding box. Default 0.5.
#' @param dim_size Integer. Cells along the longest dimension. Default 100.
#' @param time_step Numeric or `NULL`. Integration time step in minutes.
#' @param verbose Logical. Print progress messages for multi-track.
#' @param ... Additional arguments passed to methods.
#'
#' @return For single-track: a `terra::SpatRaster` (values sum to 1.0).
#'   For multi-track: a multi-layer `terra::SpatRaster` on a common grid,
#'   one named layer per track, each summing to 1.0.
#'
#' @details
#' The dBGB UD uses an anisotropic Gaussian kernel at each time step,
#' with the kernel elongated along the direction of travel. This produces
#' narrower UDs along directed segments and wider UDs where movement is
#' more random, better capturing the actual space use of the animal.
#'
#' @references
#' Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian
#' bridges: directional factorization of diffusion in Brownian bridge
#' models. *Movement Ecology*, 2(1), 5. \doi{10.1186/2051-3933-2-5}
#'
#' @seealso [mt_dbgb_variance()] to compute variance separately,
#'   [mt_dbbmm_ud()] for the isotropic variant.
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
#' ud <- mt_dbgb_ud(f1_proj, location_error = 25,
#'                   margin = 15, window_size = 31)
#' terra::plot(ud)
#' }
#'
#' @export
mt_dbgb_ud <- function(object,
                        raster = NULL,
                        location_error,
                        margin = 15,
                        window_size = 31,
                        ext = 0.5,
                        dim_size = 100,
                        time_step = NULL,
                        verbose = TRUE,
                        ...) {
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
                               time_step = NULL,
                               verbose = TRUE,
                               ...) {
  var_obj <- mt_dbgb_variance(object, location_error = location_error,
                               margin = margin, window_size = window_size)
  mt_dbgb_ud(var_obj, raster = raster, location_error = location_error,
              ext = ext, dim_size = dim_size, time_step = time_step,
              verbose = verbose)
}

#' @export
mt_dbgb_ud.list <- function(object,
                              raster = NULL,
                              location_error,
                              margin = 15,
                              window_size = 31,
                              ext = 0.5,
                              dim_size = 100,
                              time_step = NULL,
                              verbose = TRUE,
                              ...) {
  if (!all(vapply(object, inherits, logical(1), "mt_dbgb_variance"))) {
    stop("List must contain mt_dbgb_variance objects.", call. = FALSE)
  }

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

  layers <- lapply(names(object), function(nm) {
    if (verbose) message("Computing UD for track: ", nm)
    mt_dbgb_ud(object[[nm]], raster = common_raster,
                location_error = location_error,
                ext = ext, dim_size = dim_size,
                time_step = time_step, verbose = FALSE)
  })

  stk <- terra::rast(layers)
  names(stk) <- names(object)
  stk
}

#' @export
mt_dbgb_ud.mt_dbgb_variance <- function(object,
                                          raster = NULL,
                                          location_error,
                                          margin = 15,
                                          window_size = 31,
                                          ext = 0.5,
                                          dim_size = 100,
                                          time_step = NULL,
                                          verbose = TRUE,
                                          ...) {
  td <- object$track_data
  n <- td$n_locs
  location_error <- .expand_loc_error(location_error, n)

  # Match original move package
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
