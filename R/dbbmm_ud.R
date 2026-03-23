#' Dynamic Brownian Bridge Movement Model — Utilisation Distribution
#'
#' Compute a utilisation distribution (UD) from a movement track using
#' the dynamic Brownian bridge movement model.
#'
#' @param object A `move2` object (single track, projected CRS) or an
#'   `mt_dbbmm_variance` object from [mt_dbbmm_variance()].
#' @param raster A `terra::SpatRaster` defining the output grid, or a
#'   numeric scalar giving the cell size in map units, or `NULL` to
#'   auto-compute from `dim_size`.
#' @param location_error Numeric scalar or vector of location errors.
#' @param margin Integer (odd). Margin for the variance estimation window.
#' @param window_size Integer (odd). Window size for variance estimation.
#' @param ext Numeric. Extension factor for the bounding box when
#'   auto-creating the raster.
#' @param dim_size Integer. Number of cells along the longest dimension
#'   (used when `raster` is NULL).
#' @param time_step Numeric. Time step for integration (in minutes).
#'   Defaults to 1/15 of the minimum time lag.
#' @param verbose Logical. Print computational size estimate.
#'
#' @return A `terra::SpatRaster` containing the utilisation distribution
#'   (values sum to 1).
#'
#' @details
#' If `object` is a `move2` object, the variance is estimated first
#' using [mt_dbbmm_variance()], then the UD is computed. If `object`
#' is already an `mt_dbbmm_variance` object, the UD is computed directly.
#'
#' @references
#' Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K.
#' (2012). A dynamic Brownian bridge movement model to estimate utilization
#' distributions for heterogeneous animal movement. *Journal of Animal
#' Ecology*, 81(4), 738-746.
#'
#' @examples
#' \dontrun{
#' library(move2)
#' library(sf)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#' ud <- mt_dbbmm_ud(f1_proj, location_error = 25,
#'                    window_size = 31, margin = 11)
#' terra::plot(ud)
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
