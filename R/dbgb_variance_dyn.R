#' Dynamic Bivariate Gaussian Bridge Variance Estimation
#'
#' Estimate dynamic parallel and orthogonal movement variances using a
#' sliding window approach with BIC-based breakpoint detection. The
#' bivariate Gaussian bridge decomposes the movement variance into a
#' component parallel to the direction of travel and one orthogonal to
#' it, providing more detail about the movement process than the
#' isotropic dBBMM.
#'
#' @param object A `move2` object in a projected coordinate system (not
#'   longitude/latitude). Use `sf::st_transform(x, move2::mt_aeqd_crs(x))`
#'   to project. Both single-track and multi-track objects are accepted.
#' @param location_error Numeric scalar or vector of location errors in
#'   map units. Must be positive.
#' @param margin Integer (must be odd). Minimum locations on each side of
#'   a potential breakpoint.
#' @param window_size Integer (must be odd). Number of locations in the
#'   sliding window.
#' @param parallel Logical. If `TRUE`, use parallel processing for the
#'   sliding window. Defaults to `FALSE`.
#' @param cores Integer or `NULL`. Number of cores for parallel processing.
#'
#' @return For single-track input: an `mt_dbgb_variance` S3 object
#'   containing:
#'   \describe{
#'     \item{`para_sd`}{Numeric vector of parallel standard deviations.}
#'     \item{`orth_sd`}{Numeric vector of orthogonal standard deviations.}
#'     \item{`n_estim`}{Number of windows each position was estimated in.}
#'     \item{`seg_interest`}{Logical vector of fully-covered segments.}
#'     \item{`margin`, `window_size`}{Parameters used.}
#'     \item{`track_data`}{Coordinates, timestamps, CRS from the input.}
#'   }
#'
#'   For multi-track input: a named list of `mt_dbgb_variance` objects.
#'   Tracks with too few locations are skipped with a warning.
#'
#' @details
#' The method extends the dBBMM by decomposing the variance into a
#' parallel component (along the direction of travel between consecutive
#' locations) and an orthogonal component (perpendicular to it). This
#' allows distinguishing directed movement (high parallel, low orthogonal
#' variance) from random movement (similar variances in both directions).
#'
#' A directionality index can be computed from the variances:
#' `I_d = (para - orth) / (para + orth)`, where values near 0 indicate
#' Brownian motion and positive values indicate directional movement.
#'
#' @references
#' Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian
#' bridges: directional factorization of diffusion in Brownian bridge
#' models. *Movement Ecology*, 2(1), 5. \doi{10.1186/2051-3933-2-5}
#'
#' @seealso [mt_dbgb_ud()] to compute the UD, [mt_motion_variance()] to
#'   extract variances as a data.frame, [mt_dbbmm_variance()] for the
#'   isotropic variant.
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
#' var_obj <- mt_dbgb_variance(f1_proj, location_error = 25,
#'                              margin = 15, window_size = 31)
#' var_obj
#'
#' # Directionality index
#' mv <- mt_motion_variance(var_obj)
#' I_d <- (mv$para - mv$orth) / (mv$para + mv$orth)
#' plot(mt_time(f1_proj), I_d, type = "l",
#'      xlab = "Time", ylab = "Directionality index")
#' }
#'
#' @export
mt_dbgb_variance <- function(object, location_error, margin, window_size,
                              parallel = FALSE, cores = NULL) {
  .validate_move2(object)

  if (mt_n_tracks(object) > 1) {
    return(.dbgb_variance_multi(object, location_error, margin, window_size,
                                 parallel, cores))
  }

  .dbgb_variance_single(object, location_error, margin, window_size,
                          parallel, cores)
}

#' @keywords internal
.dbgb_variance_single <- function(object, location_error, margin, window_size,
                                   parallel = FALSE, cores = NULL) {
  td <- .extract_track_data(object)
  n <- td$n_locs

  location_error <- .expand_loc_error(location_error, n)

  if (n < window_size) {
    stop("window_size (", window_size, ") cannot be larger than the number of ",
         "locations (", n, ").", call. = FALSE)
  }
  if (any((c(margin, window_size) %% 2) != 1)) {
    stop("margin and window_size must both be odd.", call. = FALSE)
  }

  n_windows <- n - window_size + 1

  process_window <- function(w) {
    idx <- w:(w + window_size - 1)
    res <- bgb_var_break(
      x_coords = td$x[idx],
      y_coords = td$y[idx],
      time_mins = td$time_mins[idx],
      location_error = location_error[idx],
      margin = margin
    )
    res[nrow(res), ] <- NA
    data.frame(seg = idx, paraSd = res$paraSd, orthSd = res$orthSd)
  }

  all_results <- .run_lapply(seq_len(n_windows), process_window, parallel, cores)

  combined <- do.call(rbind, all_results)
  combined <- combined[!is.na(combined$paraSd), ]

  agg <- aggregate(cbind(paraSd, orthSd) ~ seg, data = combined,
                   FUN = function(x) sqrt(mean(x^2)))
  n_estim_df <- aggregate(paraSd ~ seg, data = combined, FUN = length)

  para_sd <- rep(NA_real_, n)
  orth_sd <- rep(NA_real_, n)
  n_estim <- rep(NA_real_, n)

  para_sd[agg$seg] <- agg$paraSd
  orth_sd[agg$seg] <- agg$orthSd
  n_estim[n_estim_df$seg] <- n_estim_df$paraSd

  seg_interest <- !is.na(n_estim) & (n_estim == max(n_estim, na.rm = TRUE))

  structure(
    list(
      para_sd = para_sd,
      orth_sd = orth_sd,
      n_estim = n_estim,
      seg_interest = seg_interest,
      margin = margin,
      window_size = window_size,
      track_data = list(
        x = td$x, y = td$y,
        time_mins = td$time_mins,
        n_locs = n,
        crs = st_crs(object),
        timestamps = mt_time(object)
      )
    ),
    class = "mt_dbgb_variance"
  )
}

#' @keywords internal
.dbgb_variance_multi <- function(object, location_error, margin, window_size,
                                  parallel, cores) {
  tracks <- .split_tracks(object)
  results <- list()
  skipped <- character(0)

  for (nm in names(tracks)) {
    trk <- tracks[[nm]]
    if (nrow(trk) < window_size) {
      skipped <- c(skipped, nm)
      next
    }
    results[[nm]] <- .dbgb_variance_single(
      trk, location_error = location_error,
      margin = margin, window_size = window_size,
      parallel = parallel, cores = cores
    )
  }

  if (length(skipped) > 0) {
    warning("Tracks skipped (fewer locations than window_size): ",
            paste(skipped, collapse = ", "), call. = FALSE)
  }
  if (length(results) == 0) {
    stop("No tracks had enough locations for the given window_size.", call. = FALSE)
  }

  results
}

#' Print method for mt_dbgb_variance
#'
#' @param x An `mt_dbgb_variance` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mt_dbgb_variance <- function(x, ...) {
  cat("Dynamic Bivariate Gaussian Bridge \u2014 variance estimate\n")
  cat(sprintf("  Locations: %d\n", x$track_data$n_locs))
  cat(sprintf("  Window size: %d, Margin: %d\n", x$window_size, x$margin))
  cat(sprintf("  Parallel SD range: %.2f \u2013 %.2f\n",
              min(x$para_sd, na.rm = TRUE), max(x$para_sd, na.rm = TRUE)))
  cat(sprintf("  Orthogonal SD range: %.2f \u2013 %.2f\n",
              min(x$orth_sd, na.rm = TRUE), max(x$orth_sd, na.rm = TRUE)))
  invisible(x)
}

#' @rdname mt_motion_variance
#' @export
mt_motion_variance.mt_dbgb_variance <- function(x, ...) {
  data.frame(para = x$para_sd^2, orth = x$orth_sd^2)
}
