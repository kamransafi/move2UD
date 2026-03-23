#' Dynamic Brownian Motion Variance Estimation
#'
#' Estimate the dynamic Brownian motion variance using a sliding window
#' approach with BIC-based breakpoint detection. Accepts single-track or
#' multi-track `move2` objects.
#'
#' @param object A `move2` object in a projected coordinate system (not
#'   longitude/latitude). Use `sf::st_transform(x, move2::mt_aeqd_crs(x))`
#'   to project. Both single-track and multi-track objects are accepted.
#'   Empty geometries must be removed beforehand.
#' @param location_error Numeric scalar or vector of location errors in map
#'   units (typically metres). If a scalar, the same error is used for all
#'   locations. Must be positive.
#' @param window_size Integer (must be odd). The number of locations in each
#'   sliding window. Larger values produce smoother variance estimates but
#'   miss short-term behavioural changes.
#' @param margin Integer (must be odd). The minimum number of locations on
#'   each side of a potential breakpoint within a window.
#'   Must satisfy `window_size >= 2 * margin + 1`.
#' @param parallel Logical. If `TRUE`, use parallel processing for the
#'   sliding window loop via [parallel::mclapply()]. Only effective on
#'   Unix/macOS; falls back to sequential on Windows with a message.
#'   Defaults to `FALSE`.
#' @param cores Integer or `NULL`. Number of cores for parallel processing.
#'   Defaults to `parallel::detectCores() - 1`.
#'
#' @return For a single-track input: an `mt_dbbmm_variance` S3 object
#'   containing:
#'   \describe{
#'     \item{`variance`}{Numeric vector of estimated BM variances per
#'       location (`NA` for positions outside the estimable range at the
#'       track margins).}
#'     \item{`in_windows`}{Numeric vector: number of overlapping windows
#'       each location was estimated in.}
#'     \item{`interest`}{Logical vector: `TRUE` for locations fully
#'       covered by the sliding window (maximum overlap).}
#'     \item{`break_list`}{Integer vector of positions where behavioural
#'       breakpoints were detected.}
#'     \item{`window_size`, `margin`}{The parameters used.}
#'     \item{`track_data`}{List with `x`, `y`, `time_mins`, `n_locs`,
#'       `crs`, `timestamps` from the input.}
#'   }
#'
#'   For multi-track input: a named list of `mt_dbbmm_variance` objects,
#'   one per track. Tracks with fewer locations than `window_size` are
#'   skipped with a warning.
#'
#' @details
#' The method estimates Brownian motion variance using a leave-one-out
#' likelihood approach within a sliding window (Horne et al. 2007). At
#' each window position, it tests whether a single variance or two
#' variances (split at a breakpoint) better fit the data, using BIC for
#' model selection (Kranstauber et al. 2012). The window slides one
#' position at a time, and the final variance for each location is the
#' mean across all windows that covered it.
#'
#' The variance estimation and breakpoint testing are implemented in C
#' (Brent's method optimizer) for performance. The sliding window loop
#' can optionally be parallelised across CPU cores.
#'
#' @references
#' Horne, J. S., Garton, E. O., Krone, S. M., & Lewis, J. S. (2007).
#' Analyzing animal movements using Brownian bridges. *Ecology*, 88(9),
#' 2354-2363. \doi{10.1890/06-0957.1}
#'
#' Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K.
#' (2012). A dynamic Brownian bridge movement model to estimate utilization
#' distributions for heterogeneous animal movement. *Journal of Animal
#' Ecology*, 81(4), 738-746. \doi{10.1111/j.1365-2656.2012.01955.x}
#'
#' @seealso [mt_dbbmm_ud()] to compute the utilisation distribution from
#'   the variance estimate, [mt_motion_variance()] to extract the variance
#'   vector, [mt_dbgb_variance()] for the directional (bivariate) variant.
#'
#' @examples
#' \donttest{
#' library(move2)
#' library(sf)
#'
#' # Load and prepare example data
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#'
#' # Single track
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#' var_obj <- mt_dbbmm_variance(f1_proj, location_error = 25,
#'                               window_size = 31, margin = 11)
#' var_obj
#' plot(mt_time(f1_proj), mt_motion_variance(var_obj),
#'      type = "l", xlab = "Time", ylab = "BM variance")
#'
#' # Multiple tracks
#' fishers_proj <- st_transform(fishers, mt_aeqd_crs(fishers))
#' var_list <- mt_dbbmm_variance(fishers_proj, location_error = 25,
#'                                window_size = 31, margin = 11)
#' names(var_list)
#' }
#'
#' @export
mt_dbbmm_variance <- function(object, location_error, window_size, margin,
                               parallel = FALSE, cores = NULL) {
  .validate_move2(object)

  if (mt_n_tracks(object) > 1) {
    return(.dbbmm_variance_multi(object, location_error, window_size, margin,
                                  parallel, cores))
  }

  .dbbmm_variance_single(object, location_error, window_size, margin,
                           parallel, cores)
}

#' @keywords internal
.dbbmm_variance_single <- function(object, location_error, window_size, margin,
                                    parallel = FALSE, cores = NULL) {
  td <- .extract_track_data(object)
  n <- td$n_locs

  time_lag <- c(diff(td$time_mins), 0)
  location_error <- .expand_loc_error(location_error, n)

  if (n < window_size) {
    stop("window_size (", window_size, ") cannot be larger than the number of ",
         "locations (", n, ").", call. = FALSE)
  }
  if (any((c(margin, window_size) %% 2) != 1)) {
    stop("margin and window_size must both be odd.", call. = FALSE)
  }
  if (window_size < 2 * margin + 1) {
    stop("window_size must be at least 2 * margin + 1.", call. = FALSE)
  }

  breaks_range <- margin:(window_size - margin + 1)
  uneven_breaks <- breaks_range[(breaks_range %% 2) == 1]
  n_windows <- n - window_size + 1

  process_window <- function(w) {
    idx <- w:(w - 1 + window_size)
    res <- .Call("bm_variance_window_c",
                 td$x[idx], td$y[idx], time_lag[idx], location_error[idx],
                 as.integer(uneven_breaks), as.integer(breaks_range))
    list(
      vars = res$vars,
      locs = w - 1 + margin:(window_size - margin),
      break_pos = if (res$has_break) w - 1 + res$break_pos else NULL
    )
  }

  results <- .run_lapply(1:n_windows, process_window, parallel, cores)

  all_vars <- do.call(rbind, lapply(results, function(r) {
    data.frame(BMvar = r$vars, loc = r$locs)
  }))
  breaks_found <- unlist(lapply(results, function(r) r$break_pos))

  agg <- aggregate(BMvar ~ loc, data = all_vars,
                   FUN = function(x) c(mean = mean(x), length = length(x)))
  means_mat <- agg$BMvar

  variance <- rep(NA_real_, n)
  in_windows <- rep(NA_real_, n)
  variance[agg$loc] <- means_mat[, "mean"]
  in_windows[agg$loc] <- means_mat[, "length"]

  interest <- rep(FALSE, n)
  interest[agg$loc] <- means_mat[, "length"] == max(means_mat[, "length"])

  if (is.null(breaks_found)) breaks_found <- integer(0)

  structure(
    list(
      variance = variance,
      in_windows = in_windows,
      interest = interest,
      break_list = breaks_found,
      window_size = window_size,
      margin = margin,
      track_data = list(
        x = td$x, y = td$y,
        time_mins = td$time_mins,
        n_locs = n,
        crs = st_crs(object),
        timestamps = mt_time(object)
      )
    ),
    class = "mt_dbbmm_variance"
  )
}

#' @keywords internal
.dbbmm_variance_multi <- function(object, location_error, window_size, margin,
                                   parallel, cores) {
  tracks <- .split_tracks(object)
  results <- list()
  skipped <- character(0)

  for (nm in names(tracks)) {
    trk <- tracks[[nm]]
    n <- nrow(trk)
    if (n < window_size) {
      skipped <- c(skipped, nm)
      next
    }
    results[[nm]] <- .dbbmm_variance_single(
      trk, location_error = location_error,
      window_size = window_size, margin = margin,
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

#' @keywords internal
.run_lapply <- function(X, FUN, parallel, cores) {
  if (parallel) {
    if (.Platform$OS.type == "windows") {
      message("Note: parallel processing uses mclapply which is not available ",
              "on Windows. Falling back to sequential processing.")
      lapply(X, FUN)
    } else {
      if (is.null(cores)) cores <- max(1, parallel::detectCores() - 1)
      parallel::mclapply(X, FUN, mc.cores = cores)
    }
  } else {
    lapply(X, FUN)
  }
}


#' Print method for mt_dbbmm_variance
#'
#' @param x An `mt_dbbmm_variance` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mt_dbbmm_variance <- function(x, ...) {
  cat("Dynamic Brownian Bridge Movement Model \u2014 variance estimate\n")
  cat(sprintf("  Locations: %d\n", x$track_data$n_locs))
  cat(sprintf("  Window size: %d, Margin: %d\n", x$window_size, x$margin))
  cat(sprintf("  Breakpoints found: %d\n", length(x$break_list)))
  cat(sprintf("  Variance range: %.2f \u2013 %.2f\n",
              min(x$variance, na.rm = TRUE), max(x$variance, na.rm = TRUE)))
  invisible(x)
}


#' Extract Motion Variance
#'
#' Extract the estimated movement variance from a dBBMM or dBGB variance
#' object.
#'
#' @param x An `mt_dbbmm_variance` or `mt_dbgb_variance` object, as
#'   returned by [mt_dbbmm_variance()] or [mt_dbgb_variance()].
#' @param ... Ignored.
#'
#' @return For `mt_dbbmm_variance`: a numeric vector of variances (one
#'   per location, `NA` at track margins).
#'
#'   For `mt_dbgb_variance`: a `data.frame` with columns `para`
#'   (parallel variance) and `orth` (orthogonal variance).
#'
#' @seealso [mt_dbbmm_variance()], [mt_dbgb_variance()]
#'
#' @examples
#' \donttest{
#' library(move2)
#' library(sf)
#' fishers <- mt_read(mt_example())
#' fishers <- fishers[!st_is_empty(fishers), ]
#' f1 <- fishers[mt_track_id(fishers) == "F1", ]
#' f1_proj <- st_transform(f1, mt_aeqd_crs(f1))
#' var_obj <- mt_dbbmm_variance(f1_proj, location_error = 25,
#'                               window_size = 31, margin = 11)
#' head(mt_motion_variance(var_obj))
#' }
#'
#' @export
mt_motion_variance <- function(x, ...) {
  UseMethod("mt_motion_variance")
}

#' @rdname mt_motion_variance
#' @export
mt_motion_variance.mt_dbbmm_variance <- function(x, ...) {
  x$variance
}
