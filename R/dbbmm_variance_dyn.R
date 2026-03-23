#' Dynamic Brownian Motion Variance Estimation
#'
#' Estimate the dynamic Brownian motion variance using a sliding window
#' approach with BIC-based breakpoint detection.
#'
#' @param object A `move2` object containing a single track. Must be in a
#'   projected coordinate system (not longitude/latitude). Use
#'   [move2::mt_aeqd_crs()] to create a suitable projection.
#' @param location_error Numeric scalar or vector of location errors (in map
#'   units) for each position.
#' @param window_size Integer (must be odd). The number of locations in each
#'   sliding window.
#' @param margin Integer (must be odd). The minimum number of locations on
#'   each side of a potential breakpoint within a window.
#' @param parallel Logical. If `TRUE`, use parallel processing for the
#'   sliding window loop. Defaults to `FALSE`. Uses [parallel::mclapply()]
#'   on Unix/macOS and sequential processing on Windows.
#' @param cores Integer. Number of cores for parallel processing. Defaults
#'   to `parallel::detectCores() - 1`.
#'
#' @return An `mt_dbbmm_variance` object (S3 list) containing:
#'   \describe{
#'     \item{variance}{Numeric vector of estimated BM variances per segment
#'       (NA for positions outside the estimable range).}
#'     \item{in_windows}{Number of windows each segment was estimated in.}
#'     \item{interest}{Logical vector indicating segments fully covered by
#'       the sliding window.}
#'     \item{break_list}{Integer vector of detected breakpoint positions.}
#'     \item{window_size}{The window size used.}
#'     \item{margin}{The margin used.}
#'     \item{track_data}{Coordinates, timestamps, and CRS from the input.}
#'   }
#'
#' @details
#' The method estimates BM variance using a leave-one-out likelihood approach
#' within a sliding window. At each window position, it tests whether a single
#' variance or two variances (split at a breakpoint) better fit the data,
#' using BIC for model selection.
#'
#' The variance estimation and breakpoint testing are implemented in C for
#' performance. The sliding window loop can optionally be parallelised across
#' cores.
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
#' var_obj <- mt_dbbmm_variance(f1_proj, location_error = 25,
#'                               window_size = 31, margin = 11)
#' }
#'
#' @export
mt_dbbmm_variance <- function(object, location_error, window_size, margin,
                               parallel = FALSE, cores = NULL) {
  td <- .extract_track_data(object)
  n <- td$n_locs

  time_lag <- c(diff(td$time_mins), 0)
  location_error <- .expand_loc_error(location_error, n)

  if (n < window_size) {
    stop("window_size cannot be larger than the number of locations.", call. = FALSE)
  }
  if (any((c(margin, window_size) %% 2) != 1)) {
    stop("margin and window_size must both be odd.", call. = FALSE)
  }

  breaks_range <- margin:(window_size - margin + 1)
  if (is.unsorted(breaks_range)) {
    stop("window_size must be at least 2 * margin.", call. = FALSE)
  }
  uneven_breaks <- breaks_range[(breaks_range %% 2) == 1]

  n_windows <- n - window_size + 1

  # Function to process one window position using C
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

  # Process all windows — parallel or sequential
  if (parallel && .Platform$OS.type != "windows") {
    if (is.null(cores)) cores <- max(1, parallel::detectCores() - 1)
    results <- parallel::mclapply(1:n_windows, process_window, mc.cores = cores)
  } else {
    results <- lapply(1:n_windows, process_window)
  }

  # Aggregate results
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

#' Extract motion variance
#'
#' @param x An `mt_dbbmm_variance` or `mt_dbgb_variance` object.
#' @param ... Ignored.
#' @return For `mt_dbbmm_variance`: a numeric vector of variances.
#'   For `mt_dbgb_variance`: a data.frame with columns `para` and `orth`.
#' @export
mt_motion_variance <- function(x, ...) {
  UseMethod("mt_motion_variance")
}

#' @export
mt_motion_variance.mt_dbbmm_variance <- function(x, ...) {
  x$variance
}
