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
#' using BIC for model selection. The window slides one position at a time
#' across the entire trajectory, and the final variance for each segment is
#' the mean across all windows that covered it.
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
mt_dbbmm_variance <- function(object, location_error, window_size, margin) {
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
  breaks_found <- c()
  bm_vars <- data.frame(BMvar = numeric(0), loc = numeric(0))

  for (w in 1:(n - window_size + 1)) {
    idx <- w:(w - 1 + window_size)
    x_sub <- td$x[idx]
    y_sub <- td$y[idx]
    tl_sub <- time_lag[idx]
    le_sub <- location_error[idx]

    whole <- bm_variance(time_lag = tl_sub, location_error = le_sub,
                         x = x_sub, y = y_sub)
    whole$BIC <- -2 * whole$cll + log(window_size)

    break_best <- list(BIC = Inf)
    for (b in uneven_breaks) {
      before <- bm_variance(
        x = x_sub[1:b], y = y_sub[1:b],
        location_error = le_sub[1:b], time_lag = tl_sub[1:b]
      )
      after <- bm_variance(
        x = x_sub[b:window_size], y = y_sub[b:window_size],
        location_error = le_sub[b:window_size], time_lag = tl_sub[b:window_size]
      )
      bic_break <- -2 * (before$cll + after$cll) + 2 * log(window_size)
      if (bic_break < break_best$BIC) {
        break_best <- list(
          w = w, b = b,
          var_before = before$BMvar, var_after = after$BMvar,
          BIC = bic_break
        )
      }
    }

    if (break_best$BIC < whole$BIC) {
      window_vars <- c(
        rep(break_best$var_before, sum(breaks_range < break_best$b)),
        rep(break_best$var_after, sum(breaks_range > break_best$b))
      )
      breaks_found <- c(breaks_found, (w - 1 + break_best$b))
    } else {
      window_vars <- rep(whole$BMvar, length(breaks_range) - 1)
    }

    bm_vars <- rbind(bm_vars, data.frame(
      BMvar = window_vars,
      loc = w - 1 + margin:(window_size - margin)
    ))
  }

  agg <- aggregate(BMvar ~ loc, data = bm_vars,
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
