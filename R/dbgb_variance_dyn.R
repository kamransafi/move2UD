#' Dynamic Bivariate Gaussian Bridge Variance Estimation
#'
#' Estimate dynamic parallel and orthogonal movement variances using a
#' sliding window approach with BIC-based breakpoint detection.
#'
#' @param object A `move2` object (single or multi-track). Must be in a
#'   projected coordinate system.
#' @param location_error Numeric scalar or vector of location errors.
#' @param margin Integer (odd). Margin for breakpoint detection.
#' @param window_size Integer (odd). Sliding window size.
#' @param parallel Logical. If `TRUE`, use parallel processing.
#' @param cores Integer. Number of cores for parallel processing.
#'
#' @return For single-track input: an `mt_dbgb_variance` object.
#'   For multi-track input: a named list of `mt_dbgb_variance` objects.
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

#' @export
mt_motion_variance.mt_dbgb_variance <- function(x, ...) {
  data.frame(para = x$para_sd^2, orth = x$orth_sd^2)
}
