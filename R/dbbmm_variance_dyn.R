#' Dynamic Brownian Motion Variance Estimation
#'
#' Estimate the dynamic Brownian motion variance using a sliding window
#' approach with BIC-based breakpoint detection.
#'
#' @param object A `move2` object (single or multi-track). Must be in a
#'   projected coordinate system. Use [move2::mt_aeqd_crs()] to project.
#' @param location_error Numeric scalar or vector of location errors (in map
#'   units) for each position.
#' @param window_size Integer (must be odd). The number of locations in each
#'   sliding window.
#' @param margin Integer (must be odd). The minimum number of locations on
#'   each side of a potential breakpoint within a window.
#' @param parallel Logical. If `TRUE`, use parallel processing for the
#'   sliding window loop. Defaults to `FALSE`.
#' @param cores Integer. Number of cores for parallel processing.
#'
#' @return For a single-track input: an `mt_dbbmm_variance` object.
#'   For multi-track input: a named list of `mt_dbbmm_variance` objects,
#'   one per track. Tracks with too few locations for the given
#'   window_size are skipped with a warning.
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

#' Helper: run lapply or mclapply depending on parallel flag
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
