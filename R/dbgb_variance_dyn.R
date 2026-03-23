#' Dynamic Bivariate Gaussian Bridge Variance Estimation
#'
#' Estimate dynamic parallel and orthogonal movement variances using a
#' sliding window approach with BIC-based breakpoint detection.
#'
#' @param object A `move2` object containing a single track. Must be in a
#'   projected coordinate system. Use [move2::mt_aeqd_crs()] to create
#'   a suitable projection.
#' @param location_error Numeric scalar or vector of location errors.
#' @param margin Integer (odd). Margin for breakpoint detection.
#' @param window_size Integer (odd). Sliding window size.
#'
#' @return An `mt_dbgb_variance` object (S3 list) containing:
#'   \describe{
#'     \item{para_sd}{Numeric vector of parallel SD per position.}
#'     \item{orth_sd}{Numeric vector of orthogonal SD per position.}
#'     \item{n_estim}{Number of windows each position was estimated in.}
#'     \item{seg_interest}{Logical vector of segments fully within the window.}
#'     \item{margin}{Margin used.}
#'     \item{window_size}{Window size used.}
#'     \item{track_data}{Coordinates, timestamps, CRS from the input.}
#'   }
#'
#' @references
#' Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian
#' bridges: directional factorization of diffusion in Brownian bridge
#' models. *Movement Ecology*, 2(1), 5.
#'
#' @export
mt_dbgb_variance <- function(object, location_error, margin, window_size) {
  td <- .extract_track_data(object)
  n <- td$n_locs

  location_error <- .expand_loc_error(location_error, n)

  if (n < window_size) {
    stop("window_size cannot be larger than the number of locations.", call. = FALSE)
  }
  if (any((c(margin, window_size) %% 2) != 1)) {
    stop("margin and window_size must both be odd.", call. = FALSE)
  }

  all_results <- vector("list", n - window_size + 1)

  for (w in seq_len(n - window_size + 1)) {
    idx <- w:(w + window_size - 1)
    res <- bgb_var_break(
      x_coords = td$x[idx],
      y_coords = td$y[idx],
      time_mins = td$time_mins[idx],
      location_error = location_error[idx],
      margin = margin
    )
    res[nrow(res), ] <- NA

    all_results[[w]] <- data.frame(
      seg = idx,
      paraSd = res$paraSd,
      orthSd = res$orthSd
    )
  }

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
