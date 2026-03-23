#' Parallel and Orthogonal Decomposition
#'
#' Decompose the displacement between a point and a reference position
#' into parallel and orthogonal components relative to a direction vector.
#'
#' @param mu Matrix of expected positions (n x 2).
#' @param direction_point Matrix or numeric vector defining the direction
#'   (typically the next location).
#' @param point Matrix of actual positions (n x 2).
#'
#' @return A matrix with columns `deltaPara` and `deltaOrth`.
#' @keywords internal
delta_para_orth <- function(mu, direction_point, point) {
  if (!is.matrix(mu)) mu <- matrix(mu, ncol = 2, nrow = nrow(point), byrow = TRUE)
  if (!is.matrix(direction_point)) {
    direction_point <- matrix(direction_point, ncol = 2, nrow = nrow(point), byrow = TRUE)
  }

  C <- sqrt(rowSums((mu - direction_point)^2))
  B <- sqrt(rowSums((point - direction_point)^2))
  A <- sqrt(rowSums((point - mu)^2))

  tmp <- ((A^2 + C^2 - B^2) / (2 * A * C))

  same_loc <- (rowSums(mu == point) == ncol(mu))
  same_loc_dir <- (rowSums(mu == direction_point) == ncol(mu))

  tmp[same_loc_dir] <- sqrt(0.5)
  if (any(same_loc_dir)) {
    warning("Brownian motion assumed, because no direction could be calculated")
  }
  tmp[same_loc] <- 0
  tmp[tmp > 1] <- 1
  tmp[tmp < -1] <- -1

  theta <- acos(tmp)
  cbind(deltaPara = A * cos(theta), deltaOrth = A * sin(theta))
}

#' BGB Variance Estimation (non-dynamic)
#'
#' Estimate bivariate Gaussian bridge variance for a piece of track.
#'
#' @param x_coords Numeric vector of x coordinates.
#' @param y_coords Numeric vector of y coordinates.
#' @param time_mins Numeric vector of timestamps in minutes.
#' @param loc_err Numeric vector of location errors.
#' @param sd_para Initial parallel SD (for optimization).
#' @param sd_orth Initial orthogonal SD (for optimization).
#'
#' @return A list with `sd_para`, `sd_orth`, and `cll` (conditional log-likelihood).
#' @keywords internal
bgb_var_single <- function(x_coords, y_coords, time_mins, loc_err,
                           sd_para = NULL, sd_orth = NULL) {
  n <- length(x_coords)
  if ((n %% 2) != 1) stop("Need an odd number of locations for BGB variance")

  is <- (1:n)[(1:n) %% 2 == 0]
  alphas <- (time_mins[is] - time_mins[is - 1]) /
    (time_mins[is + 1] - time_mins[is - 1])

  coords <- cbind(x_coords, y_coords)
  mus <- coords[is - 1, , drop = FALSE] +
    alphas * (coords[is + 1, , drop = FALSE] - coords[is - 1, , drop = FALSE])

  para_orth <- delta_para_orth(mus, coords[is + 1, , drop = FALSE],
                                coords[is, , drop = FALSE])

  if (is.null(sd_para) || is.null(sd_orth)) {
    # Optimize both
    opt <- optim(c(1, 1), function(pars) {
      errs <- alphas^2 * loc_err[is + 1]^2 +
        (1 - alphas)^2 * loc_err[is - 1]^2
      sd_mul <- alphas * (1 - alphas) * (time_mins[is + 1] - time_mins[is - 1])
      sp <- sqrt(errs + pars[1]^2 * sd_mul)
      so <- sqrt(errs + pars[2]^2 * sd_mul)
      -.Call("llBGBvar", cbind(sp, so)^2, para_orth)
    }, method = "L-BFGS-B", lower = 0, upper = 1e10)

    sd_para <- opt$par[1]
    sd_orth <- opt$par[2]
    cll <- -opt$value
  } else {
    errs <- alphas^2 * loc_err[is + 1]^2 +
      (1 - alphas)^2 * loc_err[is - 1]^2
    sd_mul <- alphas * (1 - alphas) * (time_mins[is + 1] - time_mins[is - 1])
    sp <- sqrt(errs + sd_para^2 * sd_mul)
    so <- sqrt(errs + sd_orth^2 * sd_mul)
    cll <- .Call("llBGBvar", cbind(sp, so)^2, para_orth)
  }

  list(sd_para = sd_para, sd_orth = sd_orth, cll = cll)
}


#' BGB Variance with Breakpoint Detection
#'
#' Test all possible breakpoint locations within a window for independent
#' breaks in parallel and orthogonal variance.
#'
#' @keywords internal
bgb_var_break <- function(x_coords, y_coords, time_mins, loc_err, margin) {
  n <- length(x_coords)
  coords <- cbind(x_coords, y_coords)

  is <- (1:n)[(1:n) %% 2 == 0]
  alphas <- (time_mins[is] - time_mins[is - 1]) /
    (time_mins[is + 1] - time_mins[is - 1])
  mus <- coords[is - 1, , drop = FALSE] +
    alphas * (coords[is + 1, , drop = FALSE] - coords[is - 1, , drop = FALSE])
  para_orth <- delta_para_orth(mus, coords[is + 1, , drop = FALSE],
                                coords[is, , drop = FALSE])

  errs <- alphas^2 * loc_err[is + 1]^2 + (1 - alphas)^2 * loc_err[is - 1]^2
  sd_mul <- alphas * (1 - alphas) * (time_mins[is + 1] - time_mins[is - 1])

  potential_breaks <- 2:(n - 1)
  margin_breaks <- potential_breaks[potential_breaks >= margin &
                                      potential_breaks <= (1 + n - margin) &
                                      (potential_breaks %% 2) == 1]

  para_breaks <- c(NA, margin_breaks)
  orth_breaks <- c(NA, margin_breaks)
  grid <- expand.grid(para_break = para_breaks, orth_break = orth_breaks)
  grid$res <- NA_real_

  optims <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    pb <- grid$para_break[i]
    ob <- grid$orth_break[i]

    init <- c(paraBefore = 100, orthBefore = 100)
    if (!is.na(pb)) init <- c(init, paraAfter = 100)
    if (!is.na(ob)) init <- c(init, orthAfter = 100)

    obj_fn <- function(pars) {
      para_sd <- rep(pars["paraBefore"], length(errs))
      orth_sd <- rep(pars["orthBefore"], length(errs))
      if (!is.na(pb)) {
        para_sd[seq_len(length(errs)) > floor(pb / 2)] <- pars["paraAfter"]
      }
      if (!is.na(ob)) {
        orth_sd[seq_len(length(errs)) > floor(ob / 2)] <- pars["orthAfter"]
      }
      sp <- sqrt(errs + para_sd^2 * sd_mul)
      so <- sqrt(errs + orth_sd^2 * sd_mul)
      -.Call("llBGBvar", cbind(sp, so)^2, para_orth)
    }

    opt <- optim(init, obj_fn, method = "L-BFGS-B",
                 lower = 0, upper = 1e10,
                 control = list(fnscale = 1))
    grid$res[i] <- -opt$value
    optims[[i]] <- opt
  }

  n_params <- 2 + rowSums(!is.na(grid[, c("para_break", "orth_break")]))
  grid$BIC <- -2 * grid$res + n_params * log(n)

  best_idx <- which.min(grid$BIC)
  opt <- optims[[best_idx]]

  if (any(opt$par == 0)) warning("Optimized to zero")

  # Build result vectors
  result <- cbind(
    paraSd = rep(opt$par["paraBefore"], n),
    orthSd = rep(opt$par["orthBefore"], n)
  )

  brk <- grid[best_idx, ]
  if (!is.na(brk$para_break)) {
    result[seq_len(nrow(result)) >= brk$para_break, "paraSd"] <- opt$par["paraAfter"]
  }
  if (!is.na(brk$orth_break)) {
    result[seq_len(nrow(result)) >= brk$orth_break, "orthSd"] <- opt$par["orthAfter"]
  }

  # NA out margins
  min_brk <- min(c(para_breaks, orth_breaks), na.rm = TRUE)
  max_brk <- max(c(para_breaks, orth_breaks), na.rm = TRUE)
  result[seq_len(nrow(result)) < min_brk, ] <- NA
  result[seq_len(nrow(result)) >= max_brk, ] <- NA

  as.data.frame(result)
}
