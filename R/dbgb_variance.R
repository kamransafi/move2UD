#' Parallel and Orthogonal Decomposition
#'
#' Decompose the displacement between a point and a reference position
#' into parallel and orthogonal components relative to a direction vector.
#'
#' @param mu Matrix of expected positions (n x 2).
#' @param direction_point Matrix or numeric vector defining the direction.
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


#' BGB Variance with Breakpoint Detection (optimised)
#'
#' Test breakpoint locations within a window using a sequential search:
#' first find the best parallel break, then the best orthogonal break
#' conditional on it. This reduces the search from O(n^2) to O(n).
#'
#' @keywords internal
bgb_var_break <- function(x_coords, y_coords, time_mins, location_error, margin) {
  n <- length(x_coords)
  coords <- cbind(x_coords, y_coords)

  is <- (1:n)[(1:n) %% 2 == 0]
  alphas <- (time_mins[is] - time_mins[is - 1]) /
    (time_mins[is + 1] - time_mins[is - 1])
  mus <- coords[is - 1, , drop = FALSE] +
    alphas * (coords[is + 1, , drop = FALSE] - coords[is - 1, , drop = FALSE])
  para_orth <- delta_para_orth(mus, coords[is + 1, , drop = FALSE],
                                coords[is, , drop = FALSE])

  errs <- alphas^2 * location_error[is + 1]^2 + (1 - alphas)^2 * location_error[is - 1]^2
  sd_mul <- alphas * (1 - alphas) * (time_mins[is + 1] - time_mins[is - 1])
  n_pairs <- length(errs)

  potential_breaks <- 2:(n - 1)
  margin_breaks <- potential_breaks[potential_breaks >= margin &
                                      potential_breaks <= (1 + n - margin) &
                                      (potential_breaks %% 2) == 1]

  # Helper: evaluate log-likelihood for given para/orth SDs
  eval_ll <- function(para_sd_vec, orth_sd_vec) {
    sp <- sqrt(errs + para_sd_vec^2 * sd_mul)
    so <- sqrt(errs + orth_sd_vec^2 * sd_mul)
    .Call("llBGBvar", cbind(sp, so)^2, para_orth)
  }

  # Helper: optimize para and orth SDs given break positions
  optimize_sds <- function(para_break, orth_break) {
    init <- c(paraBefore = 100, orthBefore = 100)
    if (!is.na(para_break)) init <- c(init, paraAfter = 100)
    if (!is.na(orth_break)) init <- c(init, orthAfter = 100)

    obj_fn <- function(pars) {
      para_sd <- rep(pars["paraBefore"], n_pairs)
      orth_sd <- rep(pars["orthBefore"], n_pairs)
      if (!is.na(para_break)) {
        para_sd[seq_len(n_pairs) > floor(para_break / 2)] <- pars["paraAfter"]
      }
      if (!is.na(orth_break)) {
        orth_sd[seq_len(n_pairs) > floor(orth_break / 2)] <- pars["orthAfter"]
      }
      -eval_ll(para_sd, orth_sd)
    }

    opt <- optim(init, obj_fn, method = "L-BFGS-B",
                 lower = 0, upper = 1e10,
                 control = list(fnscale = 1))
    n_params <- 2 + (!is.na(para_break)) + (!is.na(orth_break))
    bic <- -2 * (-opt$value) + n_params * log(n)
    list(opt = opt, bic = bic)
  }

  # Step 1: No break at all
  res_none <- optimize_sds(NA, NA)
  best <- res_none
  best_pb <- NA
  best_ob <- NA

  # Step 2: Try each para break (no orth break)
  for (pb in margin_breaks) {
    res <- optimize_sds(pb, NA)
    if (res$bic < best$bic) {
      best <- res
      best_pb <- pb
      best_ob <- NA
    }
  }

  # Step 3: Try each orth break (no para break)
  for (ob in margin_breaks) {
    res <- optimize_sds(NA, ob)
    if (res$bic < best$bic) {
      best <- res
      best_pb <- NA
      best_ob <- ob
    }
  }

  # Step 4: Try orth breaks conditional on the best para break, and vice versa
  if (!is.na(best_pb)) {
    for (ob in margin_breaks) {
      res <- optimize_sds(best_pb, ob)
      if (res$bic < best$bic) {
        best <- res
        best_ob <- ob
      }
    }
  }
  if (!is.na(best_ob)) {
    for (pb in margin_breaks) {
      res <- optimize_sds(pb, best_ob)
      if (res$bic < best$bic) {
        best <- res
        best_pb <- pb
      }
    }
  }

  opt <- best$opt
  if (any(opt$par == 0)) warning("Optimized to zero")

  # Build result vectors
  result <- cbind(
    paraSd = rep(opt$par["paraBefore"], n),
    orthSd = rep(opt$par["orthBefore"], n)
  )

  if (!is.na(best_pb)) {
    result[seq_len(nrow(result)) >= best_pb, "paraSd"] <- opt$par["paraAfter"]
  }
  if (!is.na(best_ob)) {
    result[seq_len(nrow(result)) >= best_ob, "orthSd"] <- opt$par["orthAfter"]
  }

  # NA out margins
  min_brk <- min(c(margin_breaks), na.rm = TRUE)
  max_brk <- max(c(margin_breaks), na.rm = TRUE)
  result[seq_len(nrow(result)) < min_brk, ] <- NA
  result[seq_len(nrow(result)) >= max_brk, ] <- NA

  as.data.frame(result)
}
