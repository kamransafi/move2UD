#' @importFrom sf st_coordinates st_crs st_is_longlat st_bbox st_as_sf st_is_empty
#' @importFrom move2 mt_time mt_track_id mt_n_tracks mt_time_lags mt_is_move2
#' @importFrom terra rast ext xFromCol yFromRow ncol nrow ncell values
#' @importFrom units set_units as_units
#' @importFrom stats aggregate optim optimize
#' @useDynLib move2UD, .registration = TRUE
"_PACKAGE"

#' Extract coordinates, timestamps, and time lags from a move2 object
#'
#' Internal helper that validates and extracts the numeric vectors
#' needed by the dBBMM/dBGB algorithms. Performs all input validation.
#'
#' @param x A `move2` object (single track).
#' @return A list with components `x`, `y`, `time_mins`, `n_locs`.
#' @keywords internal
.extract_track_data <- function(x) {
  .validate_move2(x)

  if (mt_n_tracks(x) > 1) {
    stop("Internal error: .extract_track_data called with multiple tracks.",
         call. = FALSE)
  }

  coords <- st_coordinates(x)
  if (any(!is.finite(coords))) {
    stop("Input contains non-finite coordinates (NA, NaN, or Inf).", call. = FALSE)
  }

  ts_numeric <- as.numeric(mt_time(x))
  if (is.unsorted(ts_numeric, strictly = FALSE)) {
    stop("Timestamps are not in chronological order.", call. = FALSE)
  }
  if (any(diff(ts_numeric) == 0)) {
    stop("Input contains duplicate timestamps. Remove duplicates first: ",
         "move2::mt_filter_unique(x)", call. = FALSE)
  }

  time_mins <- ts_numeric / 60

  list(
    x = coords[, 1],
    y = coords[, 2],
    time_mins = time_mins,
    n_locs = nrow(coords)
  )
}

#' Expand location error to match number of locations
#' @keywords internal
.expand_loc_error <- function(location_error, n_locs) {
  if (length(location_error) == 1) {
    location_error <- rep(location_error, n_locs)
  }
  if (length(location_error) != n_locs) {
    stop("location_error must be length 1 or equal to the number of locations.",
         call. = FALSE)
  }
  if (any(is.na(location_error))) {
    stop("location_error must not contain NAs.", call. = FALSE)
  }
  if (any(location_error <= 0)) {
    stop("location_error must be positive.", call. = FALSE)
  }
  location_error
}

#' Validate a move2 object (common checks that apply regardless of track count)
#' @keywords internal
.validate_move2 <- function(x) {
  if (!mt_is_move2(x)) {
    stop("Input must be a move2 object.", call. = FALSE)
  }
  if (st_is_longlat(x)) {
    stop("Cannot use longitude/latitude coordinates. ",
         "Transform to a projected CRS first, e.g.: ",
         "sf::st_transform(x, move2::mt_aeqd_crs(x))",
         call. = FALSE)
  }
  empties <- st_is_empty(x)
  if (any(empties)) {
    stop("Input contains ", sum(empties), " empty geometries (failed GPS fixes). ",
         "Remove them first: x <- x[!sf::st_is_empty(x), ]",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Split a multi-track move2 object into a named list of single-track objects
#' @keywords internal
.split_tracks <- function(x) {
  .validate_move2(x)
  ids <- unique(mt_track_id(x))
  tracks <- lapply(ids, function(id) {
    x[mt_track_id(x) == id, ]
  })
  names(tracks) <- as.character(ids)

  # Validate each track individually
  for (nm in names(tracks)) {
    trk <- tracks[[nm]]
    if (nrow(trk) == 0) next
    coords <- st_coordinates(trk)
    if (any(!is.finite(coords))) {
      stop("Track '", nm, "' contains non-finite coordinates.", call. = FALSE)
    }
    ts_numeric <- as.numeric(mt_time(trk))
    if (is.unsorted(ts_numeric, strictly = FALSE)) {
      stop("Track '", nm, "' has timestamps not in chronological order.", call. = FALSE)
    }
    if (any(diff(ts_numeric) == 0)) {
      stop("Track '", nm, "' contains duplicate timestamps. ",
           "Remove duplicates first: move2::mt_filter_unique(x)", call. = FALSE)
    }
  }

  tracks
}

#' Calculate extended bounding box for raster creation
#' @keywords internal
.extcalc <- function(x, ext = 0.3) {
  bb <- st_bbox(x)
  x_range <- as.numeric(bb["xmax"] - bb["xmin"])
  y_range <- as.numeric(bb["ymax"] - bb["ymin"])
  c(
    xmin = as.numeric(bb["xmin"]) - ext * x_range,
    xmax = as.numeric(bb["xmax"]) + ext * x_range,
    ymin = as.numeric(bb["ymin"]) - ext * y_range,
    ymax = as.numeric(bb["ymax"]) + ext * y_range
  )
}

#' Create a raster grid for UD computation
#'
#' @param x A move2 or sf object to derive the extent from.
#' @param cell_size Numeric cell size in map units, or NULL to auto-compute.
#' @param dim_size Number of cells along the longest dimension (used if cell_size is NULL).
#' @param ext Extension factor for the bounding box.
#' @return A `terra::SpatRaster` with the appropriate extent, resolution, and CRS.
#' @keywords internal
.make_raster <- function(x, cell_size = NULL, dim_size = 10, ext = 0.3) {
  range <- .extcalc(x, ext = ext)
  x_range <- range["xmax"] - range["xmin"]
  y_range <- range["ymax"] - range["ymin"]

  if (is.null(cell_size)) {
    cell_size <- max(x_range, y_range) / dim_size
  }

  ymin <- range["ymin"] - (ceiling(y_range / cell_size) * cell_size - y_range) / 2
  ymax <- range["ymax"] + (ceiling(y_range / cell_size) * cell_size - y_range) / 2
  xmin <- range["xmin"] - (ceiling(x_range / cell_size) * cell_size - x_range) / 2
  xmax <- range["xmax"] + (ceiling(x_range / cell_size) * cell_size - x_range) / 2

  nr <- round((ymax - ymin) / cell_size)
  nc <- round((xmax - xmin) / cell_size)

  r <- rast(
    nrows = nr, ncols = nc,
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    crs = st_crs(x)$wkt
  )
  r
}
