# move2UD

**Dynamic Brownian Bridge Movement Models for `move2` Objects**

`move2UD` provides utilisation distribution (UD) estimation from animal tracking data using dynamic Brownian bridge movement models (dBBMM) and dynamic bivariate Gaussian bridges (dBGB). It is designed to work with [`move2`](https://bartk.gitlab.io/move2/), [`sf`](https://r-spatial.github.io/sf/), and [`terra`](https://rspatial.github.io/terra/) — the modern R spatial stack.

This package migrates the dBBMM/dBGB functionality previously available in the [`move`](https://cran.r-project.org/package=move) package (which depends on the retired `sp`/`raster` ecosystem) into the current generation of spatial tools. The C computational kernels are preserved from the original implementation; the R interface layer has been rewritten.

## Installation

```r
# Install from GitHub
devtools::install_github("kamransafi/move2UD")
```

### Dependencies

All dependencies are available on CRAN:

- [`move2`](https://cran.r-project.org/package=move2) (>= 0.5.0) — input data format
- [`sf`](https://cran.r-project.org/package=sf) — spatial operations
- [`terra`](https://cran.r-project.org/package=terra) — raster output
- [`units`](https://cran.r-project.org/package=units) — physical units handling

## Quick start

```r
library(move2UD)
library(move2)
library(sf)

# Load example data (fishers from LaPoint et al.)
fishers <- mt_read(mt_example())

# Prepare: remove empty geometries, select one individual, project
fishers <- fishers[!st_is_empty(fishers), ]
f1 <- fishers[mt_track_id(fishers) == "F1", ]
f1_proj <- st_transform(f1, mt_aeqd_crs(f1))

# Compute dBBMM utilisation distribution
ud <- mt_dbbmm_ud(f1_proj, location_error = 25,
                   window_size = 31, margin = 11)

# Plot
terra::plot(ud, main = "Fisher F1 — dBBMM Utilisation Distribution")
```

## Functions

### Variance estimation

| Function | Description |
|---|---|
| `mt_dbbmm_variance()` | Estimate dynamic Brownian motion variance using a sliding window with BIC-based breakpoint detection |
| `mt_dbgb_variance()` | Estimate dynamic bivariate Gaussian bridge variance, decomposed into parallel (along direction of travel) and orthogonal components |

### Utilisation distributions

| Function | Description |
|---|---|
| `mt_dbbmm_ud()` | Compute a dBBMM utilisation distribution. Accepts a `move2` object directly (estimates variance internally) or a pre-computed `mt_dbbmm_variance` object |
| `mt_dbgb_ud()` | Compute a dBGB utilisation distribution. Same interface as `mt_dbbmm_ud()` |

### Accessors

| Function | Description |
|---|---|
| `mt_motion_variance()` | Extract the estimated variance(s). Returns a numeric vector for dBBMM, or a data.frame with `para` and `orth` columns for dBGB |

## Usage

### dBBMM — isotropic variance

The dynamic Brownian bridge movement model estimates a single movement variance parameter that is allowed to change along the trajectory. This is appropriate when the direction of movement does not matter for the variance estimate.

```r
# Step 1: Estimate variance
var_dbbmm <- mt_dbbmm_variance(f1_proj, location_error = 25,
                                window_size = 31, margin = 11)
var_dbbmm
#> Dynamic Brownian Bridge Movement Model — variance estimate
#>   Locations: 1349
#>   Window size: 31, Margin: 11
#>   ...

# Inspect the variance over time
plot(mt_time(f1_proj), mt_motion_variance(var_dbbmm),
     type = "l", xlab = "Time", ylab = "BM variance")

# Step 2: Compute UD from variance object
ud <- mt_dbbmm_ud(var_dbbmm, location_error = 25)
terra::plot(ud)
```

Or in one step:

```r
ud <- mt_dbbmm_ud(f1_proj, location_error = 25,
                   window_size = 31, margin = 11)
```

### dBGB — directional variance

The dynamic bivariate Gaussian bridge decomposes the movement variance into a parallel component (along the direction of travel) and an orthogonal component (perpendicular to it). This provides more detail about the movement process — directional movement produces high parallel variance and low orthogonal variance, while random movement produces similar values in both.

```r
# Estimate directional variance
var_dbgb <- mt_dbgb_variance(f1_proj, location_error = 25,
                              margin = 15, window_size = 31)

# Extract and inspect
mv <- mt_motion_variance(var_dbgb)
head(mv)
#>       para      orth
#> 1       NA        NA
#> 2       NA        NA
#> ...

# Directionality index: >0 means directional, ~0 means Brownian
I_d <- (mv$para - mv$orth) / (mv$para + mv$orth)

# Compute UD
ud_bgb <- mt_dbgb_ud(var_dbgb, location_error = 25)
terra::plot(ud_bgb)
```

### Key parameters

| Parameter | Description | Guidance |
|---|---|---|
| `location_error` | GPS/telemetry error in map units (metres for projected CRS) | Scalar or per-location vector. Typical GPS: 10-30m |
| `window_size` | Number of locations in the sliding window (must be odd) | Larger = smoother variance, misses short behavioural changes. Smaller = more sensitive, more noise |
| `margin` | Minimum locations on each side of a breakpoint (must be odd) | Must be < `window_size / 2`. Larger = more robust breakpoints |
| `ext` | Extension factor for auto-generated raster extent | Increase if you get "grid not large enough" errors |
| `dim_size` | Number of cells along the longest axis of the raster | Higher = finer resolution, slower computation |

### Input requirements

All functions validate their input and produce informative error messages. The requirements are:

- **Single track**: The `move2` object must contain exactly one individual. Filter first with `dplyr::filter()` or process tracks in a loop.
- **Projected CRS**: Must not be longitude/latitude. Use `move2::mt_aeqd_crs()` for an azimuthal equidistant projection centred on the data:
  ```r
  data_proj <- st_transform(data, mt_aeqd_crs(data))
  ```
- **No empty geometries**: Remove failed GPS fixes first: `data <- data[!sf::st_is_empty(data), ]`
- **Ordered, unique timestamps**: Data must be sorted by time with no duplicate timestamps. Use `move2::mt_filter_unique()` if needed.
- **Positive location error**: `location_error` must be > 0 (not zero).
- **Window/margin constraints**: `window_size` and `margin` must be odd, and `window_size >= 2 * margin + 1`.

## Output

Variance functions return S3 objects (`mt_dbbmm_variance` or `mt_dbgb_variance`) that contain the estimated variances, breakpoint information, and the original track coordinates.

UD functions return a `terra::SpatRaster` where cell values sum to 1.0, representing the probability of the animal being in each cell. Standard `terra` functions work directly:

```r
# Plot
terra::plot(ud)

# Extract contours (e.g., 95% UD)
terra::as.contour(ud, levels = 0.95)

# Write to file
terra::writeRaster(ud, "fisher_ud.tif")

# Get volume (cumulative probability)
ud_sorted <- sort(terra::values(ud), decreasing = TRUE)
ud_cumsum <- cumsum(ud_sorted)
```

## Migration from `move`

If you previously used the `move` package, the mapping is:

| `move` function | `move2UD` function |
|---|---|
| `brownian.motion.variance.dyn()` | `mt_dbbmm_variance()` |
| `brownian.bridge.dyn()` | `mt_dbbmm_ud()` |
| `dynBGBvariance()` | `mt_dbgb_variance()` |
| `dynBGB()` | `mt_dbgb_ud()` |
| `getMotionVariance()` | `mt_motion_variance()` |
| `spTransform(x, center=TRUE)` | `sf::st_transform(x, move2::mt_aeqd_crs(x))` |

Key differences:
- Input is a `move2` object (not `Move`/`MoveStack`)
- Output UD is a `terra::SpatRaster` (not `RasterLayer`)
- Variance objects are S3 lists (not S4 classes inheriting from `SpatialPointsDataFrame`)
- Parameter names use `snake_case`: `location_error` (not `location.error`), `window_size` (not `window.size`)

## References

- Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K. (2012). A dynamic Brownian bridge movement model to estimate utilization distributions for heterogeneous animal movement. *Journal of Animal Ecology*, 81(4), 738-746. [doi:10.1111/j.1365-2656.2012.01955.x](https://doi.org/10.1111/j.1365-2656.2012.01955.x)

- Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian bridges: directional factorization of diffusion in Brownian bridge models. *Movement Ecology*, 2(1), 5. [doi:10.1186/2051-3933-2-5](https://doi.org/10.1186/2051-3933-2-5)

- Horne, J. S., Garton, E. O., Krone, S. M., & Lewis, J. S. (2007). Analyzing animal movements using Brownian bridges. *Ecology*, 88(9), 2354-2363. [doi:10.1890/06-0957.1](https://doi.org/10.1890/06-0957.1)

## Authors

- **Bart Kranstauber** (maintainer) — University of Amsterdam
- **Kamran Safi** — Max Planck Institute of Animal Behavior
- **Anne K. Scharf** — Max Planck Institute of Animal Behavior

## License

GPL (>= 3)
