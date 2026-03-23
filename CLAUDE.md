# CLAUDE.md — move2UD: Dynamic Brownian Bridge Movement Models for move2

## Project overview

`move2UD` is an R package that provides dynamic Brownian Bridge Movement Model (dBBMM) and dynamic Bivariate Gaussian Bridge (dBGB) functionality for the `move2`/`sf`/`terra` ecosystem. It migrates the core UD estimation and movement variance algorithms from the `move` package (which depends on the retired `sp`/`raster` stack) into the modern spatial framework.

The package is authored by Bart Kranstauber, Kamran Safi, and Anne K. Scharf — the same team behind `move2`. The C computational kernels are preserved from the original `move` implementation; the R interface layer has been rewritten.

**GitHub repository**: https://github.com/kamransafi/move2UD

## Related project

This package was spun off from the book project "Analysing and Mapping of Animal Movement in R" (https://gitlab.com/anneks/AnimalMovementAnalysis) during the migration of chapters from `move`/`sp`/`raster` to `move2`/`sf`/`terra`. The book will use `move2UD` for dBBMM/dBGB examples once the package is stable.

## Architecture

### C layer (numerical core)

The heavy computation is in C, preserved verbatim from the `move` package:

| File | Function | Purpose |
|------|----------|---------|
| `src/bgb_bbmm.c` | `dbbmm2()` | Isotropic dBBMM grid computation — evaluates Brownian bridge probability density on a raster grid |
| `src/bgb_bbmm.c` | `bgb()` | Anisotropic BGB grid computation — decomposes into parallel/orthogonal components per grid cell |
| `src/bgbVar.c` | `llBGBvar()` | Bivariate Gaussian log-likelihood for BGB variance optimisation |
| `src/init.c` | — | Native routine registration |

These C functions take raw numeric vectors (coordinates, variances, times, grid coordinates) and return numeric matrices. They have **no dependency on any R spatial package** — all coupling to spatial objects is in the R layer.

### R layer (interface)

| File | Exports | Purpose |
|------|---------|---------|
| `R/bm_variance.R` | `bm_variance()` (internal) | Core leave-one-out likelihood BM variance estimation |
| `R/dbbmm_variance_dyn.R` | `mt_dbbmm_variance()`, `mt_motion_variance()` | Dynamic sliding-window variance with BIC breakpoint detection |
| `R/dbbmm_ud.R` | `mt_dbbmm_ud()` | dBBMM utilisation distribution computation |
| `R/dbgb_variance.R` | `delta_para_orth()`, `bgb_var_single()`, `bgb_var_break()` (internal) | BGB parallel/orthogonal decomposition and breakpoint detection |
| `R/dbgb_variance_dyn.R` | `mt_dbgb_variance()` | Dynamic BGB variance estimation |
| `R/dbgb_ud.R` | `mt_dbgb_ud()` | dBGB utilisation distribution computation |
| `R/utils.R` | internal helpers | `.extract_track_data()`, `.expand_loc_error()`, `.extcalc()`, `.make_raster()` |

### Naming conventions

All exported functions follow the `move2` naming pattern:

- **`mt_` prefix** for all functions (matches `move2` convention: `mt_speed()`, `mt_distance()`, etc.)
- **Method abbreviation**: `dbbmm` / `dbgb` (well-known in movement ecology)
- **`_variance` / `_ud`**: the computation type
- **`mt_motion_variance()`**: accessor following the `mt_` + noun pattern (like `mt_time()`, `mt_track_data()`)
- **`location_error`**: consistent parameter name across all functions (not `loc_err`)

### Object model

The old `move` package used a deep S4 class hierarchy (`dBMvariance` → `.MoveTrackSingle` → `SpatialPointsDataFrame`, `DBBMM` → `.UD` → `RasterLayer`). This package uses simple S3 lists:

- **`mt_dbbmm_variance`**: list with `variance`, `in_windows`, `interest`, `break_list`, `window_size`, `margin`, `track_data`
- **`mt_dbgb_variance`**: list with `para_sd`, `orth_sd`, `n_estim`, `seg_interest`, `margin`, `window_size`, `track_data`
- **UD output**: plain `terra::SpatRaster` (no custom class — use `terra` functions directly)

### CRS handling

Users should project their data before calling any `move2UD` function. The recommended approach uses `move2::mt_aeqd_crs()`:

```r
data_proj <- sf::st_transform(data, move2::mt_aeqd_crs(data))
```

This creates an azimuthal equidistant projection centred on the data, matching what the old `move::spTransform(x, center=TRUE)` did.

## Build and check

```bash
# Build
R CMD build .

# Check
R CMD check move2UD_*.tar.gz

# Or in R:
devtools::check()
devtools::test()
```

### Dependencies

All dependencies are on CRAN:
- `move2` — input objects
- `sf` — spatial operations, coordinate extraction, CRS handling
- `terra` — raster output for utilisation distributions
- `units` — physical units for time lags
- `methods` — for `.Call()` to C routines

Test dependencies: `testthat` (>= 3.0.0)

## Coding standards

### R code style
- **snake_case** for all function and variable names (e.g., `dbbmm_variance_dyn`, `location_error`)
- Internal functions prefixed with `.` (e.g., `.extract_track_data`)
- Use `sf::st_coordinates()`, `move2::mt_time()`, etc. — never `sp::coordinates()` or `move::timestamps()`
- Return `terra::SpatRaster` for raster outputs — never `raster::RasterLayer`
- Use S3 classes and generics — no S4

### C code
- The C kernels in `src/bgb_bbmm.c` and `src/bgbVar.c` are **preserved from the `move` package** and should not be modified unless fixing a bug or adding parallelisation (e.g., OpenMP). Any changes must preserve numerical equivalence with the original `move` implementation.
- New C/C++ code (e.g., Rcpp reimplementations) should go in new files, not modify the originals.

### Input validation
- All user-facing functions must check that the input is a `move2` object in a projected CRS (not lon/lat)
- `location_error` must be validated: positive, no NAs, correct length
- `window_size` and `margin` must be odd
- Informative error messages referencing what the user should do (e.g., "Transform to a projected CRS first using sf::st_transform()")

## Testing strategy

### Numerical equivalence tests
The primary correctness criterion is that `move2UD` produces **identical numerical results** to the original `move` package functions for the same input data. Tests should:

1. Load a test dataset as both a `move` object and a `move2` object
2. Run the `move` function (e.g., `brownian.motion.variance.dyn()`) and the `move2UD` function (`dbbmm_variance_dyn()`)
3. Compare results with `expect_equal()` using a tolerance appropriate for floating-point arithmetic (typically `tolerance = 1e-8`)

This applies to:
- Variance estimates (`dbbmm_variance_dyn` vs `brownian.motion.variance.dyn`)
- UD raster values (`dbbmm_ud` vs `brownian.bridge.dyn`)
- BGB parallel/orthogonal variances (`dbgb_variance_dyn` vs `dynBGBvariance`)
- BGB UD values (`dbgb_ud` vs `dynBGB`)

### Unit tests for components
- `bm_variance()`: known inputs → known likelihood and variance
- `delta_para_orth()`: geometric decomposition with simple test cases (0°, 45°, 90° angles)
- `.extract_track_data()`: correct extraction from move2 objects; error on lon/lat
- `.expand_loc_error()`: scalar expansion, length mismatch error, NA error
- `.make_raster()`: correct extent, resolution, CRS propagation
- `get_motion_variance()`: correct accessor for both `dbbmm_var` and `dbgb_var`

### Edge case tests
- Single-individual vs multi-individual input (multi should error or be handled)
- Very short tracks (fewer locations than window_size → informative error)
- Tracks with irregular time intervals
- Location error as scalar vs vector
- Window size equal to track length (boundary case)

### Performance regression tests
- Benchmark key functions on a standard dataset and track execution time
- Ensure performance improvements (parallelisation, Rcpp) don't change results

### Test data
Use `move2::mt_read(move2::mt_example())` as the primary test dataset (fisher data, bundled with `move2`). For numerical equivalence tests against `move`, also keep a small reference dataset with pre-computed `move` package results stored as `.rds` files in `tests/testthat/fixtures/`.

## Planned improvements (in order of priority)

1. **Parallelise the window loop** — `parallel::mclapply` or `future.apply::future_lapply` for the sliding window in `dbbmm_variance_dyn()` and `dbgb_variance_dyn()`. Each window position is independent.
2. **Rcpp reimplementation of `bm_variance()`** — the inner likelihood evaluation and optimisation are the main bottleneck. Moving from R's `optimize()` to a C++ Brent's method would give 10-50× speedup.
3. **OpenMP in C kernels** — the inner grid loops in `dbbmm2()` and `bgb()` are trivially parallelisable.
4. **Multi-track support** — accept `move2` objects with multiple tracks; process per-track and return a list or stack.
5. **Roxygen documentation** — full `@param`, `@return`, `@examples` for all exported functions.
6. **Vignette** — worked example from `move2` object to UD plot, showing both dBBMM and dBGB workflows.

## Authoritative references

### Methodology
- Horne, J. S., Garton, E. O., Krone, S. M., & Lewis, J. S. (2007). Analyzing animal movements using Brownian bridges. *Ecology*, 88(9), 2354-2363.
- Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M., & Safi, K. (2012). A dynamic Brownian bridge movement model to estimate utilization distributions for heterogeneous animal movement. *Journal of Animal Ecology*, 81(4), 738-746.
- Kranstauber, B., Safi, K., & Bartumeus, F. (2014). Bivariate Gaussian bridges: directional factorization of diffusion in Brownian bridge models. *Movement Ecology*, 2(1), 5.

### Original implementation
- `move` package source: https://cran.r-project.org/package=move (v4.2.7)
- Files migrated: `R/brownianmotionvariancedyn.R`, `R/brownianbridgedyn.R`, `R/dBGB.R`, `R/getMotionVariance.R`, `src/bgb_bbmm.c`, `src/bgbVar.c`

### Package ecosystem
- `move2`: https://bartk.gitlab.io/move2/ — input object format
- `sf`: https://r-spatial.github.io/sf/ — spatial vector data
- `terra`: https://rspatial.github.io/terra/ — raster output
- `units`: https://r-quantities.github.io/units/ — physical units
