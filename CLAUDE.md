# CLAUDE.md — move2UD: Dynamic Brownian Bridge Movement Models for move2

## Project overview

`move2UD` is an R package that provides dynamic Brownian Bridge Movement Model (dBBMM) and dynamic Bivariate Gaussian Bridge (dBGB) functionality for the `move2`/`sf`/`terra` ecosystem. It migrates the core UD estimation and movement variance algorithms from the `move` package (which depends on the retired `sp`/`raster` stack) into the modern spatial framework.

The package is authored by Bart Kranstauber, Kamran Safi, and Anne K. Scharf — the same team behind `move2`.

**GitHub repository**: https://github.com/kamransafi/move2UD

## Related project

This package was spun off from the book project "Analysing and Mapping of Animal Movement in R" (https://gitlab.com/anneks/AnimalMovementAnalysis) during the migration of chapters from `move`/`sp`/`raster` to `move2`/`sf`/`terra`. The book will use `move2UD` for dBBMM/dBGB examples once the package is stable.

## Architecture

### C layer (numerical core)

| File | Function | Purpose |
|------|----------|---------|
| `src/bgb_bbmm.c` | `dbbmm2()`, `bgb()` | Original grid kernels from `move` — preserved as reference |
| `src/bgb_bbmm_omp.c` | `dbbmm2_omp()`, `bgb_omp()` | **OpenMP-parallelised** versions used in production. Inner grid loops parallelised with `#pragma omp parallel for collapse(2)`. Guarded with `#ifdef _OPENMP` for CRAN compatibility |
| `src/bm_variance_c.c` | `bm_variance_c()`, `bm_variance_window_c()` | **C implementation of BM variance estimation**. Brent's method optimizer replaces R's `optimize()`. Processes entire sliding window (whole window + all breakpoints) in one `.Call()` |
| `src/bgbVar.c` | `llBGBvar()` | Bivariate Gaussian log-likelihood for BGB variance optimisation (from `move`) |
| `src/init.c` | — | Native routine registration for all 7 C functions |
| `src/Makevars` | — | Compiler flags: `$(SHLIB_OPENMP_CFLAGS)` for OpenMP support |

All C functions take raw numeric vectors and return numeric matrices. They have **no dependency on any R spatial package**.

### R layer (interface)

| File | Exports | Purpose |
|------|---------|---------|
| `R/utils.R` | internal helpers | `.extract_track_data()` (validates and extracts from move2), `.expand_loc_error()`, `.extcalc()`, `.make_raster()` |
| `R/bm_variance.R` | `bm_variance()` (internal) | Pure R reference implementation of BM variance estimation. Superseded by C version but retained for testing |
| `R/dbbmm_variance_dyn.R` | `mt_dbbmm_variance()`, `mt_motion_variance()` | Dynamic sliding-window variance with BIC breakpoint detection. Uses C kernel. Supports `parallel` and `cores` parameters |
| `R/dbbmm_ud.R` | `mt_dbbmm_ud()` | dBBMM utilisation distribution. Uses OpenMP kernel `dbbmm2_omp` |
| `R/dbgb_variance.R` | `delta_para_orth()`, `bgb_var_single()`, `bgb_var_break()` (internal) | BGB parallel/orthogonal decomposition and breakpoint detection. Uses sequential break search (O(n) instead of O(n²)) |
| `R/dbgb_variance_dyn.R` | `mt_dbgb_variance()` | Dynamic BGB variance. Supports `parallel` and `cores` parameters |
| `R/dbgb_ud.R` | `mt_dbgb_ud()` | dBGB utilisation distribution. Uses OpenMP kernel `bgb_omp` |

### Naming conventions

All exported functions follow the `move2` naming pattern:

- **`mt_` prefix** for all functions (matches `move2`: `mt_speed()`, `mt_distance()`, etc.)
- **Method abbreviation**: `dbbmm` / `dbgb`
- **`_variance` / `_ud`**: the computation type
- **`mt_motion_variance()`**: accessor following `mt_` + noun pattern
- **`location_error`**: consistent parameter name across all functions
- No naming collisions with `sf` (`st_` prefix), `terra` (no prefix), or `move2` (`mt_` prefix)

### Object model

Simple S3 lists replacing the old `move` S4 hierarchy:

- **`mt_dbbmm_variance`**: `variance`, `in_windows`, `interest`, `break_list`, `window_size`, `margin`, `track_data`
- **`mt_dbgb_variance`**: `para_sd`, `orth_sd`, `n_estim`, `seg_interest`, `margin`, `window_size`, `track_data`
- **UD output**: plain `terra::SpatRaster` (values sum to 1.0)

### Input validation

`.extract_track_data()` performs comprehensive validation on every call:

- Must be a `move2` object (`mt_is_move2()`)
- Must contain exactly one track (`mt_n_tracks() == 1`)
- Must be in a projected CRS (not lon/lat)
- No empty geometries (`st_is_empty()`)
- No NaN/NA/Inf coordinates
- Timestamps must be chronologically ordered
- No duplicate timestamps

Additional validation in the variance functions:
- `location_error` must be positive, no NAs, correct length
- `window_size` and `margin` must be odd
- `window_size >= 2 * margin + 1`

UD functions validate that the grid computation produces finite results.

## Build and check

```bash
R CMD build .
R CMD INSTALL move2UD_*.tar.gz
R CMD check move2UD_*.tar.gz

# Or in R:
devtools::install()
devtools::test()
devtools::check()
```

### Dependencies

All on CRAN:
- `move2` — input objects
- `sf` — spatial operations, coordinate extraction, CRS handling, empty geometry detection
- `terra` — raster output for utilisation distributions
- `units` — physical units for time lags

Suggested: `parallel` (for `mclapply`), `testthat` (>= 3.0.0)

## Coding standards

### R code style
- **snake_case** for all function and variable names
- Internal functions prefixed with `.` (e.g., `.extract_track_data`)
- Use `sf::st_coordinates()`, `move2::mt_time()`, etc. — never `sp::coordinates()` or `move::timestamps()`
- Return `terra::SpatRaster` for raster outputs — never `raster::RasterLayer`
- Use S3 classes and generics — no S4

### C code
- Original kernels in `src/bgb_bbmm.c` and `src/bgbVar.c` are **preserved verbatim from `move`** as reference. Do not modify.
- Production code uses `src/bgb_bbmm_omp.c` (OpenMP versions) and `src/bm_variance_c.c` (C variance estimation)
- All OpenMP pragmas must be guarded with `#ifdef _OPENMP` for CRAN compatibility
- `R_CheckUserInterrupt()` must not be called from within OpenMP parallel regions — call it in the outer sequential loop

### Input validation rules
- All user-facing functions must call `.extract_track_data()` which performs all input checks
- Error messages must tell the user what to do (e.g., "Remove them first: x <- x[!sf::st_is_empty(x), ]")
- Never silently accept invalid input — always `stop()` with `call. = FALSE`

## Testing

### Current state
70 tests, 0 failures, 0 warnings.

### Test files
| File | Tests | Purpose |
|------|-------|---------|
| `test-utils.R` | 6 | `.expand_loc_error`, `.extcalc`, `.make_raster` |
| `test-bm-variance.R` | 4 | Core variance: structure, stationary, tortuous, input validation |
| `test-delta-para-orth.R` | 5 | Geometric decomposition at 0°, 45°, 90°, multiple points, same location |
| `test-dbbmm-variance-dyn.R` | 6 | Variance: projected data, lon/lat rejection, even window, too-small track, accessor, equivalence with `move` |
| `test-dbbmm-ud.R` | 3 | UD: returns SpatRaster summing to 1, pre-computed variance, cell size |
| `test-dbgb.R` | 4 | dBGB: variance structure, accessor, UD normalisation, equivalence with `move` |
| `test-input-validation.R` | 5 | Multi-track, empty geometries, lon/lat, window/margin, non-move2 input |

### Numerical equivalence
Tests compare against reference fixtures (`.rds` files in `tests/testthat/fixtures/`) generated from the `move` package. Tolerance is 1% to account for projection centering differences between `sp::spTransform(center=TRUE)` and `sf::st_transform()` with custom AEQD.

To regenerate fixtures (requires both `move` and `move2` installed):
```r
source("tests/testthat/helper-generate-fixtures.R")
generate_fixtures()
```

## Performance

Benchmarks on fisher F1 (1349 locations, 8 cores):

| Operation | `move` (old) | `move2UD` | Speedup |
|---|---|---|---|
| dBBMM variance (w=31) | 5.0s | 0.46s (C + parallel) | **10.9x** |
| dBBMM variance (w=71) | 14.2s | 0.74s (C + parallel) | **19.3x** |
| dBBMM UD (dim=100) | 5.2s | 0.4s (OpenMP) | **12.9x** |
| dBGB variance (w=31) | 30.1s | 6.1s (parallel) | **5.0x** |

### Performance architecture
1. **C Brent's method** (`bm_variance_c.c`): Replaces R's `optimize()` for BM variance — 8-9x speedup on inner loop
2. **Parallel window loops**: `parallel::mclapply` for both dBBMM and dBGB variance — 2-2.5x additional
3. **Sequential break search** for dBGB: O(n) instead of O(n²) — 1.5x algorithmic improvement
4. **OpenMP grid kernels** (`bgb_bbmm_omp.c`): `#pragma omp parallel for collapse(2)` on inner cell loops — 12.9x for UD computation

## Remaining work (in priority order)

1. **Roxygen documentation** — full `@param`, `@return`, `@examples` for all exported functions. Currently missing `.Rd` files (R CMD check WARNING)
2. **S3 method signature alignment** — ensure method signatures match generics (R CMD check WARNING)
3. **Vignette** — worked example from `move2` object to UD plot
4. **Multi-track convenience** — wrapper that loops over tracks and returns a list of UDs
5. **dBGB C optimizer** — move the L-BFGS-B optimization in `bgb_var_break()` to C for further speedup (currently the dBGB bottleneck)
6. **`\donttest` examples** — convert from `\dontrun` so examples are checked during `R CMD check --run-donttest`

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
