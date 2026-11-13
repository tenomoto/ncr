# ncr: Interactive Exploratory Visualization of NetCDF in Base R

`ncr` is a lightweight, GrADS- and NCL-inspired S3 wrapper around `RNetCDF`, `fields`, and `maps` for rapid exploratory analysis and visualization of NetCDF datasets in base R graphics.

## Why `fields` + `maps` (Recommended over `sf` / `terra`)

1. **Zero C-Library Conflicts with `RNetCDF`**: On macOS and Linux, loading `terra` or `sf` alongside `RNetCDF` in the same R session can cause segmentation faults on exit or during file operations due to duplicate linked C library symbols (`libnetcdf`, `libhdf5`, `libproj`, `libgdal`). `fields` and `maps` are pure base-graphics companions with no conflicting native geospatial libraries.
2. **Unified Base R Coordinate System (`par("usr")`)**: `graphics::.filled.contour()` paired with `fields::image.plot(legend.only = TRUE)` and `maps::map(add = TRUE)` preserves the active user coordinate window (`par("usr")`), allowing contour lines (`contour(add = TRUE)`), coastlines (`maps::map(add = TRUE)`), and self-scaling vector arrows (`drawvec()`) to compose on the same canvas without layout distortion.

## Installation

Install from a local source checkout:

```r
install.packages(".", repos = NULL, type = "source")
```

Or from within R:

```r
install.packages(".", repos = NULL, type = "source")
```

Alternatively using pak:

```r
install.packages("pak") # skip if already installed
pak::pak("tenomoto/ncr")
```

## Quickstart

```r
library(ncr)

# 1. Open a NetCDF file (nc_example() creates a synthetic 4D dataset)
f  <- nc_example()
nc <- nc_open(f)

# Inspect dimensions, variables, and current global state
print(nc)
nc_meta(nc)

# 2. Instant 2D filled-contour map at default state (lev = 1, time = 1)
plot(nc$temp)
# or:
plot(nc, "temp")

# 3. Line contours (filled = FALSE) or filled contours via contour()
contour(nc, "slp", filled = FALSE)
contour(nc$temp, filled = TRUE)

# 4. Wind / current vectors and scalar + vector overlays
plot(nc, c("u", "v"), lev = 850, skip = 2)
plot(nc, "temp", uv = c("u", "v"), lev = 850, skip = 2)

# 5. Stateful GrADS-like slicing with nc_set() / nc_reset()
nc_set(lon = c(100, 220), lat = c(0, 65), lev = 500, time = "2026-01-02 00:00:00")
plot(nc$temp, uv = nc[c("u", "v")])
nc_reset()

# 6. Easy ad-hoc slicing syntax with `[`
plot(nc["temp", lon = c(110, 170), lat = c(15, 55), lev = 700, time = 4])

# 7. Map projections ("orthographic", "robinson", "mollweide", "npstere", "lambert", ...)
plot(nc, "temp", proj = "orthographic", lon0 = 140, lat0 = 30, uv = c("u", "v"), skip = 2)
plot(nc, "temp", proj = "npstere")
plot(nc, "temp", proj = "robinson")

# 8. Hovmoeller diagrams with automatic POSIXct time axis
nc_hovmoeller(nc, "temp", lat = 0)
# or via slice syntax:
plot(nc["temp", lat = 0, time = NULL], time_axis = "x")

# 9. Vertical cross-sections with dual Pressure [hPa] and Height [km] axes
nc_xsection(nc, "temp", lon = 140, zlog = TRUE)
# or via slice syntax:
plot(nc["u", lon = 140, lev = NULL])

# 10. Clean up
nc_close(nc)
```

## License

BSD 2-Clause License (`BSD_2_clause + file LICENSE`).
