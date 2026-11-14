#' Draw Self-Scaling Vector Arrows
#'
#' Plots a 2D vector field `(u, v)` on an existing base R graphics device,
#' scaling arrow lengths automatically to the grid spacing and plot window
#' (`par("usr")`). Accepts either two `"ncr_slice"` objects `drawvec(u, v)`,
#' a single vector `"ncr_slice"` object `drawvec(uv)`, or explicit numeric
#' coordinates and matrices `drawvec(x, y, u, v)`.
#'
#' @param x Numeric vector of X coordinates, an `"ncr_slice"` object for the
#'   `u` component (when `y` is an `"ncr_slice"` for `v`), or a vector
#'   `"ncr_slice"` object (created via `nc[c("u", "v")]`).
#' @param y Numeric vector of Y coordinates, an `"ncr_slice"` object for the
#'   `v` component, or `NULL` when `x` is a vector `"ncr_slice"`.
#' @param u,v Numeric matrices (or vectors) of horizontal and vertical vector
#'   components matching `length(x)` by `length(y)` (ignored when `x` is an
#'   `"ncr_slice"`).
#' @param length Arrowhead length in inches passed to [graphics::arrows()]
#'   (default `0.05`).
#' @param uvref Reference vector magnitude used for scaling. If `NA` (default),
#'   computed automatically as `max(abs(c(u, v)), na.rm = TRUE)`.
#' @param skip Integer subsampling stride (`1L` plots every grid point, `2L`
#'   plots every second point along `x` and `y`).
#' @param col Arrow color (default `"black"`).
#' @param proj,lon0,lat0 Optional map projection settings when `x` is an
#'   `"ncr_slice"` (defaults to `nc_get()`).
#' @param ref_legend Logical; if `TRUE`, draws a reference vector arrow above
#'   the `"Ref vector: ..."` label in the bottom-right margin.
#' @param ref_units Optional character string for reference vector units.
#' @param ... Additional graphical parameters passed to [graphics::arrows()].
#' @return Invisibly returns the reference magnitude `uvref`.
#' @examples
#' x <- seq(100, 160, by = 10)
#' y <- seq(10, 50, by = 10)
#' u <- outer(x, y, function(lon, lat) 10 * cos(lat * pi / 180))
#' v <- outer(x, y, function(lon, lat) 5 * sin(lon * pi / 180))
#' plot(NA, xlim = range(x), ylim = range(y), xlab = "lon", ylab = "lat")
#' drawvec(x, y, u, v)
#' @export
drawvec <- function(x, y = NULL, u = NULL, v = NULL,
                    length = 0.05, uvref = NA, skip = 1L, col = "black",
                    proj = NULL, lon0 = NULL, lat0 = NULL,
                    ref_legend = FALSE, ref_units = "", ...) {
  if (inherits(x, "ncr_slice")) {
    s <- x
    if (inherits(y, "ncr_slice")) {
      if (!identical(dim(x$z), dim(y$z))) {
        stop("'u' and 'v' slices must have identical dimensions.", call. = FALSE)
      }
      u_mat <- x$z
      v_mat <- y$z
    } else if (isTRUE(s$has_uv) && !is.null(s$u) && !is.null(s$v)) {
      u_mat <- s$u
      v_mat <- s$v
    } else {
      stop("Provide either drawvec(u_slice, v_slice) or a vector slice from nc[c('u', 'v')].",
           call. = FALSE)
    }

    if (is.null(proj)) proj <- nc_get()$proj
    if (is.null(lon0)) lon0 <- nc_get()$lon0
    if (is.null(lat0)) lat0 <- nc_get()$lat0
    proj_norm <- tolower(trimws(as.character(proj)[1L]))
    is_proj <- (s$plot_type == "lonlat") &&
      !(proj_norm %in% c("lonlat", "cylindrical", "equirectangular", "eqc", "none", ""))
    if (is_proj) {
      return(.draw_projected_vectors(s$x, s$y, u_mat, v_mat, proj = proj_norm,
                                     lon0 = lon0, lat0 = lat0,
                                     skip = skip, uvref = uvref, length = length,
                                     col = col, ref_legend = ref_legend,
                                     ref_units = if (nzchar(ref_units)) ref_units else if (!is.null(s$uv_units)) s$uv_units else s$units,
                                     ...))
    }
    x <- s$x
    y <- s$y
    u <- u_mat
    v <- v_mat
  } else if (is.null(y) || is.null(u) || is.null(v)) {
    stop("Numeric drawvec() requires x, y, u, and v (or pass ncr_slice objects).",
         call. = FALSE)
  }

  skip <- max(1L, as.integer(skip)[1L])
  if (skip > 1L && is.matrix(u) && is.matrix(v)) {
    ix <- seq.int(1L, length(x), by = skip)
    iy <- seq.int(1L, length(y), by = skip)
    x <- x[ix]
    y <- y[iy]
    u <- u[ix, iy, drop = FALSE]
    v <- v[ix, iy, drop = FALSE]
  }
  nx <- length(x)
  ny <- length(y)
  usr <- graphics::par("usr")
  d_grid <- min((usr[2L] - usr[1L]) / nx, (usr[4L] - usr[3L]) / ny)
  if (is.na(uvref)) uvref <- max(abs(c(u, v)), na.rm = TRUE)
  if (!is.finite(uvref) || uvref <= 0) return(invisible(uvref))
  d <- d_grid / uvref
  u <- 0.5 * as.vector(u) * d
  v <- 0.5 * as.vector(v) * d
  X <- rep(x, times = length(y))
  Y <- rep(y, each  = length(x))
  pin <- graphics::par("pin")
  dx_in <- (2 * u) * pin[1L] / max(abs(usr[2L] - usr[1L]), 1e-12)
  dy_in <- (2 * v) * pin[2L] / max(abs(usr[4L] - usr[3L]), 1e-12)
  len_in <- sqrt(dx_in^2 + dy_in^2)
  keep <- !is.na(u) & !is.na(v) & !is.na(X) & !is.na(Y) & (len_in >= 1e-3)
  if (any(keep)) {
    graphics::arrows(X[keep] - u[keep], Y[keep] - v[keep],
                     X[keep] + u[keep], Y[keep] + v[keep],
                     length = length, col = col, ...)
  }
  if (isTRUE(ref_legend)) {
    .draw_ref_vector(uvref, ref_len_usr = d_grid, length = length,
                     units = ref_units, col = col)
  }
  invisible(uvref)
}

.draw_ref_vector <- function(ref, ref_len_usr, length = 0.05, units = "", col = "black") {
  if (!is.finite(ref) || ref <= 0 || !is.finite(ref_len_usr) || ref_len_usr <= 0) {
    return(invisible(NULL))
  }
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  line_y <- graphics::par("csi") * abs(usr[4L] - usr[3L]) / max(pin[2L], 1e-6)
  y_arr <- usr[3L] - 2.3 * line_y
  x_end <- usr[2L]
  x_start <- x_end - ref_len_usr
  graphics::arrows(x_start, y_arr, x_end, y_arr,
                   length = length, col = col, xpd = NA)
  u_str <- if (is.character(units) && length(units) == 1L && nzchar(units)) {
    paste0(" ", units)
  } else {
    ""
  }
  graphics::mtext(sprintf("Ref vector: %s%s", signif(ref, 3), u_str),
                  side = 1, line = 3.2, adj = 1, cex = 0.78, col = "gray25")
  invisible(NULL)
}

#' Convert Between Pressure and Standard-Atmosphere Height
#'
#' Converts pressure (`hPa`) to approximate log-pressure scale height (`km`)
#' and vice-versa for dual-axis vertical cross-sections.
#'
#' @param p Numeric vector of pressure in `hPa`.
#' @param z Numeric vector of altitude/height in `km`.
#' @param p0 Reference surface pressure in `hPa` (default `1013.25`).
#' @param H Scale height in `km` (default `7.2`).
#' @return Numeric vector of height in `km` (`pres2alt`) or pressure in `hPa`
#'   (`alt2pres`).
#' @examples
#' pres2alt(c(1000, 850, 500, 200, 100))
#' alt2pres(c(0, 5, 10, 16))
#' @export
pres2alt <- function(p, p0 = 1013.25, H = 7.2) {
  -H * log(p / p0)
}

#' @rdname pres2alt
#' @export
alt2pres <- function(z, p0 = 1013.25, H = 7.2) {
  p0 * exp(-z / H)
}

#' Plot or Contour a NetCDF Dataset or Slice
#'
#' High-level S3 plotting methods for `"ncr"` and `"ncr_slice"` objects.
#' Automatically detects the slice geometry and produces:
#' * **2D Horizontal Maps (`lon`-`lat`)**: Filled or line contours, overlaid
#'   world coastlines (`maps::map()`), optional map projections (`proj`), and
#'   self-scaling wind/current vectors (`drawvec()`).
#' * **Hovmoeller Diagrams (`lon`-`time` or `lat`-`time`)**: Filled/line
#'   contours with formatted `POSIXct` time axes and degree labels.
#' * **Vertical Cross-Sections (`lon`-`lev` or `lat`-`lev`)**: Inverted
#'   log-pressure or linear-pressure axis on the left and standard-atmosphere
#'   height (`km`) on the right (or height on the left and pressure on the
#'   right).
#' * **1D Profiles / Time Series**: Line plots with formatted time or
#'   pressure-height axes.
#'
#' @param x An `"ncr"` or `"ncr_slice"` object.
#' @param var Optional variable name(s) when `x` is an `"ncr"` object:
#'   length 1 for scalar, length 2 for vector `c("u", "v")`, or length 3 for
#'   scalar + vector `c("temp", "u", "v")`.
#' @param lon,lat,lev,time Optional dimension selectors passed to [nc_slice()]
#'   (for `contour.ncr()`).
#' @param filled Logical; if `TRUE` (default for scalar `plot()`), draws filled
#'   contours with a colorbar via [fields::image.plot()].
#' @param contour Logical; if `TRUE`, draws contour lines over filled contours.
#' @param map Map overlay specification for `lon`-`lat` slices (defaults to
#'   `nc_get()$map` when `NULL`):
#'   * `TRUE`: draws coarse world coastlines via `maps::map("world")`.
#'   * `FALSE`: no map overlay.
#'   * `"hires"` or `"worldHires"`: high-resolution world coastlines via the
#'     suggested `mapdata` package.
#'   * `"japan"`, `"china"`, `"rivers"`, or any other `maps`/`mapdata` database
#'     name (e.g., `"world2"`, `"usa"`, `"state"`).
#'   * A custom `list(x, y)` of coastline coordinates (in degrees).
#' @param proj Character string specifying map projection (e.g., `"lonlat"`,
#'   `"orthographic"`, `"stereographic"`, `"npstere"`, `"spstere"`,
#'   `"robinson"`, `"mollweide"`, `"mercator"`, `"lambert"`). Defaults to the
#'   global state `nc_get()$proj`.
#' @param lon0,lat0 Optional center longitude and latitude for map projections.
#' @param nlevels Approximate number of contour levels (default `12`).
#' @param levels Optional explicit numeric vector of contour break levels.
#' @param col Color vector or palette function. Defaults to a diverging/spectral
#'   HCL palette (`hcl.colors(..., "RdYlBu", rev = TRUE)`).
#' @param zlog Logical; for pressure cross-sections, whether to space the
#'   vertical axis logarithmically in pressure (default `TRUE`).
#' @param time_axis For Hovmoeller diagrams, which axis carries time: `"y"`
#'   (default) or `"x"`.
#' @param uv Optional character vector `c("u", "v")` (when `x` is `"ncr"`) or
#'   `"ncr_slice"` vector object to overlay arrows on a scalar plot.
#' @param skip Integer subsampling factor for vector arrows (default `1L`).
#' @param uvref Reference vector length passed to [drawvec()].
#' @param length Arrowhead length in inches passed to [drawvec()] (default
#'   `0.05`).
#' @param add Logical; if `TRUE`, adds contours or vectors to the current plot.
#' @param main,sub,xlab,ylab Optional plot title, subtitle, and axis labels.
#' @param ... Additional slicing arguments (`lon`, `lat`, `lev`, `time`) when
#'   `x` is an `"ncr"` object, or graphical parameters passed to underlying plot
#'   functions.
#' @return Invisibly returns the plotted `"ncr_slice"` object.
#' @examples
#' f <- nc_example()
#' nc <- nc_open(f)
#'
#' # 1. Quick filled contour map with coastlines
#' plot(nc, "temp")
#'
#' # 2. Scalar + wind vector overlay at 850 hPa
#' plot(nc, "temp", uv = c("u", "v"), lev = 850)
#'
#' # 3. Orthographic projection
#' plot(nc, "temp", proj = "orthographic", lon0 = 140, lat0 = 30)
#'
#' # 4. Hovmoeller diagram along the equator
#' nc_hovmoeller(nc, "temp", lat = 0)
#'
#' # 5. Vertical pressure-height cross-section along 140E
#' nc_xsection(nc, "temp", lon = 140)
#'
#' nc_close(nc)
#' @export
plot.ncr <- function(x,
                     var = NULL,
                     ...,
                     uv = NULL,
                     filled = NULL,
                     contour = TRUE,
                     map = NULL,
                     proj = NULL,
                     lon0 = NULL,
                     lat0 = NULL,
                     nlevels = 12,
                     levels = NULL,
                     col = NULL,
                     zlog = TRUE,
                     time_axis = "y",
                     skip = 1L,
                     uvref = NA,
                     length = 0.05,
                     add = FALSE,
                     main = NULL,
                     sub = NULL,
                     xlab = NULL,
                     ylab = NULL) {
  uv_slice <- if (inherits(uv, "ncr_slice")) uv else NULL
  uv_names <- if (is.character(uv)) uv else NULL
  s <- nc_slice(x, var = var, uv = uv_names, ...)
  if (!is.null(uv_slice)) {
    s$u <- if (!is.null(uv_slice$u)) uv_slice$u else uv_slice$z
    s$v <- uv_slice$v
    s$has_uv <- TRUE
  }
  plot.ncr_slice(
    s,
    filled    = filled,
    contour   = contour,
    map       = map,
    proj      = proj,
    lon0      = lon0,
    lat0      = lat0,
    nlevels   = nlevels,
    levels    = levels,
    col       = col,
    zlog      = zlog,
    time_axis = time_axis,
    skip      = skip,
    uvref     = uvref,
    length    = length,
    add       = add,
    main      = main,
    sub       = sub,
    xlab      = xlab,
    ylab      = ylab
  )
}

.reslice_if_needed <- function(x, lon, lat, lev, time) {
  has_override <- !is_missing_arg(lon) || !is_missing_arg(lat) ||
    !is_missing_arg(lev) || !is_missing_arg(time)
  if (!has_override) return(x)

  if (is.null(x$source) || isTRUE(x$source$nc$env$closed)) {
    stop("Cannot re-slice this 'ncr_slice' object (dataset is closed or slice is a derived expression).",
         call. = FALSE)
  }
  src <- x$source
  nc_slice(
    src$nc,
    var  = src$var,
    lon  = if (is_missing_arg(lon))  src$lon  else lon,
    lat  = if (is_missing_arg(lat))  src$lat  else lat,
    lev  = if (is_missing_arg(lev))  src$lev  else lev,
    time = if (is_missing_arg(time)) src$time else time,
    uv   = src$uv
  )
}

#' @rdname plot.ncr
#' @export
plot.ncr_slice <- function(x,
                           ...,
                           lon = missing_arg(),
                           lat = missing_arg(),
                           lev = missing_arg(),
                           time = missing_arg(),
                           uv = NULL,
                           filled = NULL,
                           contour = TRUE,
                           map = NULL,
                           proj = NULL,
                           lon0 = NULL,
                           lat0 = NULL,
                           nlevels = 12,
                           levels = NULL,
                           col = NULL,
                           zlog = TRUE,
                           time_axis = "y",
                           skip = 1L,
                           uvref = NA,
                           length = 0.05,
                           add = FALSE,
                           main = NULL,
                           sub = NULL,
                           xlab = NULL,
                           ylab = NULL) {
  x <- .reslice_if_needed(x, lon = lon, lat = lat, lev = lev, time = time)
  if (inherits(uv, "ncr_slice")) {
    x$u <- if (!is.null(uv$u)) uv$u else uv$z
    x$v <- uv$v
    x$has_uv <- TRUE
  }
  if (is.null(filled)) {
    filled <- !isTRUE(x$is_vector)
  }
  if (is.null(map)) {
    map <- nc_get()$map
  }
  if (is.null(proj)) {
    proj <- nc_get()$proj
  }
  if (is.null(lon0)) lon0 <- nc_get()$lon0
  if (is.null(lat0)) lat0 <- nc_get()$lat0

  if (x$plot_type == "1d" || length(x$active_dims) == 1L) {
    .plot_1d_slice(x, main = main, sub = sub, xlab = xlab, ylab = ylab, zlog = zlog,
                   add = add, ...)
    return(invisible(x))
  }

  if (x$plot_type == "hovmoeller") {
    .plot_hovmoeller_slice(
      x, filled = filled, contour = contour, nlevels = nlevels,
      levels = levels, col = col, time_axis = time_axis,
      skip = skip, uvref = uvref, length = length,
      main = main, sub = sub, xlab = xlab, ylab = ylab, ...
    )
    return(invisible(x))
  }

  if (x$plot_type == "xsection") {
    .plot_xsection_slice(
      x, filled = filled, contour = contour, nlevels = nlevels,
      levels = levels, col = col, zlog = zlog,
      skip = skip, uvref = uvref, length = length,
      main = main, sub = sub, xlab = xlab, ylab = ylab, ...
    )
    return(invisible(x))
  }

  # 2D horizontal map ("lonlat") or generic 2D slice
  .plot_2d_slice(
    x, filled = filled, contour = contour, map = map,
    proj = proj, lon0 = lon0, lat0 = lat0,
    nlevels = nlevels, levels = levels, col = col,
    skip = skip, uvref = uvref, length = length,
    add = add, main = main, sub = sub, xlab = xlab, ylab = ylab, ...
  )
  invisible(x)
}

#' @rdname plot.ncr
#' @export
contour.ncr <- function(x, var = NULL,
                        lon = missing_arg(),
                        lat = missing_arg(),
                        lev = missing_arg(),
                        time = missing_arg(),
                        uv = NULL,
                        filled = FALSE,
                        ...) {
  s <- nc_slice(x, var = var, lon = lon, lat = lat, lev = lev, time = time, uv = uv)
  contour.ncr_slice(s, filled = filled, ...)
}

#' @rdname plot.ncr
#' @export
contour.ncr_slice <- function(x,
                              lon = missing_arg(),
                              lat = missing_arg(),
                              lev = missing_arg(),
                              time = missing_arg(),
                              filled = FALSE,
                              contour = !filled,
                              map = NULL,
                              proj = NULL,
                              lon0 = NULL,
                              lat0 = NULL,
                              nlevels = 12,
                              levels = NULL,
                              col = NULL,
                              add = FALSE,
                              main = NULL,
                              sub = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              ...) {
  x <- .reslice_if_needed(x, lon = lon, lat = lat, lev = lev, time = time)
  plot.ncr_slice(
    x,
    filled  = filled,
    contour = if (filled) contour else TRUE,
    map     = if (isTRUE(add)) FALSE else map,
    proj    = proj,
    lon0    = lon0,
    lat0    = lat0,
    nlevels = nlevels,
    levels  = levels,
    col     = col,
    add     = add,
    main    = main,
    sub     = sub,
    xlab    = xlab,
    ylab    = ylab,
    ...
  )
}

#' @rdname plot.ncr
#' @export
points.ncr_slice <- function(x, skip = 1L, uvref = NA, length = 0.05, col = "black",
                             proj = NULL, lon0 = NULL, lat0 = NULL, ...) {
  if (!isTRUE(x$has_uv) || is.null(x$u) || is.null(x$v)) {
    stop("'points.ncr_slice' requires a vector slice (created with c('u', 'v') or uv = c('u', 'v')).",
         call. = FALSE)
  }
  if (is.null(proj)) proj <- nc_get()$proj
  if (is.null(lon0)) lon0 <- nc_get()$lon0
  if (is.null(lat0)) lat0 <- nc_get()$lat0
  proj_norm <- tolower(trimws(as.character(proj)[1L]))
  is_proj <- (x$plot_type == "lonlat") &&
    !(proj_norm %in% c("lonlat", "cylindrical", "equirectangular", "eqc", "none", ""))
  if (is_proj) {
    .draw_projected_vectors(x$x, x$y, x$u, x$v, proj = proj_norm,
                            lon0 = lon0, lat0 = lat0,
                            skip = skip, uvref = uvref, length = length)
  } else {
    drawvec(x$x, x$y, x$u, x$v, length = length, uvref = uvref, skip = skip, col = col, ...)
  }
}

#' Hovmoeller Diagram and Vertical Cross-Section Helpers
#'
#' Convenience wrappers that automatically slice along time (`nc_hovmoeller()`)
#' or vertical levels (`nc_xsection()`) and render formatted Hovmoeller or
#' dual-axis pressure/height cross-section plots.
#'
#' @inheritParams plot.ncr
#' @param lon,lat,lev,time Dimension selectors passed to [nc_slice()].
#' @return Invisibly returns the plotted `"ncr_slice"` object.
#' @examples
#' f <- nc_example()
#' nc <- nc_open(f)
#' nc_hovmoeller(nc, "temp", lat = 0)
#' nc_xsection(nc, "temp", lon = 140)
#' nc_close(nc)
#' @export
nc_hovmoeller <- function(x,
                          var = NULL,
                          lon = missing_arg(),
                          lat = missing_arg(),
                          lev = missing_arg(),
                          time = NULL,
                          time_axis = "y",
                          filled = TRUE,
                          contour = TRUE,
                          ...) {
  if (inherits(x, "ncr")) {
    # If neither lon nor lat was fixed by caller or global state, default to lat = 0 (equator)
    st <- nc_get()
    lon_eff <- if (is_missing_arg(lon)) st$lon else lon
    lat_eff <- if (is_missing_arg(lat)) st$lat else lat
    if (is.null(lon_eff) && is.null(lat_eff)) {
      lat <- 0
    }
    s <- nc_slice(x, var = var, lon = lon, lat = lat, lev = lev, time = time)
  } else if (inherits(x, "ncr_slice")) {
    s <- x
  } else {
    stop("'x' must be an 'ncr' or 'ncr_slice' object.", call. = FALSE)
  }
  plot.ncr_slice(s, filled = filled, contour = contour, time_axis = time_axis, ...)
}

#' @rdname nc_hovmoeller
#' @export
nc_xsection <- function(x,
                        var = NULL,
                        lon = missing_arg(),
                        lat = missing_arg(),
                        lev = NULL,
                        time = missing_arg(),
                        zlog = TRUE,
                        filled = TRUE,
                        contour = TRUE,
                        ...) {
  if (inherits(x, "ncr")) {
    st <- nc_get()
    lon_eff <- if (is_missing_arg(lon)) st$lon else lon
    lat_eff <- if (is_missing_arg(lat)) st$lat else lat
    if (is.null(lon_eff) && is.null(lat_eff)) {
      lon <- stats::median(x$meta$dims[[x$meta$axis_map$lon]]$vals, na.rm = TRUE)
    }
    s <- nc_slice(x, var = var, lon = lon, lat = lat, lev = lev, time = time)
  } else if (inherits(x, "ncr_slice")) {
    s <- x
  } else {
    stop("'x' must be an 'ncr' or 'ncr_slice' object.", call. = FALSE)
  }
  plot.ncr_slice(s, filled = filled, contour = contour, zlog = zlog, ...)
}

# ---- Internal Plotting Helpers ----

.default_title <- function(x) {
  u <- if (nzchar(x$units)) sprintf(" [%s]", x$units) else ""
  base_title <- paste0(x$long_name, u)
  if (length(x$fixed_dims) > 0L) {
    fx <- vapply(x$fixed_dims, function(d) {
      if (d$axis == "time" && !is.null(d$time_posix)) {
        format(d$time_posix[1L], "%Y-%m-%d %H:%M")
      } else {
        paste0(d$name, "=", signif(d$vals[1L], 4),
               if (nzchar(d$units)) paste0(" ", d$units) else "")
      }
    }, character(1))
    paste0(base_title, "\n(", paste(fx, collapse = ", "), ")")
  } else {
    base_title
  }
}

.compute_levels_cols <- function(z, nlevels = 12, levels = NULL, col = NULL) {
  zrng <- range(z, finite = TRUE)
  if (!all(is.finite(zrng))) zrng <- c(0, 1)
  if (zrng[1L] == zrng[2L]) zrng <- zrng + c(-0.5, 0.5)
  if (is.null(levels)) {
    levels <- pretty(zrng, n = nlevels)
  }
  nbins <- max(1L, length(levels) - 1L)
  if (is.null(col)) {
    cols <- grDevices::hcl.colors(nbins, "RdYlBu", rev = TRUE)
  } else if (is.function(col)) {
    cols <- col(nbins)
  } else if (length(col) == 1L && !startsWith(col, "#") &&
             col %in% grDevices::hcl.pals()) {
    cols <- grDevices::hcl.colors(nbins, col, rev = TRUE)
  } else {
    cols <- col
  }
  list(levels = levels, cols = cols)
}

.pretty_deg_ticks <- function(lim, axis_type = c("lon", "lat"), target_n = 6L) {
  axis_type <- match.arg(axis_type)
  lim <- sort(range(lim, na.rm = TRUE))
  span <- diff(lim)
  if (!is.finite(span) || span <= 0) return(lim[1L])

  steps <- if (axis_type == "lon") {
    c(1, 2, 5, 10, 15, 30, 45, 60, 90)
  } else {
    c(1, 2, 5, 10, 15, 30, 45, 90)
  }

  counts <- vapply(steps, function(s) {
    lo <- ceiling((lim[1L] - 1e-9) / s) * s
    hi <- floor((lim[2L] + 1e-9) / s) * s
    if (hi < lo) return(0L)
    as.integer(round((hi - lo) / s)) + 1L
  }, integer(1))

  valid <- which(counts >= 3L & counts <= 9L)
  best <- if (length(valid) > 0L) {
    valid[which.min(abs(counts[valid] - target_n))]
  } else {
    which.min(abs(counts - target_n))
  }
  step <- steps[best]
  lo <- ceiling((lim[1L] - 1e-9) / step) * step
  hi <- floor((lim[2L] + 1e-9) / step) * step
  if (hi < lo) return(pretty(lim, n = target_n))
  seq(lo, hi, by = step)
}

.format_deg_axis <- function(side, at, axis_type = c("lon", "lat")) {
  axis_type <- match.arg(axis_type)
  labs <- vapply(at, function(v) {
    if (is.na(v)) return("")
    if (axis_type == "lon") {
      v_norm <- ((v + 180) %% 360) - 180
      if (abs(v_norm) < 1e-6 || abs( abs(v_norm) - 180 ) < 1e-6) {
        paste0(round(abs(v_norm), 1), "\u00b0")
      } else if (v_norm > 0) {
        paste0(round(v_norm, 1), "\u00b0E")
      } else {
        paste0(round(abs(v_norm), 1), "\u00b0W")
      }
    } else {
      if (abs(v) < 1e-6) {
        "0\u00b0"
      } else if (v > 0) {
        paste0(round(v, 1), "\u00b0N")
      } else {
        paste0(round(abs(v), 1), "\u00b0S")
      }
    }
  }, character(1))
  graphics::axis(side, at = at, labels = labs)
}

.overlay_world_map <- function(xlim, ylim, map = TRUE, proj = "lonlat",
                               lon0 = NULL, lat0 = NULL, lat_lim = NULL,
                               col = "gray25", lwd = 0.8) {
  if (is.null(map) || isFALSE(map)) return(invisible(NULL))

  # Custom list(x, y) polyline
  if (is.list(map) && !is.null(map$x) && !is.null(map$y)) {
    m <- list(x = map$x, y = map$y)
    if (proj == "lonlat") {
      graphics::lines(m$x, m$y, col = col, lwd = lwd)
    } else {
      if (!is.null(lat_lim)) {
        out <- !is.na(m$y) & (m$y < lat_lim[1L] | m$y > lat_lim[2L])
        m$x[out] <- NA_real_
        m$y[out] <- NA_real_
      }
      pm <- .project_map_lines(m, proj = proj, lon0 = lon0, lat0 = lat0)
      graphics::lines(pm$x, pm$y, col = col, lwd = lwd)
    }
    return(invisible(NULL))
  }

  dbs <- if (isTRUE(map)) "world" else as.character(map)
  if (length(dbs) == 1L && tolower(dbs[1L]) == "japan") {
    # mapdata::japan only contains Japan; layer worldHires (excluding Japan) for neighbors
    dbs <- c("worldHires", "japan")
  }
  has_japan_db <- any(tolower(dbs) %in% c("japan", "mapdata::japan"))

  mapdata_dbs <- c("worldHires", "world2Hires", "japan", "china", "rivers")
  wrap <- if (max(xlim, na.rm = TRUE) > 180.5) c(0, 360) else c(-180, 180)

  for (db in dbs) {
    if (tolower(db) %in% c("hires", "high", "worldhires")) db <- "worldHires"
    is_mapdata <- (db %in% mapdata_dbs) || startsWith(db, "mapdata::")
    if (is_mapdata) {
      if (!requireNamespace("mapdata", quietly = TRUE)) {
        stop(sprintf("Package 'mapdata' is required for map = '%s'. Install it with install.packages('mapdata').", db),
             call. = FALSE)
      }
      if (!startsWith(db, "mapdata::")) db <- paste0("mapdata::", db)
    }

    show_interior <- endsWith(db, "japan") || endsWith(db, "state") || endsWith(db, "county")
    is_world_db <- grepl("(world|worldHires|world2|world2Hires)$", db, ignore.case = TRUE)
    regions_arg <- "."
    exact_arg <- FALSE
    if (has_japan_db && is_world_db) {
      regs <- tryCatch({
        if (proj == "lonlat" && max(xlim, na.rm = TRUE) <= 180.5) {
          maps::map(db, xlim = xlim, ylim = ylim, namesonly = TRUE, plot = FALSE)
        } else {
          maps::map(db, namesonly = TRUE, plot = FALSE)
        }
      }, error = function(e) character(0))
      regs <- regs[!grepl("^Japan($|:)", regs, ignore.case = TRUE)]
      if (length(regs) == 0L) next
      regions_arg <- regs
      exact_arg <- TRUE
    }

    is_sub_border <- endsWith(db, "japan") && length(dbs) > 1L
    col_i <- col
    lwd_i <- if (is_sub_border) lwd * 0.75 else lwd

    if (proj == "lonlat") {
      tryCatch(
        maps::map(db, regions = regions_arg, exact = exact_arg,
                  xlim = xlim, ylim = ylim, wrap = wrap, add = TRUE,
                  col = col_i, lwd = lwd_i, interior = show_interior),
        error = function(e) NULL
      )
    } else {
      m <- tryCatch(
        maps::map(db, regions = regions_arg, exact = exact_arg,
                  plot = FALSE, fill = FALSE, interior = show_interior),
        error = function(e) NULL
      )
      if (!is.null(m)) {
        if (!is.null(lat_lim)) {
          out <- !is.na(m$y) & (m$y < lat_lim[1L] | m$y > lat_lim[2L])
          m$x[out] <- NA_real_
          m$y[out] <- NA_real_
        }
        pm <- .project_map_lines(m, proj = proj, lon0 = lon0, lat0 = lat0)
        graphics::lines(pm$x, pm$y, col = col_i, lwd = lwd_i)
      }
    }
  }
}

.plot_2d_slice <- function(x, filled, contour, map, proj, lon0, lat0,
                           nlevels, levels, col, skip, uvref, length,
                           add, main, sub = NULL, xlab, ylab, ...) {
  proj_norm <- tolower(trimws(as.character(proj)[1L]))
  is_proj <- (x$plot_type == "lonlat") &&
    !(proj_norm %in% c("lonlat", "cylindrical", "equirectangular", "eqc", "none", ""))

  if (is_proj) {
    .plot_projected_2d(
      x, filled = filled, contour = contour, map = map,
      proj = proj_norm, lon0 = lon0, lat0 = lat0,
      nlevels = nlevels, levels = levels, col = col,
      skip = skip, uvref = uvref, length = length,
      main = main, sub = sub, ...
    )
    return(invisible(x))
  }

  xv <- x$x
  yv <- x$y
  z  <- x$z
  lc <- .compute_levels_cols(z, nlevels = nlevels, levels = levels, col = col)
  if (is.null(main)) main <- .default_title(x)
  if (is.null(xlab)) xlab <- if (x$x_dim$axis == "lon") "Longitude" else x$x_dim$name
  if (is.null(ylab)) ylab <- if (x$y_dim$axis == "lat") "Latitude"  else x$y_dim$name

  if (isTRUE(add)) {
    if (isTRUE(filled)) {
      graphics::.filled.contour(xv, yv, z, levels = lc$levels, col = lc$cols)
    }
    if (isTRUE(contour) && !isTRUE(x$is_vector)) {
      graphics::contour(xv, yv, z, levels = lc$levels, add = TRUE,
                        col = if (filled) "gray20" else if (!is.null(col) && length(col) == 1L) col else "black",
                        ...)
    }
    if (isTRUE(x$has_uv) && !is.null(x$u) && !is.null(x$v)) {
      drawvec(xv, yv, x$u, x$v, length = length, uvref = uvref, skip = skip)
    }
    return(invisible(x))
  }

  old_mar <- graphics::par("mar")
  if (isTRUE(filled)) {
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    mar_new <- old_mar
    mar_new[4L] <- max(mar_new[4L], 5.5)
    graphics::par(mar = mar_new)
  }

  xlim <- range(xv, na.rm = TRUE)
  ylim <- range(yv, na.rm = TRUE)
  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                 main = main, sub = sub, xaxs = "i", yaxs = "i", axes = FALSE, ...)

  if (isTRUE(filled)) {
    graphics::.filled.contour(xv, yv, z, levels = lc$levels, col = lc$cols)
  }

  if (isTRUE(contour) && (!isTRUE(x$is_vector) || !isTRUE(filled))) {
    if (!isTRUE(x$is_vector)) {
      c_col <- if (isTRUE(filled)) "gray25" else if (!is.null(col) && length(col) == 1L) col else "black"
      graphics::contour(xv, yv, z, levels = lc$levels, add = TRUE, col = c_col,
                        lwd = if (filled) 0.6 else 1.0)
    }
  }

  if (!isFALSE(map) && !is.null(map) && x$plot_type == "lonlat") {
    .overlay_world_map(xlim, ylim, map = map, proj = "lonlat")
  }

  if (isTRUE(x$has_uv) && !is.null(x$u) && !is.null(x$v)) {
    uv_u <- if (!is.null(x$uv_units)) x$uv_units else x$units
    drawvec(xv, yv, x$u, x$v, length = length, uvref = uvref, skip = skip,
            ref_legend = TRUE, ref_units = uv_u)
  }

  if (x$x_dim$axis == "lon") {
    .format_deg_axis(1, .pretty_deg_ticks(xlim, "lon"), "lon")
  } else {
    graphics::axis(1)
  }
  if (x$y_dim$axis == "lat") {
    .format_deg_axis(2, .pretty_deg_ticks(ylim, "lat"), "lat")
  } else {
    graphics::axis(2)
  }
  graphics::box()

  if (isTRUE(filled)) {
    fields::image.plot(
      legend.only = TRUE,
      zlim        = range(lc$levels),
      breaks      = lc$levels,
      col         = lc$cols,
      legend.mar  = 4.2
    )
  }
}

.plot_projected_2d <- function(x, filled, contour, map, proj, lon0, lat0,
                               nlevels, levels, col, skip, uvref, length,
                               main, sub = NULL, ...) {
  lon <- x$x
  lat <- x$y
  z   <- x$z
  nx  <- length(lon)
  ny  <- length(lat)

  lc <- .compute_levels_cols(z, nlevels = nlevels, levels = levels, col = col)
  if (is.null(main)) main <- .default_title(x)

  # Cell edges for filled quad rendering
  dlon <- if (nx > 1L) diff(lon) else 1
  dlat <- if (ny > 1L) diff(lat) else 1
  lon_e <- c(lon[1L] - dlon[1L] / 2, lon[-nx] + dlon / 2, lon[nx] + dlon[nx - 1L] / 2)
  lat_e <- c(lat[1L] - dlat[1L] / 2, lat[-ny] + dlat / 2, lat[ny] + dlat[ny - 1L] / 2)
  lat_e <- pmax(pmin(lat_e, 89.9), -89.9)

  # Project grid cell corners
  LON_E <- rep(lon_e, times = ny + 1L)
  LAT_E <- rep(lat_e, each  = nx + 1L)
  pe <- nc_project(LON_E, LAT_E, proj = proj, lon0 = lon0, lat0 = lat0)
  lon0 <- pe$lon0
  lat0 <- pe$lat0
  PX <- matrix(pe$x, nrow = nx + 1L, ncol = ny + 1L)
  PY <- matrix(pe$y, nrow = nx + 1L, ncol = ny + 1L)

  is_polar <- (pe$proj == "stereographic") && (abs(lat0) >= 89)
  lat_lim <- NULL
  if (is_polar) {
    lat_bound <- if (lat0 >= 89) {
      max(min(lat, na.rm = TRUE), 0)
    } else {
      min(max(lat, na.rm = TRUE), 0)
    }
    lat_lim <- if (lat0 >= 89) c(lat_bound, 90) else c(-90, lat_bound)
    phi_b <- lat_bound * pi / 180
    r_bound <- 2 * cos(phi_b) / (1 + sin(abs(phi_b)))
    xlim <- c(-1.12, 1.12) * r_bound
    ylim <- c(-1.12, 1.12) * r_bound
  } else if (pe$proj == "orthographic") {
    xlim <- c(-1.05, 1.05)
    ylim <- c(-1.05, 1.05)
  } else {
    xlim <- if (any(is.finite(PX))) range(PX, finite = TRUE) else c(-1, 1)
    ylim <- if (any(is.finite(PY))) range(PY, finite = TRUE) else c(-1, 1)
  }

  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  has_sub <- !is.null(sub) && nzchar(sub)
  mar_bot <- if (isTRUE(x$has_uv) || has_sub) 3.8 else 1.4
  mar_top <- if (grepl("\n", main)) 3.6 else 2.8
  mar_new <- c(
    mar_bot,
    1.4,
    mar_top,
    if (isTRUE(filled)) 5.2 else 1.4
  )
  graphics::par(mar = mar_new)

  graphics::plot(NA, xlim = xlim, ylim = ylim, asp = 1,
                 xlab = "", ylab = "", main = main, sub = sub, axes = FALSE, ...)

  if (isTRUE(filled)) {
    # Vectorized polygon construction for projected quads
    i_idx <- rep(seq_len(nx), times = ny)
    j_idx <- rep(seq_len(ny), each  = nx)
    x1 <- PX[cbind(i_idx,     j_idx)]
    x2 <- PX[cbind(i_idx + 1L, j_idx)]
    x3 <- PX[cbind(i_idx + 1L, j_idx + 1L)]
    x4 <- PX[cbind(i_idx,     j_idx + 1L)]

    y1 <- PY[cbind(i_idx,     j_idx)]
    y2 <- PY[cbind(i_idx + 1L, j_idx)]
    y3 <- PY[cbind(i_idx + 1L, j_idx + 1L)]
    y4 <- PY[cbind(i_idx,     j_idx + 1L)]

    z_vec <- as.vector(z)
    c_idx <- findInterval(z_vec, lc$levels, all.inside = TRUE)
    cell_cols <- lc$cols[c_idx]

    # Filter out quads that cross the projection seam or contain NA vertices
    span_x <- diff(xlim)
    span_y <- diff(ylim)
    max_dx <- pmax(abs(x2 - x1), abs(x3 - x2), abs(x4 - x3), abs(x1 - x4))
    max_dy <- pmax(abs(y2 - y1), abs(y3 - y2), abs(y4 - y3), abs(y1 - y4))
    valid <- !is.na(z_vec) & !is.na(max_dx) & !is.na(max_dy) &
      (max_dx < 0.25 * span_x) & (max_dy < 0.25 * span_y)

    if (any(valid)) {
      vx <- rbind(x1[valid], x2[valid], x3[valid], x4[valid], NA_real_)
      vy <- rbind(y1[valid], y2[valid], y3[valid], y4[valid], NA_real_)
      graphics::polygon(as.vector(vx), as.vector(vy),
                        col = cell_cols[valid], border = NA)
    }
  }

  if (isTRUE(contour) && !isTRUE(x$is_vector)) {
    clines <- grDevices::contourLines(lon, lat, z, levels = lc$levels)
    c_col <- if (isTRUE(filled)) "gray25" else if (!is.null(col) && length(col) == 1L) col else "black"
    c_lwd <- if (isTRUE(filled)) 0.6 else 1.0
    for (cl in clines) {
      pc <- .project_map_lines(list(x = cl$x, y = cl$y),
                               proj = pe$proj, lon0 = lon0, lat0 = lat0)
      if (any(!is.na(pc$x)) && any(!is.na(pc$y))) {
        graphics::lines(pc$x, pc$y, col = c_col, lwd = c_lwd)
      }
    }
  }

  .draw_graticule(pe$proj, lon0 = lon0, lat0 = lat0, lat_lim = lat_lim)

  if (!isFALSE(map) && !is.null(map)) {
    .overlay_world_map(range(lon), range(lat), map = map, proj = pe$proj,
                       lon0 = lon0, lat0 = lat0, lat_lim = lat_lim)
  }

  if (isTRUE(x$has_uv) && !is.null(x$u) && !is.null(x$v)) {
    uv_u <- if (!is.null(x$uv_units)) x$uv_units else x$units
    .draw_projected_vectors(lon, lat, x$u, x$v, proj = pe$proj,
                            lon0 = lon0, lat0 = lat0,
                            skip = skip, uvref = uvref, length = length,
                            ref_legend = TRUE, ref_units = uv_u)
  }

  if (isTRUE(filled)) {
    y_nfc <- graphics::grconvertY(ylim, from = "user", to = "nfc")
    plt   <- graphics::par("plt")
    sp <- c(
      min(plt[2L] + 0.02, 0.90),
      min(plt[2L] + 0.055, 0.94),
      max(y_nfc[1L], 0.02),
      min(y_nfc[2L], 0.98)
    )
    if (sp[4L] > sp[3L] + 0.05 && sp[2L] > sp[1L] + 0.01) {
      fields::image.plot(
        legend.only = TRUE,
        zlim        = range(lc$levels),
        breaks      = lc$levels,
        col         = lc$cols,
        smallplot   = sp
      )
    } else {
      fields::image.plot(
        legend.only = TRUE,
        zlim        = range(lc$levels),
        breaks      = lc$levels,
        col         = lc$cols,
        legend.mar  = 4.2
      )
    }
  }
}

.draw_projected_vectors <- function(lon, lat, u, v, proj, lon0, lat0,
                                    skip = 1L, uvref = NA, length = 0.05,
                                    col = "black", ref_legend = FALSE,
                                    ref_units = "", ...) {
  skip <- max(1L, as.integer(skip)[1L])
  iy <- seq.int(1L, length(lat), by = skip)
  # Exclude exact poles where zonal tangent direction is singular
  iy <- iy[abs(lat[iy]) <= 89.5]
  if (length(iy) == 0L) return(invisible(NA_real_))

  # Latitude-adaptive longitude thinning so arrows do not crowd near the poles
  LON <- numeric(0)
  LAT <- numeric(0)
  u_s <- numeric(0)
  v_s <- numeric(0)
  for (j in iy) {
    cos_j <- max(cos(lat[j] * pi / 180), 0.08)
    step_lon <- max(skip, as.integer(round(skip / cos_j)))
    ix_j <- seq.int(1L, length(lon), by = step_lon)
    LON <- c(LON, lon[ix_j])
    LAT <- c(LAT, rep(lat[j], length(ix_j)))
    u_s <- c(u_s, u[ix_j, j])
    v_s <- c(v_s, v[ix_j, j])
  }

  p0 <- nc_project(LON, LAT, proj = proj, lon0 = lon0, lat0 = lat0)
  # Small displacement in degrees to compute local projected tangent orientation
  eps <- 0.25
  cos_lat <- pmax(cos(LAT * pi / 180), 0.05)
  p1 <- nc_project(LON + eps * u_s / cos_lat, LAT + eps * v_s,
                   proj = proj, lon0 = lon0, lat0 = lat0)
  dx <- p1$x - p0$x
  dy <- p1$y - p0$y
  dir_norm <- sqrt(dx^2 + dy^2)
  spd <- sqrt(u_s^2 + v_s^2)
  if (is.na(uvref)) uvref <- max(spd, na.rm = TRUE)
  if (!is.finite(uvref) || uvref <= 0) return(invisible(uvref))

  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  # Base arrow scale on visible latitude rows across the plot height
  n_vis_y <- max(length(unique(LAT[!is.na(p0$x) & !is.na(p0$y)])), 6L)
  cell_d <- min(abs(usr[2L] - usr[1L]), abs(usr[4L] - usr[3L])) / n_vis_y
  scale_len <- 0.55 * cell_d * (spd / uvref)

  ux_all <- (dx / dir_norm) * scale_len
  vy_all <- (dy / dir_norm) * scale_len
  dx_in <- (2 * ux_all) * pin[1L] / max(abs(usr[2L] - usr[1L]), 1e-12)
  dy_in <- (2 * vy_all) * pin[2L] / max(abs(usr[4L] - usr[3L]), 1e-12)
  len_in <- sqrt(dx_in^2 + dy_in^2)

  keep <- !is.na(p0$x) & !is.na(p0$y) & !is.na(dir_norm) &
    (dir_norm > 1e-9) & (spd > 1e-9) & is.finite(len_in) & (len_in >= 1e-3)
  if (any(keep)) {
    graphics::arrows(p0$x[keep] - ux_all[keep], p0$y[keep] - vy_all[keep],
                     p0$x[keep] + ux_all[keep], p0$y[keep] + vy_all[keep],
                     length = length, col = col, ...)
  }
  if (isTRUE(ref_legend)) {
    .draw_ref_vector(uvref, ref_len_usr = 1.1 * cell_d, length = length,
                     units = ref_units, col = col)
  }
  invisible(uvref)
}

.plot_hovmoeller_slice <- function(x, filled, contour, nlevels, levels, col,
                                   time_axis = "y", skip, uvref, length,
                                   main, sub = NULL, xlab, ylab, ...) {
  time_on_y <- tolower(time_axis)[1L] != "x"
  # In .slice_single_var, active_dims[[1]] is spatial (lon/lat) and active_dims[[2]] is time
  s_dim <- x$x_dim
  t_dim <- x$y_dim
  s_vals <- x$x
  t_vals <- if (!is.null(t_dim$time_posix)) as.numeric(t_dim$time_posix) else x$y
  z <- x$z

  if (!time_on_y) {
    xv <- t_vals
    yv <- s_vals
    z  <- t(z)
    if (is.null(xlab)) xlab <- "Time"
    if (is.null(ylab)) ylab <- if (s_dim$axis == "lon") "Longitude" else "Latitude"
  } else {
    xv <- s_vals
    yv <- t_vals
    if (is.null(xlab)) xlab <- if (s_dim$axis == "lon") "Longitude" else "Latitude"
    if (is.null(ylab)) ylab <- "Time"
  }

  lc <- .compute_levels_cols(z, nlevels = nlevels, levels = levels, col = col)
  if (is.null(main)) main <- paste("Hovmoeller:", .default_title(x))

  old_mar <- graphics::par("mar")
  mar_new <- old_mar
  need_mar <- FALSE
  if (isTRUE(filled)) {
    mar_new[4L] <- max(mar_new[4L], 5.5)
    need_mar <- TRUE
  }
  ylab_plot <- ylab
  ylab_line <- NA_real_
  if (time_on_y && !is.null(t_dim$time_posix)) {
    t_info <- .time_axis_info(t_dim$time_posix)
    w_lines <- if (grDevices::dev.cur() > 1L) {
      max(graphics::strwidth(t_info$labels, units = "inches"), na.rm = TRUE) /
        max(graphics::par("csi"), 1e-3)
    } else {
      max(nchar(t_info$labels), 4L) * 0.48
    }
    left_needed <- max(mar_new[2L], ceiling(w_lines) + 2.4)
    if (left_needed > mar_new[2L]) {
      mar_new[2L] <- left_needed
      need_mar <- TRUE
    }
    ylab_line <- mar_new[2L] - 1.3
    ylab_plot <- ""
  }
  if (need_mar) {
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    graphics::par(mar = mar_new)
  }

  xlim <- range(xv, na.rm = TRUE)
  ylim <- range(yv, na.rm = TRUE)
  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab_plot,
                 main = main, sub = sub, xaxs = "i", yaxs = "i", axes = FALSE, ...)
  if (!is.na(ylab_line) && nzchar(ylab)) {
    graphics::mtext(ylab, side = 2, line = ylab_line)
  }

  if (isTRUE(filled)) {
    graphics::.filled.contour(xv, yv, z, levels = lc$levels, col = lc$cols)
  }
  if (isTRUE(contour)) {
    graphics::contour(xv, yv, z, levels = lc$levels, add = TRUE,
                      col = if (filled) "gray25" else "black",
                      lwd = if (filled) 0.6 else 1.0)
  }

  # Format axes cleanly
  if (time_on_y) {
    if (s_dim$axis %in% c("lon", "lat")) {
      .format_deg_axis(1, .pretty_deg_ticks(xlim, s_dim$axis), s_dim$axis)
    } else {
      graphics::axis(1)
    }
    if (!is.null(t_dim$time_posix)) {
      .draw_nice_time_axis(2, t_dim$time_posix)
    } else {
      graphics::axis(2)
    }
  } else {
    if (!is.null(t_dim$time_posix)) {
      .draw_nice_time_axis(1, t_dim$time_posix)
    } else {
      graphics::axis(1)
    }
    if (s_dim$axis %in% c("lon", "lat")) {
      .format_deg_axis(2, .pretty_deg_ticks(ylim, s_dim$axis), s_dim$axis)
    } else {
      graphics::axis(2)
    }
  }
  graphics::box()

  if (isTRUE(filled)) {
    fields::image.plot(
      legend.only = TRUE,
      zlim        = range(lc$levels),
      breaks      = lc$levels,
      col         = lc$cols,
      legend.mar  = 4.2
    )
  }
}

.time_axis_info <- function(t_posix) {
  span_sec <- diff(range(as.numeric(t_posix), na.rm = TRUE))
  fmt <- if (span_sec <= 3 * 86400) {
    "%b %d %H:%M"
  } else if (span_sec <= 180 * 86400) {
    "%b %d"
  } else if (span_sec <= 3 * 365 * 86400) {
    "%Y-%m"
  } else {
    "%Y"
  }
  at_t <- pretty(t_posix, n = 6)
  list(at = at_t, fmt = fmt, labels = format(at_t, fmt))
}

.draw_nice_time_axis <- function(side, t_posix) {
  info <- .time_axis_info(t_posix)
  graphics::axis.POSIXct(side, at = info$at, format = info$fmt,
                         las = if (side == 2) 1 else 0)
}

.is_pressure_dim <- function(d) {
  un <- tolower(d$units)
  pos <- tolower(d$positive)
  if (pos == "down" || grepl("^(hpa|pa|mb|millibar|mbar|dbar)$", un) ||
      grepl("pressure|plev|isobaric", tolower(d$name))) {
    return(TRUE)
  }
  if (pos == "up" || grepl("^(m|km|meter|meters|kilometers)$", un) ||
      grepl("height|alt|altitude", tolower(d$name))) {
    return(FALSE)
  }
  # Fallback heuristic: standard atmospheric pressure levels in 10..1100
  isTRUE(d$reversed) || (max(d$vals, na.rm = TRUE) >= 500 && min(d$vals, na.rm = TRUE) >= 1)
}

.plot_xsection_slice <- function(x, filled, contour, nlevels, levels, col,
                                 zlog = TRUE, skip, uvref, length,
                                 main, sub = NULL, xlab, ylab, ...) {
  s_dim <- x$x_dim
  v_dim <- x$y_dim
  xv <- x$x
  lev <- x$y
  z <- x$z

  is_pres <- .is_pressure_dim(v_dim)
  # Convert Pa to hPa for display if needed
  pres_scale <- if (tolower(v_dim$units) == "pa" && max(lev, na.rm = TRUE) > 2000) 0.01 else 1.0
  lev_disp <- lev * pres_scale

  use_log <- is_pres && isTRUE(zlog) && all(lev_disp > 0)
  if (use_log) {
    # y_coord = -log(p) is strictly increasing when p is sorted decreasingly!
    # Since lev_disp is currently sorted in ascending order (e.g. 100..1000),
    # -log(lev_disp) is decreasing; reverse it so y_coord is strictly increasing
    # from bottom (1000 hPa, -log(1000)) to top (100 hPa, -log(100))!
    ord <- rev(seq_along(lev_disp))
    lev_disp <- lev_disp[ord]
    z <- z[, ord, drop = FALSE]
    if (isTRUE(x$has_uv) && !is.null(x$u) && !is.null(x$v)) {
      x$u <- x$u[, ord, drop = FALSE]
      x$v <- x$v[, ord, drop = FALSE]
    }
    yv <- -log(lev_disp)
    ylim <- range(yv, na.rm = TRUE) # bottom = -log(max_p), top = -log(min_p)
  } else if (is_pres) {
    yv <- lev_disp
    ylim <- rev(range(yv, na.rm = TRUE)) # bottom = max_p, top = min_p
  } else {
    # Height coordinate (m or km)
    h_scale <- if (tolower(v_dim$units) %in% c("m", "meter", "meters") && max(lev, na.rm = TRUE) > 500) 0.001 else 1.0
    yv <- lev * h_scale
    ylim <- range(yv, na.rm = TRUE)
  }

  lc <- .compute_levels_cols(z, nlevels = nlevels, levels = levels, col = col)
  if (is.null(main)) main <- paste("Cross-Section:", .default_title(x))
  if (is.null(xlab)) xlab <- if (s_dim$axis == "lon") "Longitude" else if (s_dim$axis == "lat") "Latitude" else s_dim$name
  if (is.null(ylab)) ylab <- if (is_pres) "Pressure [hPa]" else "Height [km]"

  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar_new <- old_mar
  # Extra right margin for both secondary height/pressure axis and colorbar
  mar_new[4L] <- max(mar_new[4L], if (isTRUE(filled)) 7.8 else 4.8)
  graphics::par(mar = mar_new)

  xlim <- range(xv, na.rm = TRUE)
  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                 main = main, xaxs = "i", yaxs = "i", axes = FALSE)

  if (isTRUE(filled)) {
    graphics::.filled.contour(xv, yv, z, levels = lc$levels, col = lc$cols)
  }
  if (isTRUE(contour)) {
    graphics::contour(xv, yv, z, levels = lc$levels, add = TRUE,
                      col = if (filled) "gray25" else "black",
                      lwd = if (filled) 0.6 else 1.0)
  }

  if (isTRUE(x$has_uv) && !is.null(x$u) && !is.null(x$v)) {
    drawvec(xv, yv, x$u, x$v, length = length, uvref = uvref, skip = skip)
  }

  # Horizontal axis
  if (s_dim$axis %in% c("lon", "lat")) {
    .format_deg_axis(1, .pretty_deg_ticks(xlim, s_dim$axis), s_dim$axis)
  } else {
    graphics::axis(1)
  }

  # Dual Vertical Axes: Left = Pressure (hPa) & Right = Height (km) (or vice-versa)
  if (is_pres) {
    std_p <- c(1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10)
    prng <- range(lev_disp, na.rm = TRUE)
    p_ticks <- std_p[std_p >= prng[1L] * 0.98 & std_p <= prng[2L] * 1.02]
    if (length(p_ticks) < 3L) p_ticks <- pretty(prng, n = 6)
    p_ticks <- p_ticks[p_ticks > 0]
    y_at_p <- if (use_log) -log(p_ticks) else p_ticks
    graphics::axis(2, at = y_at_p, labels = format(p_ticks, trim = TRUE), las = 1)

    # Right axis: Standard-atmosphere Height [km]
    z_rng <- pres2alt(prng)
    z_ticks <- pretty(z_rng, n = 6)
    p_of_z <- alt2pres(z_ticks)
    y_at_z <- if (use_log) -log(p_of_z) else p_of_z
    graphics::axis(4, at = y_at_z, labels = format(z_ticks, trim = TRUE), las = 1)
    graphics::mtext("Height [km]", side = 4, line = 2.6, cex = 0.9)
  } else {
    # Primary coordinate is Height [km]
    z_ticks <- pretty(ylim, n = 6)
    graphics::axis(2, at = z_ticks, las = 1)

    # Secondary axis: Pressure [hPa]
    std_p <- c(1000, 850, 700, 500, 300, 200, 100, 50, 20, 10)
    z_of_p <- pres2alt(std_p)
    keep_p <- z_of_p >= ylim[1L] & z_of_p <= ylim[2L]
    if (any(keep_p)) {
      graphics::axis(4, at = z_of_p[keep_p], labels = std_p[keep_p], las = 1)
      graphics::mtext("Pressure [hPa]", side = 4, line = 2.6, cex = 0.9)
    }
  }
  graphics::box()

  if (isTRUE(filled)) {
    fields::image.plot(
      legend.only = TRUE,
      zlim        = range(lc$levels),
      breaks      = lc$levels,
      col         = lc$cols,
      legend.mar  = 3.8
    )
  }
}

.plot_1d_slice <- function(x, main, sub = NULL, xlab, ylab, zlog = TRUE, add = FALSE, ...) {
  d <- x$x_dim
  vals <- as.vector(x$z)
  if (is.null(main)) main <- .default_title(x)
  u_str <- if (nzchar(x$units)) sprintf(" [%s]", x$units) else ""

  if (d$axis == "lev") {
    is_pres <- .is_pressure_dim(d)
    lev <- d$vals
    use_log <- is_pres && isTRUE(zlog) && all(lev > 0)
    yv <- if (use_log) -log(lev) else lev
    ylim <- if (use_log) range(yv, na.rm = TRUE) else if (is_pres) rev(range(yv, na.rm = TRUE)) else range(yv, na.rm = TRUE)
    if (is.null(xlab)) xlab <- paste0(x$var, u_str)
    if (is.null(ylab)) ylab <- if (is_pres) "Pressure [hPa]" else "Height [km]"
    if (isTRUE(add)) {
      graphics::lines(vals, yv, ...)
    } else {
      graphics::plot(vals, yv, type = "b", pch = 19, ylim = ylim,
                     xlab = xlab, ylab = ylab, main = main, axes = FALSE, ...)
      graphics::axis(1)
      if (is_pres) {
        p_ticks <- pretty(lev, n = 6)
        p_ticks <- p_ticks[p_ticks > 0]
        graphics::axis(2, at = if (use_log) -log(p_ticks) else p_ticks,
                       labels = p_ticks, las = 1)
      } else {
        graphics::axis(2, las = 1)
      }
      graphics::box()
      graphics::grid()
    }
  } else if (d$axis == "time" && !is.null(d$time_posix)) {
    tv <- d$time_posix
    if (is.null(xlab)) xlab <- "Time"
    if (is.null(ylab)) ylab <- paste0(x$var, u_str)
    if (isTRUE(add)) {
      graphics::lines(tv, vals, ...)
    } else {
      graphics::plot(tv, vals, type = "l", lwd = 1.5, xlab = xlab, ylab = ylab,
                     main = main, xaxt = "n", ...)
      .draw_nice_time_axis(1, tv)
      graphics::grid()
    }
  } else {
    if (is.null(xlab)) xlab <- d$name
    if (is.null(ylab)) ylab <- paste0(x$var, u_str)
    if (isTRUE(add)) {
      graphics::lines(d$vals, vals, ...)
    } else {
      graphics::plot(d$vals, vals, type = "l", lwd = 1.5,
                     xlab = xlab, ylab = ylab, main = main, ...)
      graphics::grid()
    }
  }
}
