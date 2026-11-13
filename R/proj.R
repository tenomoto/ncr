#' Project Geographic Coordinates to Map Projection Coordinates
#'
#' Transforms longitude and latitude (in degrees) into projected `(x, y)`
#' coordinates in base R without external C projection libraries, avoiding binary
#' symbol conflicts with `RNetCDF`.
#'
#' @param lon,lat Numeric vectors of longitude and latitude in degrees.
#' @param proj Character string specifying the projection:
#'   * `"lonlat"` (or `"eqc"`): Equidistant cylindrical (unprojected lon/lat).
#'   * `"mercator"`: Mercator projection.
#'   * `"robinson"`: Robinson pseudocylindrical projection.
#'   * `"mollweide"`: Mollweide equal-area pseudocylindrical projection.
#'   * `"orthographic"` (or `"ortho"`): Orthographic globe projection.
#'   * `"stereographic"` (or `"stere"`): Azimuthal stereographic projection.
#'   * `"npstere"`: North Polar stereographic (`lat0 = 90`).
#'   * `"spstere"`: South Polar stereographic (`lat0 = -90`).
#'   * `"lambert"` (or `"lcc"`): Lambert Conformal Conic projection.
#' @param lon0,lat0 Numeric center longitude and latitude in degrees. If `NULL`,
#'   sensible defaults are chosen from the data domain.
#' @return A list with numeric vectors `x` and `y` (invisible points outside the
#'   projection horizon are set to `NA`), plus resolved `proj`, `lon0`, and
#'   `lat0`.
#' @examples
#' pts <- nc_project(c(120, 140, 160), c(20, 35, 50), proj = "orthographic",
#'                   lon0 = 140, lat0 = 35)
#' str(pts)
#' @export
nc_project <- function(lon, lat, proj = "lonlat", lon0 = NULL, lat0 = NULL) {
  proj <- tolower(trimws(as.character(proj)[1L]))
  if (proj %in% c("cylindrical", "equirectangular", "eqc", "none", "")) {
    proj <- "lonlat"
  } else if (proj == "ortho") {
    proj <- "orthographic"
  } else if (proj == "stere") {
    proj <- "stereographic"
  } else if (proj == "lcc") {
    proj <- "lambert"
  }

  if (proj == "npstere") {
    proj <- "stereographic"
    if (is.null(lat0)) lat0 <- 90
  } else if (proj == "spstere") {
    proj <- "stereographic"
    if (is.null(lat0)) lat0 <- -90
  }

  if (is.null(lon0)) {
    lon0 <- if (any(!is.na(lon))) mean(range(lon, na.rm = TRUE)) else 0
  }
  if (is.null(lat0)) {
    lat0 <- if (proj == "stereographic") {
      if (any(!is.na(lat)) && mean(lat, na.rm = TRUE) < 0) -90 else 90
    } else if (proj == "lambert") {
      if (any(!is.na(lat))) {
        m <- mean(range(lat, na.rm = TRUE))
        if (abs(m) < 10) 35 else m
      } else {
        35
      }
    } else if (proj == "orthographic") {
      if (any(!is.na(lat))) mean(range(lat, na.rm = TRUE)) else 20
    } else {
      0
    }
  }

  n <- max(length(lon), length(lat))
  if (length(lon) != n) lon <- rep_len(lon, n)
  if (length(lat) != n) lat <- rep_len(lat, n)

  if (proj == "lonlat") {
    return(list(x = lon, y = lat, proj = proj, lon0 = lon0, lat0 = lat0))
  }

  deg2rad <- pi / 180
  # Wrap (lon - lon0) to [-180, 180]
  dlon_deg <- ((lon - lon0 + 180) %% 360) - 180
  lam <- dlon_deg * deg2rad
  phi <- pmax(pmin(lat, 89.999), -89.999) * deg2rad
  phi0 <- lat0 * deg2rad

  if (proj == "mercator") {
    phi_m <- pmax(pmin(lat, 85), -85) * deg2rad
    x <- lam
    y <- log(tan(pi / 4 + phi_m / 2))
  } else if (proj == "robinson") {
    # Standard Robinson projection lookup table (0 to 90 deg by 5 deg)
    lat_tab <- seq(0, 90, by = 5)
    xl_tab  <- c(1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
                 0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
                 0.7186, 0.6732, 0.6213, 0.5722, 0.5322)
    yl_tab  <- c(0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
                 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
                 0.8435, 0.8936, 0.9394, 0.9761, 1.0000)
    abs_lat <- pmin(abs(lat), 90)
    xl <- stats::approx(lat_tab, xl_tab, xout = abs_lat, rule = 2)$y
    yl <- stats::approx(lat_tab, yl_tab, xout = abs_lat, rule = 2)$y
    x <- 0.8487 * xl * lam
    y <- 1.3523 * sign(lat) * yl
  } else if (proj == "mollweide") {
    theta <- phi
    for (iter in 1:6) {
      denom <- 2 + 2 * cos(2 * theta)
      denom[abs(denom) < 1e-12] <- 1e-12
      theta <- theta - (2 * theta + sin(2 * theta) - pi * sin(phi)) / denom
    }
    x <- (2 * sqrt(2) / pi) * lam * cos(theta)
    y <- sqrt(2) * sin(theta)
  } else if (proj == "orthographic") {
    cos_c <- sin(phi0) * sin(phi) + cos(phi0) * cos(phi) * cos(lam)
    vis <- !is.na(cos_c) & (cos_c >= -1e-6)
    x <- ifelse(vis, cos(phi) * sin(lam), NA_real_)
    y <- ifelse(vis, cos(phi0) * sin(phi) - sin(phi0) * cos(phi) * cos(lam), NA_real_)
  } else if (proj == "stereographic") {
    cos_c <- sin(phi0) * sin(phi) + cos(phi0) * cos(phi) * cos(lam)
    # Clip points within 20 deg of the antipode to prevent blowup
    vis <- !is.na(cos_c) & (cos_c > -0.85)
    k <- ifelse(vis, 2 / (1 + cos_c), NA_real_)
    x <- k * cos(phi) * sin(lam)
    y <- k * (cos(phi0) * sin(phi) - sin(phi0) * cos(phi) * cos(lam))
  } else if (proj == "lambert") {
    # Standard parallels +/- 10 deg around lat0
    phi1 <- max(min(lat0 - 10, 80), -80) * deg2rad
    phi2 <- max(min(lat0 + 10, 85), -85) * deg2rad
    if (abs(phi1 + phi2) < 1e-4) {
      phi1 <- 20 * deg2rad
      phi2 <- 50 * deg2rad
      if (abs(phi0) < 1e-4) phi0 <- 35 * deg2rad
    }
    if (abs(phi1 - phi2) < 1e-6) {
      n_cone <- sin(phi1)
    } else {
      n_cone <- log(cos(phi1) / cos(phi2)) /
        log(tan(pi / 4 + phi2 / 2) / tan(pi / 4 + phi1 / 2))
    }
    # Clip hemisphere opposite to cone apex
    vis <- if (n_cone >= 0) lat > -75 else lat < 75
    phi_c <- ifelse(vis, phi, NA_real_)
    F_const <- (cos(phi1) * (tan(pi / 4 + phi1 / 2))^n_cone) / n_cone
    rho  <- F_const / (tan(pi / 4 + phi_c / 2))^n_cone
    rho0 <- F_const / (tan(pi / 4 + phi0 / 2))^n_cone
    theta <- n_cone * lam
    x <- rho * sin(theta)
    y <- rho0 - rho * cos(theta)
  } else {
    stop(sprintf("Unsupported map projection '%s'.", proj), call. = FALSE)
  }

  list(x = x, y = y, proj = proj, lon0 = lon0, lat0 = lat0)
}

# Break map polyline segments that jump across the antimeridian / projection seam
.project_map_lines <- function(m, proj, lon0, lat0) {
  p <- nc_project(m$x, m$y, proj = proj, lon0 = lon0, lat0 = lat0)
  px <- p$x
  py <- p$y
  n <- length(px)
  if (n >= 2L) {
    dx <- abs(diff(px))
    dy <- abs(diff(py))
    span_x <- diff(range(px, na.rm = TRUE))
    span_y <- diff(range(py, na.rm = TRUE))
    thresh_x <- max(span_x * 0.22, 1e-6)
    thresh_y <- max(span_y * 0.22, 1e-6)
    jumps <- which((!is.na(dx) & dx > thresh_x) | (!is.na(dy) & dy > thresh_y))
    if (length(jumps) > 0L) {
      # Insert NA at segment breaks by masking the second endpoint of jumping segments
      px[jumps + 1L] <- NA_real_
      py[jumps + 1L] <- NA_real_
    }
  }
  list(x = px, y = py)
}

.draw_graticule <- function(proj, lon0, lat0, col = "gray75", lty = 3) {
  lons <- seq(-180, 150, by = 30)
  lats <- seq(-60, 60, by = 30)

  lat_fine <- seq(-85, 85, length.out = 90)
  for (ln in lons) {
    p <- .project_map_lines(list(x = rep(ln, length(lat_fine)), y = lat_fine),
                            proj = proj, lon0 = lon0, lat0 = lat0)
    graphics::lines(p$x, p$y, col = col, lty = lty, lwd = 0.6)
  }

  lon_fine <- seq(lon0 - 179.5, lon0 + 179.5, length.out = 180)
  for (lt in lats) {
    p <- .project_map_lines(list(x = lon_fine, y = rep(lt, length(lon_fine))),
                            proj = proj, lon0 = lon0, lat0 = lat0)
    graphics::lines(p$x, p$y, col = col, lty = lty, lwd = 0.6)
  }

  if (proj == "orthographic") {
    ang <- seq(0, 2 * pi, length.out = 180)
    graphics::lines(cos(ang), sin(ang), col = "gray40", lwd = 1)
  }
}
