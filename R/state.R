# Global GrADS-like state environment for ncr
.ncr_state <- new.env(parent = emptyenv())

.init_ncr_state <- function() {
  .ncr_state$lon  <- NULL
  .ncr_state$lat  <- NULL
  .ncr_state$lev  <- I(1L)
  .ncr_state$time <- 1L
  .ncr_state$proj <- "lonlat"
  .ncr_state$lon0 <- NULL
  .ncr_state$lat0 <- NULL
  .ncr_state$map  <- TRUE
}

.init_ncr_state()

#' Get, Set, or Reset Global Slicing State
#'
#' Manages the GrADS-inspired global state (`lon`, `lat`, `lev`, `time`, `proj`)
#' used for automatic slicing of NetCDF variables when dimensions are not
#' explicitly specified in `[` or [nc_slice()].
#'
#' @param lon Numeric value (nearest grid point), length-2 numeric vector
#'   `c(min, max)` (coordinate range), `I(idx)` (1-based indices), or `NULL`
#'   (full range).
#' @param lat Numeric value, length-2 numeric vector `c(min, max)`, `I(idx)`,
#'   or `NULL` (full range).
#' @param lev Numeric coordinate value (e.g., `500` for 500 hPa), length-2
#'   numeric vector `c(min, max)` (or `c(1000, 100)`), `I(idx)` for 1-based
#'   index, or `NULL` for all levels. Defaults to `I(1L)` (first level).
#' @param time Integer index (default `1L`), character timestamp or length-2
#'   range (e.g., `"2020-01-01"` or `c("2020-01-01", "2020-01-10")`), `POSIXt`
#'   / `Date`, or `NULL` for all time steps.
#' @param proj Character string naming the map projection (`"lonlat"`,
#'   `"mercator"`, `"robinson"`, `"mollweide"`, `"orthographic"`,
#'   `"stereographic"`, `"npstere"`, `"spstere"`, `"lambert"`).
#' @param lon0,lat0 Optional projection center longitude and latitude.
#' @param map Default map overlay (`TRUE` for `"world"`, `FALSE` for none,
#'   `"hires"` / `"worldHires"` / `"japan"` via `mapdata`, or a custom
#'   `list(x, y)`).
#' @return A named list of the current global state settings (`nc_set` and
#'   `nc_reset` return the updated state invisibly; `nc_get` returns it
#'   visibly).
#' @examples
#' nc_get()
#' nc_set(lon = c(100, 180), lat = c(0, 60), lev = 500, time = 1)
#' nc_reset()
#' @export
nc_set <- function(lon = missing_arg(),
                   lat = missing_arg(),
                   lev = missing_arg(),
                   time = missing_arg(),
                   proj = missing_arg(),
                   lon0 = missing_arg(),
                   lat0 = missing_arg(),
                   map = missing_arg()) {
  if (!is_missing_arg(lon))  .ncr_state$lon  <- lon
  if (!is_missing_arg(lat))  .ncr_state$lat  <- lat
  if (!is_missing_arg(lev))  .ncr_state$lev  <- lev
  if (!is_missing_arg(time)) .ncr_state$time <- time
  if (!is_missing_arg(proj)) .ncr_state$proj <- proj
  if (!is_missing_arg(lon0)) .ncr_state$lon0 <- lon0
  if (!is_missing_arg(lat0)) .ncr_state$lat0 <- lat0
  if (!is_missing_arg(map))  .ncr_state$map  <- map
  invisible(nc_get())
}

#' @rdname nc_set
#' @export
nc_get <- function() {
  list(
    lon  = .ncr_state$lon,
    lat  = .ncr_state$lat,
    lev  = .ncr_state$lev,
    time = .ncr_state$time,
    proj = .ncr_state$proj,
    lon0 = .ncr_state$lon0,
    lat0 = .ncr_state$lat0,
    map  = .ncr_state$map
  )
}

#' @rdname nc_set
#' @export
nc_reset <- function() {
  .init_ncr_state()
  invisible(nc_get())
}

missing_arg <- function() {
  structure(list(), class = "ncr_missing")
}

is_missing_arg <- function(x) {
  inherits(x, "ncr_missing")
}
