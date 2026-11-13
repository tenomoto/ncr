#' Create a Synthetic 4D Atmospheric NetCDF Dataset
#'
#' Generates a CF-compliant 4D (`lon`, `lat`, `lev`, `time`) NetCDF file using
#' `RNetCDF` with realistic synthetic meteorological fields: air temperature
#' (`temp`), zonal wind (`u`), meridional wind (`v`), and sea-level pressure
#' (`slp`). Useful for interactive exploration, examples, and tests.
#'
#' @param file Destination file path (defaults to a temporary `.nc` file).
#' @return The file path `file` (invisibly or visibly as a character string).
#' @examples
#' f <- nc_example()
#' nc <- nc_open(f)
#' print(nc)
#' nc_close(nc)
#' @export
nc_example <- function(file = tempfile(pattern = "ncr_demo_", fileext = ".nc")) {
  lon  <- seq(0, 350, by = 10)                       # 36 longitudes (0..350 E)
  lat  <- seq(-80, 80, by = 10)                      # 17 latitudes (-80..80 N)
  lev  <- c(1000, 850, 700, 500, 300, 200, 100)      # 7 pressure levels (hPa)
  time <- seq(0, by = 6, length.out = 12)            # 12 6-hourly steps (3 days)

  nlon <- length(lon)
  nlat <- length(lat)
  nlev <- length(lev)
  nt   <- length(time)

  nc <- RNetCDF::create.nc(file, format = "classic")
  on.exit(try(RNetCDF::close.nc(nc), silent = TRUE), add = TRUE)

  RNetCDF::dim.def.nc(nc, "lon",  nlon)
  RNetCDF::dim.def.nc(nc, "lat",  nlat)
  RNetCDF::dim.def.nc(nc, "lev",  nlev)
  RNetCDF::dim.def.nc(nc, "time", unlim = TRUE)

  # Coordinate variables
  RNetCDF::var.def.nc(nc, "lon",  "NC_DOUBLE", "lon")
  RNetCDF::att.put.nc(nc, "lon",  "units", "NC_CHAR", "degrees_east")
  RNetCDF::att.put.nc(nc, "lon",  "long_name", "NC_CHAR", "Longitude")
  RNetCDF::att.put.nc(nc, "lon",  "axis", "NC_CHAR", "X")
  RNetCDF::var.put.nc(nc, "lon",  lon)

  RNetCDF::var.def.nc(nc, "lat",  "NC_DOUBLE", "lat")
  RNetCDF::att.put.nc(nc, "lat",  "units", "NC_CHAR", "degrees_north")
  RNetCDF::att.put.nc(nc, "lat",  "long_name", "NC_CHAR", "Latitude")
  RNetCDF::att.put.nc(nc, "lat",  "axis", "NC_CHAR", "Y")
  RNetCDF::var.put.nc(nc, "lat",  lat)

  RNetCDF::var.def.nc(nc, "lev",  "NC_DOUBLE", "lev")
  RNetCDF::att.put.nc(nc, "lev",  "units", "NC_CHAR", "hPa")
  RNetCDF::att.put.nc(nc, "lev",  "long_name", "NC_CHAR", "Pressure Level")
  RNetCDF::att.put.nc(nc, "lev",  "axis", "NC_CHAR", "Z")
  RNetCDF::att.put.nc(nc, "lev",  "positive", "NC_CHAR", "down")
  RNetCDF::var.put.nc(nc, "lev",  lev)

  RNetCDF::var.def.nc(nc, "time", "NC_DOUBLE", "time")
  RNetCDF::att.put.nc(nc, "time", "units", "NC_CHAR", "hours since 2026-01-01 00:00:00")
  RNetCDF::att.put.nc(nc, "time", "long_name", "NC_CHAR", "Time")
  RNetCDF::att.put.nc(nc, "time", "axis", "NC_CHAR", "T")
  RNetCDF::att.put.nc(nc, "time", "calendar", "NC_CHAR", "standard")
  RNetCDF::var.put.nc(nc, "time", time)

  # Data variables
  dims4 <- c("lon", "lat", "lev", "time")
  dims3 <- c("lon", "lat", "time")

  RNetCDF::var.def.nc(nc, "temp", "NC_FLOAT", dims4)
  RNetCDF::att.put.nc(nc, "temp", "units", "NC_CHAR", "K")
  RNetCDF::att.put.nc(nc, "temp", "long_name", "NC_CHAR", "Air Temperature")

  RNetCDF::var.def.nc(nc, "u", "NC_FLOAT", dims4)
  RNetCDF::att.put.nc(nc, "u", "units", "NC_CHAR", "m/s")
  RNetCDF::att.put.nc(nc, "u", "long_name", "NC_CHAR", "Zonal Wind")

  RNetCDF::var.def.nc(nc, "v", "NC_FLOAT", dims4)
  RNetCDF::att.put.nc(nc, "v", "units", "NC_CHAR", "m/s")
  RNetCDF::att.put.nc(nc, "v", "long_name", "NC_CHAR", "Meridional Wind")

  RNetCDF::var.def.nc(nc, "slp", "NC_FLOAT", dims3)
  RNetCDF::att.put.nc(nc, "slp", "units", "NC_CHAR", "hPa")
  RNetCDF::att.put.nc(nc, "slp", "long_name", "NC_CHAR", "Sea Level Pressure")

  RNetCDF::att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR",
                      "Synthetic 4D Baroclinic Atmosphere Dataset for ncr")
  RNetCDF::att.put.nc(nc, "NC_GLOBAL", "Conventions", "NC_CHAR", "CF-1.8")

  # Populate synthetic 4D and 3D arrays
  temp <- array(0, dim = c(nlon, nlat, nlev, nt))
  u    <- array(0, dim = c(nlon, nlat, nlev, nt))
  v    <- array(0, dim = c(nlon, nlat, nlev, nt))
  slp  <- array(0, dim = c(nlon, nlat, nt))

  deg2rad <- pi / 180
  LON <- array(rep(lon, times = nlat), dim = c(nlon, nlat)) * deg2rad
  LAT <- array(rep(lat, each  = nlon), dim = c(nlon, nlat)) * deg2rad

  for (it in seq_len(nt)) {
    phase <- (time[it] / 72) * 2 * pi
    wave  <- sin(4 * LON - phase) * cos(LAT)^2
    slp[, , it] <- 1013.25 - 12 * sin(2 * LAT)^2 + 10 * wave

    for (k in seq_len(nlev)) {
      p_frac <- lev[k] / 1000
      # Temperature decreases with altitude and toward poles, plus eastward wave
      temp[, , k, it] <- 215 + (85 * cos(LAT) * p_frac^0.286) +
        6 * sin(3 * LON - phase - (1 - p_frac)) * cos(LAT)
      # Jet stream peaking near 250 hPa and +/- 35 deg latitude
      jet_v <- exp(-((lev[k] - 250) / 280)^2)
      u[, , k, it] <- (15 + 30 * jet_v) * sin(2 * LAT)^2 +
        8 * cos(4 * LON - phase) * cos(LAT)
      v[, , k, it] <- 12 * (1 + 0.8 * jet_v) * cos(4 * LON - phase) * cos(LAT)^2
    }
  }

  RNetCDF::var.put.nc(nc, "temp", temp)
  RNetCDF::var.put.nc(nc, "u",    u)
  RNetCDF::var.put.nc(nc, "v",    v)
  RNetCDF::var.put.nc(nc, "slp",  slp)
  RNetCDF::sync.nc(nc)

  file
}
