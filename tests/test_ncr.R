library(ncr)

# 1. Generate synthetic 4D dataset and open
f <- nc_example()
stopifnot(file.exists(f))

nc <- nc_open(f)
stopifnot(inherits(nc, "ncr"))
stopifnot(all(c("temp", "u", "v", "slp") %in% names(nc)))

meta <- nc_meta(nc)
stopifnot(inherits(meta, "ncr_meta"))
print(nc)
summary(nc)

# 2. Global state defaults and updates
nc_reset()
st <- nc_get()
stopifnot(identical(st$time, 1L))
stopifnot(identical(st$proj, "lonlat"))

# 3. Slicing syntax: $, [, nc_slice
s1 <- nc$temp
s2 <- nc["temp"]
stopifnot(inherits(s1, "ncr_slice"))
stopifnot(identical(dim(s1$z), c(36L, 17L)))
stopifnot(identical(s1$z, s2$z))
print(s1)

# Slice at 500 hPa and time step 3
s500 <- nc["temp", lev = 500, time = 3]
stopifnot(s500$fixed_dims[[1]]$vals == 500)

# Slice by timestamp string
s_time <- nc["temp", time = "2026-01-02 00:00:00"]
stopifnot(inherits(s_time, "ncr_slice"))

# 3b. Arithmetic, Math, and Summary operations on ncr_slice
t_c <- nc$temp - 273.15
stopifnot(inherits(t_c, "ncr_slice"))
stopifnot(identical(dim(t_c$z), dim(s1$z)))
stopifnot(max(abs(t_c$z - (s1$z - 273.15))) < 1e-10)
stopifnot(identical(range(t_c), range(t_c$z)))

dt <- nc["temp", lev = 500] - nc["temp", lev = 850]
stopifnot(inherits(dt, "ncr_slice") && identical(dt$units, "K"))

spd <- sqrt(nc$u^2 + nc$v^2)
stopifnot(inherits(spd, "ncr_slice") && all(spd$z >= 0))

# Stateful slicing via nc_set
nc_set(lon = c(100, 200), lat = c(0, 60), lev = 850, time = 2)
s_sub <- nc$temp
stopifnot(all(s_sub$x >= 100 & s_sub$x <= 200))
stopifnot(all(s_sub$y >= 0 & s_sub$y <= 60))
nc_reset()

# Vector slice and scalar + vector slice
vec <- nc[c("u", "v"), lev = 850]
stopifnot(isTRUE(vec$is_vector) && isTRUE(vec$has_uv))
sv <- nc[c("temp", "u", "v"), lev = 850]
stopifnot(!isTRUE(sv$is_vector) && isTRUE(sv$has_uv))

# 4. Projections and pressure/altitude conversion
z_km <- pres2alt(c(1000, 500, 100))
p_back <- alt2pres(z_km)
stopifnot(max(abs(p_back - c(1000, 500, 100))) < 1e-8)

for (p in c("lonlat", "mercator", "robinson", "mollweide",
            "orthographic", "stereographic", "npstere", "spstere", "lambert")) {
  pr <- nc_project(c(120, 140, 160), c(20, 35, 50), proj = p, lon0 = 140, lat0 = 35)
  stopifnot(length(pr$x) == 3L && length(pr$y) == 3L)
}

# 5. Graphics rendering tests on offscreen PDF device
pdf_file <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf_file, width = 8, height = 5)

# 2D filled & line contour maps
plot(nc, "temp")
plot(nc, "temp", map = FALSE)
plot(nc, "temp", map = list(x = c(130, 145, 145, 130), y = c(30, 30, 45, 45)))
contour(nc, "temp", lev = 500, filled = FALSE)
contour(nc$temp, lev = 500, filled = FALSE)
contour(nc$temp, filled = TRUE)

# Vectors and scalar + vector overlay
plot(nc, c("u", "v"), lev = 850, skip = 2)
plot(nc, "temp", uv = c("u", "v"), lev = 850, skip = 2)
plot(nc$slp, contour = TRUE)
points(nc[c("u", "v"), lev = 1000], skip = 2)
drawvec(nc["u", lev = 850], nc["v", lev = 850], skip = 2, col = "navy")
drawvec(nc[c("u", "v"), lev = 500], skip = 2)

# Map projections
plot(nc, "temp", proj = "orthographic", lon0 = 140, lat0 = 30, uv = c("u", "v"), skip = 2)
plot(nc, "temp", proj = "robinson")
plot(nc, "temp", proj = "npstere")
plot(nc, "temp", proj = "npstere", lat = c(30, 90))
plot(nc, "temp", proj = "spstere", lat = c(-90, -30))
plot(nc, "temp", proj = "lambert", lon = c(100, 180), lat = c(15, 65))

# Hovmoeller diagrams (time on Y and time on X)
nc_hovmoeller(nc, "temp", lat = 0, time_axis = "y")
nc_hovmoeller(nc, "temp", lat = 0, time_axis = "x")

# Vertical cross-sections (log-pressure and linear pressure, with dual pressure/height axes)
nc_xsection(nc, "temp", lon = 140, zlog = TRUE)
nc_xsection(nc, "u", lon = 140, zlog = FALSE)

# 1D vertical profile and time series
plot(nc["temp", lon = 140, lat = 35, lev = NULL])
plot(nc["temp", lon = 140, lat = 35, time = NULL])

grDevices::dev.off()
stopifnot(file.exists(pdf_file) && file.info(pdf_file)$size > 0)

# 6. Close handle cleanly
nc_close(nc)
stopifnot(isTRUE(nc$env$closed))
unlink(c(f, pdf_file))
cat("All ncr tests passed!\n")
