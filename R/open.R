#' Open, Close, and Inspect a NetCDF Dataset
#'
#' `nc_open()` opens a NetCDF dataset using `RNetCDF` and returns an S3 object
#' of class `"ncr"` containing parsed dimension coordinates, CF metadata, and
#' decoded timestamps. `nc_close()` closes the underlying NetCDF handle, and
#' `nc_meta()` returns a structured metadata summary of class `"ncr_meta"`.
#'
#' @param path Path to a NetCDF file (for `nc_open()`).
#' @param write Logical; if `TRUE`, opens the file for writing (default `FALSE`).
#' @param x An `"ncr"` object returned by `nc_open()`.
#' @param object An `"ncr"` object (for `summary()`).
#' @param ... Additional arguments passed to print/summary methods.
#' @return
#' * `nc_open()` returns an S3 object of class `"ncr"`.
#' * `nc_close()` invisibly returns `NULL`.
#' * `nc_meta()` returns an S3 object of class `"ncr_meta"`.
#' @examples
#' f <- nc_example()
#' nc <- nc_open(f)
#' print(nc)
#' nc_meta(nc)
#' nc_close(nc)
#' @export
nc_open <- function(path, write = FALSE) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("'path' must be a non-empty character string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("NetCDF file not found: '%s'", path), call. = FALSE)
  }

  handle <- RNetCDF::open.nc(path, write = isTRUE(write))
  env <- new.env(parent = emptyenv())
  env$handle <- handle
  env$closed <- FALSE
  env$path   <- normalizePath(path, mustWork = FALSE)

  reg.finalizer(env, function(e) {
    if (!isTRUE(e$closed) && !is.null(e$handle)) {
      try(RNetCDF::close.nc(e$handle), silent = TRUE)
      e$closed <- TRUE
    }
  }, onexit = TRUE)

  meta <- .parse_nc_meta(handle, env$path)

  obj <- structure(
    list(
      env  = env,
      path = env$path,
      meta = meta
    ),
    class = "ncr"
  )
  obj
}

#' @rdname nc_open
#' @export
nc_close <- function(x) {
  if (!inherits(x, "ncr")) {
    stop("'x' must be an 'ncr' object.", call. = FALSE)
  }
  if (!isTRUE(x$env$closed) && !is.null(x$env$handle)) {
    RNetCDF::close.nc(x$env$handle)
    x$env$closed <- TRUE
    x$env$handle <- NULL
  }
  invisible(NULL)
}

#' @rdname nc_open
#' @export
nc_meta <- function(x) {
  if (!inherits(x, "ncr")) {
    stop("'x' must be an 'ncr' object.", call. = FALSE)
  }
  x$meta
}

#' @rdname nc_open
#' @export
names.ncr <- function(x) {
  names(x$meta$vars)
}

#' @rdname nc_open
#' @export
print.ncr <- function(x, ...) {
  status <- if (isTRUE(x$env$closed)) "closed" else "open"
  cat(sprintf("<ncr NetCDF [%s]> %s\n", status, basename(x$path)))
  dims <- x$meta$dims
  if (length(dims) > 0L) {
    dim_str <- vapply(dims, function(d) {
      sprintf("%s(%d)[%s]", d$name, d$length, d$axis)
    }, character(1))
    cat("  Dimensions:", paste(dim_str, collapse = ", "), "\n")
  }
  vars <- x$meta$vars
  if (length(vars) > 0L) {
    cat("  Variables:\n")
    for (v in vars) {
      dnames <- paste(v$dimnames, collapse = ", ")
      unit_str <- if (nzchar(v$units)) sprintf(" [%s]", v$units) else ""
      lname <- if (nzchar(v$long_name)) sprintf(" - %s", v$long_name) else ""
      cat(sprintf("    %-10s (%s)%s%s\n", v$name, dnames, unit_str, lname))
    }
  }
  st <- nc_get()
  map_str <- if (is.list(st$map)) "<custom>" else paste(as.character(st$map), collapse = "+")
  col_str <- if (is.function(st$col)) "<fn>" else if (is.null(st$col)) "default" else paste(as.character(st$col), collapse = "+")
  cat(sprintf("  State:      lon=%s, lat=%s, lev=%s, time=%s, proj=%s, map=%s, col=%s\n",
              .fmt_state_val(st$lon),
              .fmt_state_val(st$lat),
              .fmt_state_val(st$lev),
              .fmt_state_val(st$time),
              as.character(st$proj),
              map_str,
              col_str))
  invisible(x)
}

#' @rdname nc_open
#' @export
summary.ncr <- function(object, ...) {
  print.ncr_meta(object$meta, ...)
  invisible(object$meta)
}

#' @rdname nc_open
#' @export
print.ncr_meta <- function(x, ...) {
  cat("NetCDF File:", x$path, "\n")
  if (length(x$gatts) > 0L) {
    cat("Global Attributes:\n")
    for (nm in names(x$gatts)) {
      val <- paste(as.character(x$gatts[[nm]]), collapse = ", ")
      if (nchar(val) > 65L) val <- paste0(substr(val, 1L, 62L), "...")
      cat(sprintf("  %s: %s\n", nm, val))
    }
  }
  cat("\nDimensions:\n")
  for (d in x$dims) {
    rng <- if (d$axis == "time" && !is.null(d$time_posix)) {
      paste(format(range(d$time_posix, na.rm = TRUE), "%Y-%m-%d %H:%M"), collapse = " .. ")
    } else if (length(d$vals) > 0L && is.numeric(d$vals)) {
      paste(signif(range(d$vals, na.rm = TRUE), 5), collapse = " .. ")
    } else {
      ""
    }
    u <- if (nzchar(d$units)) sprintf(" (%s)", d$units) else ""
    cat(sprintf("  %-10s len=%-5d axis=%-5s range=[%s]%s\n",
                d$name, d$length, d$axis, rng, u))
  }
  cat("\nData Variables:\n")
  for (v in x$vars) {
    dnames <- paste(v$dimnames, collapse = ", ")
    u <- if (nzchar(v$units)) sprintf(" [%s]", v$units) else ""
    ln <- if (nzchar(v$long_name)) sprintf(" : %s", v$long_name) else ""
    cat(sprintf("  %-12s (%s) <%s>%s%s\n", v$name, dnames, v$type, u, ln))
  }
  invisible(x)
}

.fmt_state_val <- function(v) {
  if (is.null(v)) return("all")
  if (inherits(v, "AsIs")) return( paste0("I(", paste(v, collapse = ":"), ")") )
  paste(as.character(v), collapse = ":")
}

.safe_att_get <- function(nc, var, att) {
  tryCatch(
    RNetCDF::att.get.nc(nc, var, att),
    error = function(e) NULL
  )
}

.detect_axis <- function(name, units, axis_att, std_name, positive) {
  nm <- tolower(name)
  un <- tolower(if (is.null(units)) "" else as.character(units))
  ax <- toupper(if (is.null(axis_att)) "" else as.character(axis_att))
  sn <- tolower(if (is.null(std_name)) "" else as.character(std_name))

  if (ax == "X" || sn %in% c("longitude", "grid_longitude") ||
      grepl("degree.*(east|_e|e$)", un) ||
      nm %in% c("lon", "longitude", "lons", "x", "i", "rlon", "nav_lon")) {
    return("lon")
  }
  if (ax == "Y" || sn %in% c("latitude", "grid_latitude") ||
      grepl("degree.*(north|_n|n$)", un) ||
      nm %in% c("lat", "latitude", "lats", "y", "j", "rlat", "nav_lat")) {
    return("lat")
  }
  if (ax == "T" || sn == "time" ||
      grepl("(seconds|minutes|hours|days|months|years)\\s+since", un) ||
      nm %in% c("time", "t", "date", "datetime", "valid_time")) {
    return("time")
  }
  if (ax == "Z" || !is.null(positive) ||
      sn %in% c("air_pressure", "altitude", "height", "depth", "atmosphere_ln_pressure_coordinate") ||
      grepl("^(hpa|pa|mb|millibar|mbar|dbar|m|km|level|layer|sigma|eta)$", un) ||
      nm %in% c("lev", "level", "levels", "plev", "pressure", "isobaric",
                "z", "k", "height", "depth", "alt", "altitude", "bottom_top")) {
    return("lev")
  }
  "other"
}

.decode_time <- function(vals, units) {
  if (is.null(units) || !nzchar(units) || !is.numeric(vals)) return(NULL)
  if (!grepl("since", units, ignore.case = TRUE)) return(NULL)
  tryCatch({
    res <- RNetCDF::utcal.nc(units, vals, type = "c")
    if (inherits(res, "POSIXt")) {
      return(as.POSIXct(res, tz = "UTC"))
    }
    cal <- RNetCDF::utcal.nc(units, vals, type = "n")
    if (is.matrix(cal) && ncol(cal) >= 6L) {
      iso <- sprintf("%04d-%02d-%02d %02d:%02d:%02d",
                     cal[, 1], cal[, 2], cal[, 3],
                     cal[, 4], cal[, 5], floor(cal[, 6]))
      as.POSIXct(iso, tz = "UTC")
    } else {
      NULL
    }
  }, error = function(e) NULL)
}

.parse_nc_meta <- function(nc, path) {
  finq <- RNetCDF::file.inq.nc(nc)

  # Global attributes
  gatts <- list()
  if (finq$ngatts > 0L) {
    for (i in seq.int(0L, finq$ngatts - 1L)) {
      ainq <- tryCatch(RNetCDF::att.inq.nc(nc, "NC_GLOBAL", i), error = function(e) NULL)
      if (!is.null(ainq)) {
        gatts[[ainq$name]] <- .safe_att_get(nc, "NC_GLOBAL", ainq$name)
      }
    }
  }

  # All variable names first so we know which dimensions have coordinate variables
  all_var_info <- vector("list", finq$nvars)
  var_names <- character(finq$nvars)
  for (v in seq_len(finq$nvars)) {
    vinq <- RNetCDF::var.inq.nc(nc, v - 1L)
    all_var_info[[v]] <- vinq
    var_names[v] <- vinq$name
  }

  # Dimensions
  dims <- list()
  dim_names_by_id <- character(finq$ndims)
  for (d in seq_len(finq$ndims)) {
    dinq <- RNetCDF::dim.inq.nc(nc, d - 1L)
    dname <- dinq$name
    dim_names_by_id[d] <- dname

    has_cvar <- dname %in% var_names
    vals <- if (has_cvar && dinq$length > 0L) {
      tryCatch(as.vector(RNetCDF::var.get.nc(nc, dname)), error = function(e) seq_len(dinq$length))
    } else {
      seq_len(dinq$length)
    }
    units    <- if (has_cvar) .safe_att_get(nc, dname, "units") else NULL
    axis_att <- if (has_cvar) .safe_att_get(nc, dname, "axis") else NULL
    std_name <- if (has_cvar) .safe_att_get(nc, dname, "standard_name") else NULL
    positive <- if (has_cvar) .safe_att_get(nc, dname, "positive") else NULL
    long_nm  <- if (has_cvar) .safe_att_get(nc, dname, "long_name") else NULL

    axis_role <- .detect_axis(dname, units, axis_att, std_name, positive)
    time_posix <- if (axis_role == "time") .decode_time(vals, units) else NULL

    dims[[dname]] <- list(
      id         = dinq$id,
      name       = dname,
      length     = dinq$length,
      unlim      = isTRUE(dinq$unlim),
      vals       = vals,
      units      = if (is.null(units)) "" else as.character(units)[1L],
      long_name  = if (is.null(long_nm)) dname else as.character(long_nm)[1L],
      positive   = if (is.null(positive)) "" else as.character(positive)[1L],
      axis       = axis_role,
      time_posix = time_posix
    )
  }

  # Data variables (variables that are not pure 1D coordinate variables, though
  # if a file only has 1D variables we keep non-coordinate ones too)
  vars <- list()
  coord_names <- names(dims)
  for (v in seq_len(finq$nvars)) {
    vinq <- all_var_info[[v]]
    vname <- vinq$name
    if (vname %in% coord_names && vinq$ndims <= 1L) next
    dids <- vinq$dimids
    dnames <- if (length(dids) > 0L && !all(is.na(dids))) {
      dim_names_by_id[dids + 1L]
    } else {
      character(0)
    }
    axes <- vapply(dnames, function(dn) dims[[dn]]$axis, character(1))
    units   <- .safe_att_get(nc, vname, "units")
    long_nm <- .safe_att_get(nc, vname, "long_name")
    std_nm  <- .safe_att_get(nc, vname, "standard_name")

    vars[[vname]] <- list(
      id            = vinq$id,
      name          = vname,
      type          = vinq$type,
      ndims         = vinq$ndims,
      dimids        = dids,
      dimnames      = dnames,
      axes          = unname(axes),
      units         = if (is.null(units)) "" else as.character(units)[1L],
      long_name     = if (is.null(long_nm)) vname else as.character(long_nm)[1L],
      standard_name = if (is.null(std_nm)) "" else as.character(std_nm)[1L]
    )
  }

  # Map canonical axis names ("lon", "lat", "lev", "time") to dimension names
  axis_map <- list()
  for (d in dims) {
    if (d$axis %in% c("lon", "lat", "lev", "time") && is.null(axis_map[[d$axis]])) {
      axis_map[[d$axis]] <- d$name
    }
  }

  structure(
    list(
      path     = path,
      dims     = dims,
      vars     = vars,
      gatts    = gatts,
      axis_map = axis_map
    ),
    class = "ncr_meta"
  )
}
