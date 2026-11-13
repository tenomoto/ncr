#' Slice a NetCDF Dataset Using Global State or Explicit Bounds
#'
#' Extracts a 1D or 2D slice (scalar, vector pair `(u, v)`, or scalar with
#' vector overlay) from an `"ncr"` dataset. Unspecified dimensions automatically
#' inherit the global state managed by [nc_set()] / [nc_get()].
#'
#' @param x An `"ncr"` object returned by [nc_open()].
#' @param var Character vector of variable name(s):
#'   * Length 1 (`"temp"`): extracts a scalar slice.
#'   * Length 2 (`c("u", "v")`): extracts a vector field `(u, v)` (plus speed
#'     magnitude `z = sqrt(u^2 + v^2)`).
#'   * Length 3 (`c("temp", "u", "v")`): extracts scalar `"temp"` with vector
#'     overlay `("u", "v")`.
#' @param name Variable name for `$` extraction (e.g., `nc$temp`).
#' @param lon,lat,lev,time Optional dimension selectors overriding the current
#'   global state for this slice:
#'   * `NULL` or `NA`: full extent of the dimension.
#'   * Single value (e.g., `lev = 500`, `lat = 35`, `time = "2020-01-02"`):
#'     nearest coordinate value (or 1-based time index when `time` is an integer
#'     in `1..N`).
#'   * Length-2 vector `c(min, max)`: closed coordinate interval `[min, max]`.
#'   * `I(indices)`: explicit 1-based integer index or index vector (e.g.,
#'     `I(1:5)`).
#' @param uv Optional character vector of length 2 (`c("u", "v")`) to attach a
#'   vector field overlay to a scalar slice.
#' @param drop Logical; if `TRUE` (default), drops singleton dimensions.
#' @param e1,e2 Operands for `Ops.ncr_slice()` (an `"ncr_slice"` object or a
#'   numeric scalar/array).
#' @param na.rm Logical; whether to remove `NA` values in `Summary.ncr_slice()`.
#' @param ... Additional dimension selectors matched by exact dimension name in
#'   the NetCDF file.
#' @return An S3 object of class `"ncr_slice"`.
#' @examples
#' f <- nc_example()
#' nc <- nc_open(f)
#'
#' # Default 2D lon-lat slice at lev=1, time=1
#' s1 <- nc["temp"]
#' s2 <- nc$temp
#'
#' # Slice at 500 hPa
#' s500 <- nc["temp", lev = 500]
#'
#' # Vector slice (u, v)
#' vec <- nc[c("u", "v"), lev = 850]
#'
#' # Vertical cross-section along 140E (all levels)
#' xs <- nc["temp", lon = 140, lev = NULL]
#'
#' # Hovmoeller slice along the equator (all times)
#' hov <- nc["temp", lat = 0, time = NULL]
#'
#' nc_close(nc)
#' @export
nc_slice <- function(x,
                     var = NULL,
                     lon = missing_arg(),
                     lat = missing_arg(),
                     lev = missing_arg(),
                     time = missing_arg(),
                     uv = NULL,
                     drop = TRUE,
                     ...) {
  if (!inherits(x, "ncr")) {
    stop("'x' must be an 'ncr' object.", call. = FALSE)
  }
  if (isTRUE(x$env$closed) || is.null(x$env$handle)) {
    stop("The NetCDF dataset is closed.", call. = FALSE)
  }

  avail_vars <- names(x$meta$vars)
  if (length(avail_vars) == 0L) {
    stop("No data variables found in NetCDF dataset.", call. = FALSE)
  }
  if (is.null(var)) {
    var <- avail_vars[1L]
  }
  if (!is.character(var) || length(var) < 1L || length(var) > 3L) {
    stop("'var' must be a character vector of length 1 (scalar), 2 (u, v), or 3 (scalar, u, v).",
         call. = FALSE)
  }
  if (length(var) == 3L) {
    uv  <- var[2:3]
    var <- var[1L]
  }

  if (length(var) == 2L) {
    su <- .slice_single_var(x, var[1L], lon = lon, lat = lat, lev = lev,
                            time = time, drop = drop, ...)
    sv <- .slice_single_var(x, var[2L], lon = lon, lat = lat, lev = lev,
                            time = time, drop = drop, ...)
    spd <- sqrt(su$z^2 + sv$z^2)
    su$u <- su$z
    su$v <- sv$z
    su$z <- spd
    su$var <- paste(var, collapse = ",")
    su$u_var <- var[1L]
    su$v_var <- var[2L]
    su$long_name <- sprintf("Vector (%s, %s)", var[1L], var[2L])
    su$is_vector <- TRUE
    su$has_uv    <- TRUE
    return(su)
  }

  s <- .slice_single_var(x, var[1L], lon = lon, lat = lat, lev = lev,
                         time = time, drop = drop, ...)
  if (!is.null(uv)) {
    if (!is.character(uv) || length(uv) != 2L) {
      stop("'uv' must be a character vector of length 2 naming c(u, v).", call. = FALSE)
    }
    su <- .slice_single_var(x, uv[1L], lon = lon, lat = lat, lev = lev,
                            time = time, drop = drop, ...)
    sv <- .slice_single_var(x, uv[2L], lon = lon, lat = lat, lev = lev,
                            time = time, drop = drop, ...)
    s$u <- su$z
    s$v <- sv$z
    s$u_var <- uv[1L]
    s$v_var <- uv[2L]
    s$has_uv <- TRUE
  }
  s
}

#' @rdname nc_slice
#' @export
`[.ncr` <- function(x,
                    var = NULL,
                    lon = missing_arg(),
                    lat = missing_arg(),
                    lev = missing_arg(),
                    time = missing_arg(),
                    uv = NULL,
                    drop = TRUE,
                    ...) {
  nc_slice(x, var = var, lon = lon, lat = lat, lev = lev, time = time,
           uv = uv, drop = drop, ...)
}

#' @rdname nc_slice
#' @export
`$.ncr` <- function(x, name) {
  if (name %in% c("env", "path", "meta")) {
    return(.subset2(x, name))
  }
  nc_slice(x, var = name)
}

#' @rdname nc_slice
#' @export
print.ncr_slice <- function(x, ...) {
  kind <- if (isTRUE(x$is_vector)) {
    sprintf("vector (%s, %s)", x$u_var, x$v_var)
  } else if (isTRUE(x$has_uv)) {
    sprintf("scalar '%s' + vector (%s, %s)", x$var, x$u_var, x$v_var)
  } else {
    sprintf("scalar '%s'", x$var)
  }
  u_str <- if (nzchar(x$units)) sprintf(" [%s]", x$units) else ""
  cat(sprintf("<ncr_slice: %s%s (%s)>\n", kind, u_str, x$plot_type))
  if (length(x$active_dims) > 0L) {
    for (d in x$active_dims) {
      rng <- if (d$axis == "time" && !is.null(d$time_posix)) {
        paste(format(range(d$time_posix, na.rm = TRUE), "%Y-%m-%d %H:%M"), collapse = " .. ")
      } else {
        paste(signif(range(d$vals, na.rm = TRUE), 5), collapse = " .. ")
      }
      cat(sprintf("  Axis %-6s (%s): len=%d [%s] %s\n",
                  d$name, d$axis, length(d$vals), rng, d$units))
    }
  }
  if (length(x$fixed_dims) > 0L) {
    fx <- vapply(x$fixed_dims, function(d) {
      val_str <- if (d$axis == "time" && !is.null(d$time_posix)) {
        format(d$time_posix[1L], "%Y-%m-%d %H:%M")
      } else {
        paste0(signif(d$vals[1L], 5), if (nzchar(d$units)) paste0(" ", d$units) else "")
      }
      sprintf("%s=%s", d$name, val_str)
    }, character(1))
    cat("  Fixed:     ", paste(fx, collapse = ", "), "\n")
  }
  if (!is.null(x$z) && any(!is.na(x$z))) {
    rng <- signif(range(x$z, na.rm = TRUE), 5)
    cat(sprintf("  Data range: [%s, %s]\n", rng[1L], rng[2L]))
  }
  invisible(x)
}

#' @rdname nc_slice
#' @export
Ops.ncr_slice <- function(e1, e2) {
  op <- .Generic
  if (missing(e2)) {
    res <- e1
    res$z <- do.call(op, list(e1$z))
    res$var <- sprintf("%s(%s)", op, e1$var)
    res$long_name <- sprintf("%s(%s)", op, e1$long_name)
    res$is_vector <- FALSE
    res$has_uv <- FALSE
    res$u <- NULL
    res$v <- NULL
    if (!(op %in% c("+", "-"))) res$units <- ""
    return(res)
  }

  s1 <- inherits(e1, "ncr_slice")
  s2 <- inherits(e2, "ncr_slice")
  if (s1 && s2) {
    if (!identical(dim(e1$z), dim(e2$z))) {
      stop("Incompatible slice dimensions for '", op, "'.", call. = FALSE)
    }
    res <- e1
    res$z <- do.call(op, list(e1$z, e2$z))
    res$var <- sprintf("(%s %s %s)", e1$var, op, e2$var)
    res$long_name <- sprintf("(%s %s %s)", e1$long_name, op, e2$long_name)
    res$units <- if ((op %in% c("+", "-")) && identical(e1$units, e2$units)) e1$units else ""
  } else if (s1) {
    res <- e1
    res$z <- do.call(op, list(e1$z, e2))
    lbl <- if (is.numeric(e2) && length(e2) == 1L) format(e2) else "expr"
    res$var <- sprintf("(%s %s %s)", e1$var, op, lbl)
    res$long_name <- sprintf("(%s %s %s)", e1$long_name, op, lbl)
    res$units <- ""
  } else {
    res <- e2
    res$z <- do.call(op, list(e1, e2$z))
    lbl <- if (is.numeric(e1) && length(e1) == 1L) format(e1) else "expr"
    res$var <- sprintf("(%s %s %s)", lbl, op, e2$var)
    res$long_name <- sprintf("(%s %s %s)", lbl, op, e2$long_name)
    res$units <- ""
  }

  res$is_vector <- FALSE
  res$has_uv <- FALSE
  res$u <- NULL
  res$v <- NULL
  res
}

#' @rdname nc_slice
#' @export
Math.ncr_slice <- function(x, ...) {
  res <- x
  res$z <- get(.Generic, mode = "function")(x$z, ...)
  res$var <- sprintf("%s(%s)", .Generic, x$var)
  res$long_name <- sprintf("%s(%s)", .Generic, x$long_name)
  res$units <- ""
  res$is_vector <- FALSE
  res$has_uv <- FALSE
  res$u <- NULL
  res$v <- NULL
  res
}

#' @rdname nc_slice
#' @export
Summary.ncr_slice <- function(..., na.rm = FALSE) {
  args <- list(...)
  vals <- lapply(args, function(a) if (inherits(a, "ncr_slice")) a$z else a)
  do.call(get(.Generic, mode = "function"), c(vals, list(na.rm = na.rm)))
}

.slice_single_var <- function(x, var, lon, lat, lev, time, drop = TRUE, ...) {
  vmeta <- x$meta$vars[[var]]
  if (is.null(vmeta)) {
    stop(sprintf("Variable '%s' not found. Available variables: %s",
                 var, paste(names(x$meta$vars), collapse = ", ")),
         call. = FALSE)
  }

  st <- nc_get()
  req <- list(
    lon  = if (is_missing_arg(lon))  st$lon  else lon,
    lat  = if (is_missing_arg(lat))  st$lat  else lat,
    lev  = if (is_missing_arg(lev))  st$lev  else lev,
    time = if (is_missing_arg(time)) st$time else time
  )
  dots <- list(...)

  ndims <- vmeta$ndims
  if (ndims == 0L) {
    val <- RNetCDF::var.get.nc(x$env$handle, var, unpack = TRUE)
    return(structure(
      list(var = var, z = val, units = vmeta$units, long_name = vmeta$long_name,
           active_dims = list(), fixed_dims = list(), plot_type = "scalar",
           is_vector = FALSE, has_uv = FALSE),
      class = "ncr_slice"
    ))
  }

  idx_list <- vector("list", ndims)
  dim_slices <- vector("list", ndims)

  for (i in seq_len(ndims)) {
    dname <- vmeta$dimnames[i]
    dmeta <- x$meta$dims[[dname]]
    axis  <- dmeta$axis

    # Priority: exact dim name in `...`, then canonical axis in (lon, lat, lev, time)
    spec <- if (dname %in% names(dots)) {
      dots[[dname]]
    } else if (axis %in% names(req)) {
      req[[axis]]
    } else {
      I(1L)
    }

    idx <- .resolve_dim_indices(dmeta, spec)
    idx_list[[i]] <- idx

    sub_vals <- dmeta$vals[idx]
    sub_time <- if (!is.null(dmeta$time_posix)) dmeta$time_posix[idx] else NULL
    dim_slices[[i]] <- list(
      name       = dname,
      axis       = axis,
      units      = dmeta$units,
      long_name  = dmeta$long_name,
      positive   = dmeta$positive,
      vals       = sub_vals,
      time_posix = sub_time,
      indices    = idx
    )
  }

  # Read minimal bounding slab from RNetCDF and subset if non-contiguous
  start <- vapply(idx_list, min, integer(1))
  count <- vapply(idx_list, function(ix) max(ix) - min(ix) + 1L, integer(1))

  raw <- RNetCDF::var.get.nc(x$env$handle, var, start = start, count = count,
                             collapse = FALSE, unpack = TRUE)

  # Subset exact relative indices if any dimension is non-contiguous or reordered
  rel_idx <- lapply(seq_len(ndims), function(i) idx_list[[i]] - start[i] + 1L)
  arr <- do.call(`[`, c(list(raw), rel_idx, list(drop = FALSE)))

  # Identify active (length > 1) vs fixed (length == 1) dimensions
  lens <- vapply(idx_list, length, integer(1))
  active_pos <- which(lens > 1L)
  fixed_pos  <- which(lens == 1L)

  # If all dimensions are length 1, keep the first one as active
  if (length(active_pos) == 0L) {
    active_pos <- 1L
    fixed_pos  <- setdiff(seq_len(ndims), active_pos)
  }
  if (length(active_pos) > 2L && isTRUE(drop)) {
    stop(sprintf(
      "Slice of '%s' has %d varying dimensions (%s). Fix %d more dimension(s) via nc_set() or slice arguments.",
      var, length(active_pos),
      paste(vmeta$dimnames[active_pos], collapse = ", "),
      length(active_pos) - 2L
    ), call. = FALSE)
  }

  # Collapse singleton dimensions into 1D or 2D array
  if (isTRUE(drop)) {
    arr <- array(arr, dim = lens[active_pos])
  }

  active_dims <- dim_slices[active_pos]
  fixed_dims  <- dim_slices[fixed_pos]

  # Ensure canonical orientation for 2D slices:
  # - lon-lat: x = lon, y = lat
  # - xsection: x = lon or lat, y = lev
  # - hovmoeller: x = lon or lat, y = time (by default)
  plot_type <- "1d"
  if (length(active_dims) == 2L) {
    ax <- c(active_dims[[1L]]$axis, active_dims[[2L]]$axis)
    swap <- FALSE
    if (setequal(ax, c("lon", "lat"))) {
      plot_type <- "lonlat"
      if (ax[1L] == "lat" && ax[2L] == "lon") swap <- TRUE
    } else if ("lev" %in% ax && any(c("lon", "lat") %in% ax)) {
      plot_type <- "xsection"
      if (ax[1L] == "lev") swap <- TRUE
    } else if ("time" %in% ax && any(c("lon", "lat") %in% ax)) {
      plot_type <- "hovmoeller"
      if (ax[1L] == "time") swap <- TRUE
    } else {
      plot_type <- "2d"
    }

    if (swap) {
      arr <- t(arr)
      active_dims <- active_dims[c(2L, 1L)]
    }

    # Ensure both x and y coordinates are strictly increasing for contour / .filled.contour
    x_vals <- active_dims[[1L]]$vals
    y_vals <- active_dims[[2L]]$vals
    if (length(x_vals) > 1L && x_vals[1L] > x_vals[length(x_vals)]) {
      ord_x <- rev(seq_along(x_vals))
      arr <- arr[ord_x, , drop = FALSE]
      active_dims[[1L]]$vals <- x_vals[ord_x]
      if (!is.null(active_dims[[1L]]$time_posix)) {
        active_dims[[1L]]$time_posix <- active_dims[[1L]]$time_posix[ord_x]
      }
      active_dims[[1L]]$reversed <- TRUE
    }
    if (length(y_vals) > 1L && y_vals[1L] > y_vals[length(y_vals)]) {
      ord_y <- rev(seq_along(y_vals))
      arr <- arr[, ord_y, drop = FALSE]
      active_dims[[2L]]$vals <- y_vals[ord_y]
      if (!is.null(active_dims[[2L]]$time_posix)) {
        active_dims[[2L]]$time_posix <- active_dims[[2L]]$time_posix[ord_y]
      }
      active_dims[[2L]]$reversed <- TRUE
    }
  } else if (length(active_dims) == 1L) {
    x_vals <- active_dims[[1L]]$vals
    if (length(x_vals) > 1L && x_vals[1L] > x_vals[length(x_vals)]) {
      ord_x <- rev(seq_along(x_vals))
      arr <- arr[ord_x]
      active_dims[[1L]]$vals <- x_vals[ord_x]
      if (!is.null(active_dims[[1L]]$time_posix)) {
        active_dims[[1L]]$time_posix <- active_dims[[1L]]$time_posix[ord_x]
      }
      active_dims[[1L]]$reversed <- TRUE
    }
  }

  x_dim <- active_dims[[1L]]
  y_dim <- if (length(active_dims) >= 2L) active_dims[[2L]] else NULL

  structure(
    list(
      var         = var,
      z           = arr,
      u           = NULL,
      v           = NULL,
      x           = x_dim$vals,
      y           = if (!is.null(y_dim)) y_dim$vals else NULL,
      x_dim       = x_dim,
      y_dim       = y_dim,
      active_dims = active_dims,
      fixed_dims  = fixed_dims,
      units       = vmeta$units,
      long_name   = vmeta$long_name,
      plot_type   = plot_type,
      is_vector   = FALSE,
      has_uv      = FALSE,
      meta        = x$meta
    ),
    class = "ncr_slice"
  )
}

.resolve_dim_indices <- function(dmeta, spec) {
  n <- dmeta$length
  vals <- dmeta$vals

  if (is.null(spec) || (length(spec) == 1L && is.na(spec))) {
    return(seq_len(n))
  }

  # Explicit 1-based indices via I(...)
  if (inherits(spec, "AsIs")) {
    idx <- as.integer(unclass(spec))
    idx <- idx[idx >= 1L & idx <= n]
    if (length(idx) == 0L) {
      stop(sprintf("Index out of bounds for dimension '%s' (length %d).", dmeta$name, n),
           call. = FALSE)
    }
    return(idx)
  }

  # Time axis resolution
  if (dmeta$axis == "time") {
    if (inherits(spec, c("POSIXt", "Date")) || is.character(spec)) {
      if (is.null(dmeta$time_posix)) {
        stop(sprintf("Dimension '%s' does not have decoded CF timestamps.", dmeta$name),
             call. = FALSE)
      }
      t_target <- as.POSIXct(spec, tz = "UTC")
      t_coords <- dmeta$time_posix
      if (length(t_target) == 1L) {
        return(which.min(abs(as.numeric(t_coords) - as.numeric(t_target))))
      } else {
        r <- range(as.numeric(t_target), na.rm = TRUE)
        idx <- which(as.numeric(t_coords) >= r[1L] & as.numeric(t_coords) <= r[2L])
        if (length(idx) == 0L) {
          idx <- seq.int(
            which.min(abs(as.numeric(t_coords) - r[1L])),
            which.min(abs(as.numeric(t_coords) - r[2L]))
          )
        }
        return(sort(unique(idx)))
      }
    }
    if (is.numeric(spec)) {
      # Integer indices 1..n are standard for time (e.g. default time = 1)
      if (all(spec >= 1 & spec <= n & spec == floor(spec))) {
        if (length(spec) == 2L && spec[2L] > spec[1L] + 1L) {
          return(seq.int(as.integer(spec[1L]), as.integer(spec[2L])))
        }
        return(as.integer(spec))
      }
    }
  }

  # Numeric coordinate matching (with longitude wrap support)
  if (is.numeric(spec)) {
    cvals <- vals
    if (dmeta$axis == "lon") {
      # Align spec with coordinate convention (0..360 vs -180..180)
      cmin <- min(cvals, na.rm = TRUE)
      cmax <- max(cvals, na.rm = TRUE)
      if (cmin >= 0 && any(spec < 0)) {
        spec <- ifelse(spec < 0, spec + 360, spec)
      } else if (cmax <= 180 && any(spec > 180)) {
        spec <- ifelse(spec > 180, spec - 360, spec)
      }
    }

    if (length(spec) == 1L) {
      # For lev, if spec == 1 and 1 is outside coordinate range, treat as index 1
      if (dmeta$axis == "lev" && spec == 1 && (spec < min(cvals, na.rm = TRUE) || spec > max(cvals, na.rm = TRUE))) {
        return(1L)
      }
      return(which.min(abs(cvals - spec)))
    } else if (length(spec) == 2L) {
      r <- range(spec, na.rm = TRUE)
      idx <- which(cvals >= r[1L] & cvals <= r[2L])
      if (length(idx) == 0L) {
        i1 <- which.min(abs(cvals - r[1L]))
        i2 <- which.min(abs(cvals - r[2L]))
        idx <- seq.int(min(i1, i2), max(i1, i2))
      }
      return(sort(unique(idx)))
    } else {
      idx <- unique(vapply(spec, function(v) which.min(abs(cvals - v)), integer(1)))
      return(sort(idx))
    }
  }

  stop(sprintf("Unsupported slice specification for dimension '%s'.", dmeta$name),
       call. = FALSE)
}
