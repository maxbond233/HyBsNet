# Internal utilities ---------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.hybs_abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

.assert_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    .hybs_abort("`", name, "` must be a data.frame.")
  }
  invisible(x)
}

.assert_columns <- function(data, required, source_name = deparse(substitute(data))) {
  .assert_data_frame(data, source_name)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    .hybs_abort(
      source_name, " is missing required column(s): ",
      paste(missing, collapse = ", "), "."
    )
  }
  invisible(data)
}

.assert_scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < lower || x > upper) {
    .hybs_abort(
      "`", name, "` must be one finite number in [", lower, ", ", upper, "]."
    )
  }
  invisible(x)
}

.assert_flag <- function(x, name) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    .hybs_abort("`", name, "` must be TRUE or FALSE.")
  }
  invisible(x)
}

.assert_nonempty_string <- function(x, name) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    .hybs_abort("`", name, "` must be one non-empty string.")
  }
  invisible(x)
}

.with_seed <- function(seed, code) {
  .assert_scalar_number(seed, "seed")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.atomic_write_csv <- function(data, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    .hybs_abort("Refusing to overwrite existing file: ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data, tmp, row.names = FALSE, na = "")
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(tmp, path)) {
    .hybs_abort("Could not move temporary file into place: ", path)
  }
  invisible(path)
}

.scaled_vector <- function(x) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(finite) == 0L) {
    return(out)
  }
  sx <- stats::sd(x[finite])
  if (!is.finite(sx) || sx == 0) {
    out[finite] <- 0
  } else {
    out[finite] <- (x[finite] - mean(x[finite])) / sx
  }
  out
}
