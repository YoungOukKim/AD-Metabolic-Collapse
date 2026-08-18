# ---------------------------------------------------------------------
# 00_utils.R - shared helpers.
#
# Two conventions are enforced throughout this repository:
#   (1) every script runs a self-test before it touches data and exits
#       non-zero if the self-test fails;
#   (2) verdict thresholds are declared before the data are read and are
#       never adjusted afterwards. A failed gate is reported as a failure.
# ---------------------------------------------------------------------

# Partial correlation of x and y given the columns of Z.
# Kept verbatim across scripts so that every reported partial r comes from
# a single implementation.
pc <- function(x, y, Z) {
  Z <- as.matrix(Z)
  o <- complete.cases(x, y, Z)
  x <- x[o]; y <- y[o]; Z <- Z[o, , drop = FALSE]
  r <- cor(residuals(lm(x ~ Z)), residuals(lm(y ~ Z)))
  n <- sum(o); k <- ncol(Z)
  c(r = r, p = 2 * pt(-abs(r * sqrt((n - k - 2) / (1 - r^2))), n - k - 2), n = n)
}

# Tolerance comparison with a relative epsilon.
# Without the epsilon, a deviation numerically equal to the tolerance can fall
# either side of the test: abs(-0.0747 - -0.0737) evaluates to
# 0.0010000000000000009 in IEEE double, which is greater than 0.0010.
gate <- function(obs, ref, tol) {
  d <- abs(obs - ref)
  list(pass   = all(is.finite(d)) && all(d <= tol * (1 + 1e-9)),
       maxdev = if (all(is.finite(d))) max(d) else NA_real_)
}

# Decode an AnnData categorical obs column.
# AnnData stores a missing category as code -1. In R, cats[codes + 1] silently
# drops index 0, which shortens the vector and misaligns it against every other
# obs column. The length assertion makes that failure loud instead of silent.
obs_cat <- function(h5, col, n_expect = NULL) {
  p1 <- paste0("obs/", col)
  dec <- function(codes, cats) {
    out <- rep(NA_character_, length(codes))
    ok <- !is.na(codes) & codes >= 0
    out[ok] <- as.character(cats)[codes[ok] + 1L]
    out
  }
  v <- if (h5$has(paste0(p1, "/codes")))
         dec(as.integer(h5$all(paste0(p1, "/codes"))), h5$all(paste0(p1, "/categories")))
       else if (h5$has(paste0("obs/__categories/", col)))
         dec(as.integer(h5$all(p1)), h5$all(paste0("obs/__categories/", col)))
       else
         as.character(h5$all(p1))
  if (!is.null(n_expect) && length(v) != n_expect)
    stop(sprintf("obs/%s has length %d but the file has %d cells - alignment is broken.",
                 col, length(v), n_expect))
  v
}

# Minimal HDF5 accessor. hdf5r is preferred; rhdf5 is a fallback.
# hdf5r's H5File$exists() raises an H5Lexists error on a multi-level path whose
# parent is absent, so presence is tested one level at a time with names().
make_h5 <- function(path) {
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    f <- hdf5r::H5File$new(path, "r")
    nm <- function(g) tryCatch(if (g == "") names(f) else names(f[[g]]),
                               error = function(e) character(0))
    return(list(
      backend = "hdf5r",
      ls    = nm,
      has   = function(n) {
        p <- strsplit(n, "/")[[1]]
        for (i in seq_along(p)) {
          par <- if (i == 1) "" else paste(p[1:(i - 1)], collapse = "/")
          if (!(p[i] %in% nm(par))) return(FALSE)
        }
        TRUE
      },
      all   = function(n) f[[n]][],
      slice = function(n, a, b) f[[n]][a:b],
      close = function() f$close_all()))
  }
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    L <- rhdf5::h5ls(path)
    full <- gsub("^/+", "", paste0(L$group, "/", L$name))
    return(list(
      backend = "rhdf5",
      ls    = function(g) basename(full[grepl(paste0("^", gsub("^/+", "", g), "/[^/]+$"), full)]),
      has   = function(n) gsub("^/+", "", n) %in% full,
      all   = function(n) rhdf5::h5read(path, n, bit64conversion = "double"),
      slice = function(n, a, b) rhdf5::h5read(path, n, index = list(a:b), bit64conversion = "double"),
      close = function() rhdf5::h5closeAll()))
  }
  stop("Either hdf5r or rhdf5 is required.")
}

# Read one block of cells from a CSR sparse matrix and return, per cell, the
# panel gene values and the total expression sum. The per-cell step is
# vectorised over genes; there is no inner gene loop.
read_block <- function(h5, indptr, a, b, gmap, npanel) {
  s <- indptr[a] + 1; e <- indptr[b + 1]
  ncell <- b - a + 1L
  if (e < s)
    return(list(panel = matrix(0, ncell, npanel), csum = numeric(ncell),
                gi = integer(0), gx = numeric(0), cid = integer(0)))
  gi  <- as.integer(h5$slice("X/indices", s, e)) + 1L
  gx  <- as.numeric(h5$slice("X/data",    s, e))
  p   <- indptr[a:(b + 1)] - indptr[a]
  cid <- rep.int(seq_len(ncell), as.integer(diff(p)))

  csum <- numeric(ncell)
  if (length(gx)) {
    agg <- rowsum(gx, cid, reorder = FALSE)
    csum[as.integer(rownames(agg))] <- as.numeric(agg)
  }
  panel <- matrix(0, ncell, npanel)
  pos <- gmap[gi]; ok <- !is.na(pos)
  if (any(ok)) {
    lin  <- (pos[ok] - 1L) * ncell + cid[ok]   # column-major index into ncell x npanel
    agg2 <- rowsum(gx[ok], lin, reorder = FALSE)
    panel[as.integer(rownames(agg2))] <- as.numeric(agg2)
  }
  list(panel = panel, csum = csum, gi = gi, gx = gx, cid = cid)
}

# CPS slope for every gene at once, via one QR decomposition rather than
# tens of thousands of calls to lm().
slopes_all <- function(G, cps, gm) qr.coef(qr(cbind(1, cps, gm)), G)[2, ]

# Percentile of a statistic within an expression-matched set of genes.
matched_pctl <- function(stat, target_idx, expr, tol = 0.20, min_n = 200) {
  e0  <- expr[target_idx]
  idx <- which(expr >= e0 * (1 - tol) & expr <= e0 * (1 + tol))
  if (length(idx) < min_n) idx <- head(order(abs(expr - e0)), min_n)
  idx <- setdiff(idx, target_idx)
  list(pctl = mean(stat[idx] < stat[target_idx], na.rm = TRUE), n = length(idx))
}

# Mediation at the donor level. Proportion mediated = (b_total - b_direct) / b_total.
# No covariates: this is the specification that reproduces the published table.
mediate <- function(d, outcome, mediator = "MCT4") {
  # unname() matters: coef() carries names such as "Estimate", which would
  # otherwise propagate into the result and break lookup by "prop".
  tot <- unname(coef(summary(lm(as.formula(paste(outcome, "~ cps")), data = d)))["cps", c(1, 4)])
  adj <- unname(coef(summary(lm(as.formula(paste(outcome, "~ cps +", mediator)), data = d)))["cps", c(1, 4)])
  c(b = tot[1], p = tot[2], b_adj = adj[1], p_adj = adj[2],
    prop = 100 * (tot[1] - adj[1]) / tot[1])
}

# Build a synthetic CSR matrix and expose it through the h5 accessor interface,
# so that read_block can be checked against a dense ground truth without a file.
mock_h5 <- function(D) {
  nz <- which(D != 0, arr.ind = TRUE)
  nz <- nz[order(nz[, 1], nz[, 2]), , drop = FALSE]
  store <- list(indices = nz[, 2] - 1L, data = D[nz],
                indptr  = c(0, cumsum(tabulate(nz[, 1], nrow(D)))))
  list(backend = "mock", has = function(n) TRUE,
       ls = function(g) character(0),
       all   = function(n) store[[basename(n)]],
       slice = function(n, a, b) store[[basename(n)]][a:b],
       close = function() invisible(NULL))
}
