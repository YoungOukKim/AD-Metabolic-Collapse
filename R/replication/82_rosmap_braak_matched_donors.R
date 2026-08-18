# =============================================================================
# 82_rosmap_braak_matched_donors.R      (v2)
#
# WHY v2
#   v1 tried to read astrocytes.h5Seurat through SeuratDisk. That package is not
#   installed and is unmaintained. It is also unnecessary: .h5Seurat is an HDF5
#   file, and 47_braak_common_axis.R in this project already reads it directly
#   with rhdf5. This version reuses that reader verbatim, so no new dependency
#   is introduced and the per-donor matrix is built exactly as before.
#
# WHAT THIS IS FOR
#   The manuscript said, of ROSMAP DLPFC:
#       "In the same donors and the same nuclei, however, MCT4 showed no
#        association with Braak stage ... Changing the axis, and nothing else,
#        removes the discrepancy."
#   Table 4 does not support that as written: the pseudotime row is n = 380 and
#   the Braak row is n = 430. Fifty donors differ.
#
#   This script recomputes the ROSMAP Braak slope on EXACTLY the donors that
#   have a pseudotime value. If it stays null, the sentence becomes literally
#   true and we restore it. If not, we report the restricted value instead.
#
#   It changes nothing else. SEA-AD is not touched, pseudotime is not refitted.
#
# GATES (the script stops if one fails, so a wrong input cannot quietly produce
# a plausible answer)
#   G1  the two published ROSMAP values reproduce
#         pseudotime  beta = +0.0214, p = 9.59e-08, n = 380
#         Braak       beta = -0.0007, p = 0.46,     n = 430
#   G2  the pseudotime donor set sits inside the Braak donor set
#   G3  the restricted set has n = 380
#
# HOW TO RUN
#   setwd("D:/work")
#   source("82_rosmap_braak_matched_donors.R")
#
#   Needs only: rhdf5, Matrix, data.table  (all already used by this project)
#   If rhdf5 is missing:
#       if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#       BiocManager::install("rhdf5")
#
#   Output:  D:/work/out_82/82_rosmap_braak_matched.csv
#            D:/work/out_82/rosmap_donor_mct4.csv     (per-donor means, reusable)
#            D:/work/out_82/82_log.txt
#   Send the CSV and the log back.
#
# set.seed(42).
# =============================================================================

suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(data.table)
})
set.seed(42)

# ----------------------------------------------------------------- PATHS ----
WORK <- "D:/work"
ROS  <- file.path(WORK, "ROSMAP", "astrocytes.h5Seurat")
RCL  <- file.path(WORK, "ROSMAP", "ROSMAP_Green_clinical_merged.csv")
OUT  <- file.path(WORK, "out_82")

DRIVER <- "SLC16A3"
BS     <- 5000L                      # block size, as in 47_braak_common_axis.R
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "82_log.txt"), split = TRUE)

cat("=== 82 v2: ROSMAP Braak slope on the pseudotime donor set ===\n")
cat(R.version.string, "|", format(Sys.time()), "\n\n")

if (!file.exists(ROS) || !file.exists(RCL)) {
  cat("A path is wrong. Contents of", file.path(WORK, "ROSMAP"), ":\n  ",
      paste(list.files(file.path(WORK, "ROSMAP")), collapse = " | "), "\n")
  cat("Contents of", WORK, ":\n  ",
      paste(head(list.files(WORK), 60), collapse = " | "), "\n")
  sink(); stop("fix the PATHS block at the top")
}

# --------------------------------------------- 0. show the HDF5 layout ------
cat("HDF5 layout (top levels):\n")
print(h5ls(ROS, recursive = 2))
cat("\n")

# ------------------------------------------------- 1. CLINICAL --------------
CL <- fread(RCL); setnames(CL, tolower(names(CL)))
need <- c("individualid", "pseudotime", "braaksc")
miss <- setdiff(need, names(CL))
if (length(miss)) {
  cat("clinical columns:\n  ", paste(names(CL), collapse = " | "), "\n")
  sink(); stop(paste("missing clinical column(s):", paste(miss, collapse = ", ")))
}
cat("clinical rows:", nrow(CL), "\n")
cat("braaksc: "); print(table(CL$braaksc, useNA = "ifany"))
cat("pseudotime non-missing:", sum(is.finite(CL$pseudotime)), "\n\n")

# --------------------------- 2. PER-DONOR MCT4  (reader from script 47) -----
reuse <- file.path(OUT, "rosmap_donor_mct4.csv")
if (file.exists(reuse)) {
  DM <- fread(reuse); cat("reusing", reuse, "\n")
} else {
  gR  <- as.vector(h5read(ROS, "assays/RNA/features")); nGR <- length(gR)
  stopifnot(DRIVER %in% gR)
  gcol <- match(DRIVER, gR)
  ipR <- as.numeric(h5read(ROS, "assays/RNA/counts/indptr", bit64conversion = "double"))
  fl  <- c("individualID", "nCount_RNA")
  MDR <- as.data.table(lapply(fl, function(f) as.vector(h5read(ROS, paste0("meta.data/", f)))))
  setnames(MDR, fl); MDR[, cell := 1:.N]; nCR <- nrow(MDR)
  cat("nuclei:", nCR, " genes:", nGR, " MCT4 column:", gcol, "\n")

  MDR <- merge(MDR, unique(CL[, .(individualid, pseudotime, braaksc)], by = "individualid"),
               by.x = "individualID", by.y = "individualid", all.x = TRUE, sort = FALSE)
  setorder(MDR, cell)

  kR  <- is.finite(MDR$braaksc) & is.finite(MDR$nCount_RNA) & MDR$nCount_RNA > 0
  dnR <- sort(unique(MDR$individualID[kR])); nDR <- length(dnR)
  diR <- match(MDR$individualID, dnR); kR <- kR & !is.na(diR)
  cat(sprintf("astrocytes used %d | donors with Braak %d\n", sum(kR), nDR))
  ncR <- as.numeric(MDR$nCount_RNA)

  sumD <- numeric(nDR); cntD <- numeric(nDR)
  t0 <- Sys.time()
  for (s0 in seq(1L, nCR, by = BS)) {
    e0 <- min(s0 + BS - 1L, nCR)
    gi <- (s0:e0)[kR[s0:e0]]; if (!length(gi)) next
    sp <- ipR[s0]; cnt <- ipR[e0 + 1L] - sp; if (cnt <= 0) next
    ci <- h5read(ROS, "assays/RNA/counts/indices", start = sp + 1, count = cnt,
                 bit64conversion = "double")
    cd <- as.numeric(h5read(ROS, "assays/RNA/counts/data", start = sp + 1, count = cnt))
    np  <- as.integer(diff(ipR[s0:(e0 + 1L)]))
    cof <- rep(s0:e0, times = np)
    hit <- (as.integer(ci) + 1L) == gcol            # this gene only
    m   <- hit & kR[cof]
    # every kept cell contributes a denominator, even when MCT4 is zero
    d_all <- diR[gi]
    cntD  <- cntD + as.numeric(tabulate(d_all, nbins = nDR))
    if (any(m)) {
      v   <- log1p(cd[m] / ncR[cof[m]] * 1e4)
      agg <- tapply(v, factor(diR[cof[m]], levels = 1:nDR), sum)
      agg[is.na(agg)] <- 0
      sumD <- sumD + as.numeric(agg)
    }
    rm(ci, cd); cat(sprintf("    %d / %d\r", e0, nCR))
  }
  cat(sprintf("\ndone (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  DM <- data.table(donor = dnR, mct4 = sumD / pmax(cntD, 1), n_nuclei = cntD)
  fwrite(DM, reuse); cat("wrote", reuse, "\n")
}

DM <- as.data.table(DM)
DM[, donor := as.character(donor)]
CLd <- unique(CL[, .(donor = as.character(individualid), pseudotime, braaksc)], by = "donor")
D   <- merge(DM, CLd, by = "donor", all.x = TRUE)

set_braak <- D[is.finite(braaksc) & is.finite(mct4)]
set_pseu  <- D[is.finite(pseudotime) & is.finite(mct4)]
cat(sprintf("\nn with Braak %d | n with pseudotime %d\n", nrow(set_braak), nrow(set_pseu)))

# ---------------------------------------------------------------- GATES ----
fit <- function(df, axis) {
  m <- lm(df$mct4 ~ df[[axis]])
  s <- summary(m)$coefficients
  list(beta = unname(s[2, 1]), p = unname(s[2, 4]), n = nrow(df))
}
gP <- fit(set_pseu, "pseudotime"); gB <- fit(set_braak, "braaksc")

cat("\n--- G1: do the published values reproduce? ---\n")
cat(sprintf("  pseudotime  beta % .5f (pub +0.0214)  p %.3g (pub 9.59e-08)  n %d (pub 380)\n",
            gP$beta, gP$p, gP$n))
cat(sprintf("  Braak       beta % .5f (pub -0.0007)  p %.3g (pub 0.46)      n %d (pub 430)\n",
            gB$beta, gB$p, gB$n))
G1 <- abs(gP$beta - 0.0214) <= 5e-4 && abs(gB$beta + 0.0007) <= 5e-4 &&
      gP$n == 380 && gB$n == 430
if (!G1) {
  cat("\nG1 FAILED. Do not use any number below. Send this log back.\n")
  sink(); stop("G1 failed")
}
cat("  G1 passed\n")

cat("\n--- G2: pseudotime donors inside the Braak donors? ---\n")
extra <- setdiff(set_pseu$donor, set_braak$donor)
cat("  with pseudotime but no Braak:", length(extra), "\n")
if (length(extra)) { sink(); stop("G2 failed") }
cat("  G2 passed\n")

matched <- set_braak[donor %in% set_pseu$donor]
cat("\n--- G3: restricted set size ---\n  n =", nrow(matched), "(expected 380)\n")
if (nrow(matched) != 380) { sink(); stop("G3 failed") }
cat("  G3 passed\n")

# ----------------------------------------------- THE NUMBER WE CAME FOR -----
gM <- fit(matched, "braaksc")
cat("\n=== RESULT: ROSMAP Braak slope on the 380 pseudotime donors ===\n")
cat(sprintf("  beta = % .5f    p = %.4g    n = %d\n", gM$beta, gM$p, gM$n))
cat(sprintf("  (all 430 donors: beta = % .5f, p = %.4g)\n", gB$beta, gB$p))
cat(if (gM$p >= 0.05)
      "\n  Still null on the matched donors. 'the same donors ... changing the axis,\n  and nothing else' becomes literally true and the sentence can be restored.\n"
    else
      "\n  NOT null on the matched donors. The sentence stays withdrawn and this\n  restricted value is reported instead.\n")

dropped <- set_braak[!(donor %in% set_pseu$donor)]
cat(sprintf("\n  dropped donors %d | Braak mean %.2f (kept %.2f)\n",
            nrow(dropped), mean(dropped$braaksc), mean(matched$braaksc)))

out <- data.table(
  cohort = "ROSMAP DLPFC", region = "DLPFC",
  axis   = c("pseudotime", "Braak (all donors)", "Braak (pseudotime donors)"),
  donors = c(gP$n, gB$n, gM$n),
  beta   = c(gP$beta, gB$beta, gM$beta),
  p      = c(gP$p, gB$p, gM$p))
fwrite(out, file.path(OUT, "82_rosmap_braak_matched.csv"))
cat("\nwrote", file.path(OUT, "82_rosmap_braak_matched.csv"), "\n")
print(out)
sink()
