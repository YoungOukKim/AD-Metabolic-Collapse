# =============================================================================
# Table2_CSF_Proteomics.R
#
# Table 2: V-ATPase subunit dissociation and diagnostic group comparison
#          in CSF proteomics (ADNI Emory TMT-MS)
#
# Part A: Subunit Axis Characterization
#   - V1A (neuronal-Tau axis) vs V1E1 (microglial-iron axis)
#   - Detection rates, DX group differences
#   - Correlations: MAPT, GFAP, TREM2, TFRC
#   - Partial r with LAMP2 and CTSB after Tau control
#   - Inter-subunit correlation (V1A ↔ V1E1)
#
# Part B: Diagnostic group changes (ANCOVA, age/sex adjusted)
#   - APP, LCN2, TFRC, MAPT, V1A
#   - KW p, ANCOVA p, eta², robustness
#
# Input:  Emory TMT-MS CSF proteomics + DXSUM.rda + ADSL.rda
#         (requires ADNI data access: https://adni.loni.usc.edu)
# Output: output/tables/Table2_CSF_Proteomics.csv (Part A + Part B)
# =============================================================================
source("R/figures/utils.R")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ── Set your local paths here ─────────────────────────────────────────────────
EMORY_PATH     <- "D:/work/emory_results/"
ADNIMERGE_PATH <- "D:/work/ADNIMERGE2/ADNIMERGE2/data/"
# ─────────────────────────────────────────────────────────────────────────────

em <- load_adni_proteomics(EMORY_PATH, ADNIMERGE_PATH)

# ── Protein column mapping ────────────────────────────────────────────────────
P <- list(
  V1A  = "ATP6V1A_P38606",  V1E1 = "ATP6V1E1_P36543",
  MAPT = "MAPT_P10636",     GFAP = "GFAP_P14136",
  TREM2= "TREM2_Q9NZC2",    TFRC = "TFRC_P02786",
  LAMP2= "LAMP2_P13473",    CTSB = "CTSB_P07858",
  APP  = "APP_P05067",      LCN2 = "LCN2_P80188"
)
P <- P[sapply(P, function(x) x %in% names(em))]

# ── PART A: Subunit axis characterization ────────────────────────────────────
cat("\n=== Table 2 Part A: Subunit Axis Characterization ===\n")

subunits <- c("V1A", "V1E1")
targets  <- c("MAPT","GFAP","TREM2","TFRC")

partA <- data.frame()

for (su in subunits) {
  if (!su %in% names(P)) next
  n_detect <- sum(!is.na(em[[P[[su]]]]))
  pct_detect <- round(n_detect / nrow(em) * 100)

  # KW test vs DX
  kw <- kruskal.test(em[[P[[su]]]] ~ em$DX)

  # Zero-order correlations
  cor_rows <- do.call(rbind, lapply(targets[targets %in% names(P)], function(tg) {
    ct <- cor.test(em[[P[[su]]]], em[[P[[tg]]]], use = "complete.obs")
    data.frame(
      Subunit    = su,
      Feature    = paste0(tg, " correlation"),
      Value      = sprintf("r = %+.3f%s",
                           ct$estimate,
                           ifelse(ct$p.value < 0.001, "***",
                                  ifelse(ct$p.value < 0.01, "**",
                                         ifelse(ct$p.value < 0.05, "*", " (n.s.)"))))
    )
  }))

  # Partial r with LAMP2, CTSB | MAPT
  lyso_rows <- do.call(rbind, lapply(c("LAMP2","CTSB")[c("LAMP2","CTSB") %in% names(P)],
    function(tg) {
      pr <- partial_cor(em[[P[[su]]]], em[[P[[tg]]]], em[[P$MAPT]])
      data.frame(
        Subunit = su,
        Feature = paste0(tg, " partial r (|MAPT)"),
        Value   = sprintf("%+.3f%s", pr$r,
                          ifelse(is.na(pr$p), "",
                                 ifelse(pr$p < 0.001, "***",
                                        ifelse(pr$p < 0.01, "**",
                                               ifelse(pr$p < 0.05, "*", " (n.s.)")))))
      )
  }))

  header <- data.frame(
    Subunit = su,
    Feature = c("Detection rate", "Diagnostic group (DX) difference"),
    Value   = c(sprintf("%d/%d (%d%%)", n_detect, nrow(em), pct_detect),
                sprintf("KW p = %.3f", kw$p.value))
  )

  partA <- rbind(partA, header, cor_rows, lyso_rows)
}

# Inter-subunit correlation
if (all(c("V1A","V1E1") %in% names(P))) {
  n_both <- sum(!is.na(em[[P$V1A]]) & !is.na(em[[P$V1E1]]))
  ct_su  <- cor.test(em[[P$V1A]], em[[P$V1E1]], use = "complete.obs")
  partA  <- rbind(partA, data.frame(
    Subunit = "V1A ↔ V1E1",
    Feature = "Inter-subunit correlation",
    Value   = sprintf("r = %+.3f (%s, n=%d)",
                      ct_su$estimate,
                      ifelse(ct_su$p.value < 0.05, "p < 0.05", "n.s."),
                      n_both)
  ))
}

print(partA, row.names = FALSE)

# ── PART B: Diagnostic group comparison (ANCOVA) ─────────────────────────────
cat("\n=== Table 2 Part B: Diagnostic Group Changes (ANCOVA) ===\n")

proteins_B <- c("APP","LCN2","TFRC","MAPT","V1A")
partB <- do.call(rbind, lapply(proteins_B[proteins_B %in% names(P)], function(prot) {

  vals <- em[[P[[prot]]]]
  ok   <- !is.na(vals) & !is.na(em$DX)
  sub  <- em[ok, ]

  # Kruskal-Wallis
  kw <- kruskal.test(vals[ok] ~ em$DX[ok])

  # ANCOVA with age + sex
  has_age <- "AGE" %in% names(sub)
  has_sex <- "SEX" %in% names(sub) || "PTGENDER" %in% names(sub)
  sex_col <- if ("SEX" %in% names(sub)) "SEX" else "PTGENDER"

  if (has_age && has_sex) {
    formula_str <- paste0(P[[prot]], " ~ DX + AGE + ", sex_col)
    fit  <- tryCatch(lm(as.formula(formula_str), data = sub), error = function(e) NULL)
    if (!is.null(fit)) {
      aov_res <- anova(fit)
      dx_p    <- aov_res["DX", "Pr(>F)"]
      # Partial eta-squared
      ss_dx   <- aov_res["DX", "Sum Sq"]
      ss_total <- sum(aov_res[,"Sum Sq"])
      eta2    <- ss_dx / ss_total
    } else {
      dx_p <- NA; eta2 <- NA
    }
  } else {
    fit     <- lm(as.formula(paste0(P[[prot]], " ~ DX")), data = sub)
    aov_res <- anova(fit)
    dx_p    <- aov_res["DX","Pr(>F)"]
    eta2    <- aov_res["DX","Sum Sq"] / sum(aov_res[,"Sum Sq"])
  }

  robust <- if (is.na(dx_p)) "NA"
  else if (dx_p < 0.05 && kw$p.value < 0.05) "YES"
  else if (dx_p < 0.05 && kw$p.value >= 0.05) "Gained (ANCOVA)"
  else if (dx_p >= 0.05 && kw$p.value < 0.05) "Lost (ANCOVA)"
  else "n.s."

  data.frame(
    Protein    = prot,
    KW_p       = formatC(kw$p.value, format="f", digits=3),
    ANCOVA_DX_p= if(!is.na(dx_p)) formatC(dx_p, format="f", digits=3) else "NA",
    eta_sq     = if(!is.na(eta2)) round(eta2, 3) else NA,
    Robust     = robust,
    stringsAsFactors = FALSE
  )
}))

print(partB, row.names = FALSE)

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(partA, "output/tables/Table2_PartA_SubunitAxis.csv",       row.names=FALSE)
write.csv(partB, "output/tables/Table2_PartB_DiagnosticGroups.csv",  row.names=FALSE)
cat("\nSaved:\n")
cat("  output/tables/Table2_PartA_SubunitAxis.csv\n")
cat("  output/tables/Table2_PartB_DiagnosticGroups.csv\n")
