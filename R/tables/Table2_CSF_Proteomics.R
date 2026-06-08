suppressMessages(library(car))
# =============================================================================
# Table2_CSF_Proteomics.R
#
# Table 2: CSF V-ATPase subunits and diagnostic group comparison
#          (ADNI Emory TMT-MS), distribution-robust.
#
# Part A: Subunit comparison (Spearman rho primary; raw Pearson in parentheses)
#   - V1A (ATP6V1A) and V1E1 (ATP6V1E1): detection, DX group difference
#   - Spearman rho with MAPT, TREM2, TFRC; vs immunoassay (Elecsys) Tau
#   - Partial rho with LAMP2 / CTSB after Tau control (rank-based)
#   - Inter-subunit Spearman rho (V1A <-> V1E1)
#   - Interpretation: no robust subunit-specific axis
#
# Part B: Diagnostic group changes (ANCOVA, age/sex adjusted; Type II SS)
#   - APP, LCN2, TFRC, MAPT, V1A: KW p, ANCOVA DX p, eta^2, robustness
#
# Input:  Emory TMT-MS CSF proteomics + DXSUM.rda + ADSL.rda
#         (optional) Roche Elecsys immunoassay, SomaScan matrix
# Output: output/tables/Table2_PartA_Subunits.csv
#         output/tables/Table2_PartB_DiagnosticGroups.csv
# =============================================================================
source("R/figures/utils.R")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ── Set your local paths here ─────────────────────────────────────────────────
EMORY_PATH     <- "path/to/emory_results/"
ADNIMERGE_PATH <- "path/to/ADNIMERGE2/data/"
ADNI_AUX_PATH  <- "path/to/adni_aux/"   # Elecsys + SomaScan (optional)
# ─────────────────────────────────────────────────────────────────────────────

em   <- load_adni_proteomics(EMORY_PATH, ADNIMERGE_PATH)
elec <- tryCatch(load_elecsys_tau(ADNI_AUX_PATH), error = function(e) NULL)
if (!is.null(elec)) em <- merge(em, elec, by = "RID", all.x = TRUE)
has_elec <- !is.null(elec) && "ElecsysTau" %in% names(em)

P <- list(
  V1A  = "ATP6V1A_P38606",  V1E1 = "ATP6V1E1_P36543",
  MAPT = "MAPT_P10636",     GFAP = "GFAP_P14136",
  TREM2= "TREM2_Q9NZC2",    TFRC = "TFRC_P02786",
  LAMP2= "LAMP2_P13473",    CTSB = "CTSB_P07858",
  APP  = "APP_P05067",      LCN2 = "LCN2_P80188"
)
P <- P[sapply(P, function(x) x %in% names(em))]

fmt <- function(rho, raw = NULL) {
  if (is.null(raw)) sprintf("%+.2f", rho) else sprintf("%+.2f (%+.2f)", rho, raw)
}

# ── PART A: Subunit comparison (Spearman primary) ─────────────────────────────
cat("\n=== Table 2 Part A: Subunit comparison (Spearman rho primary; raw Pearson in parens) ===\n")
subs    <- c("V1A","V1E1")
targets <- c("MAPT","TREM2","TFRC")
partA <- data.frame()

for (su in subs) {
  if (!su %in% names(P)) next
  n_detect <- sum(!is.na(em[[P[[su]]]]))
  kw <- kruskal.test(em[[P[[su]]]] ~ em$DX)$p.value

  rows <- data.frame(
    Subunit = su,
    Feature = c("Detection rate", "Diagnostic group (DX) difference"),
    Value   = c(sprintf("%d/%d (%d%%)", n_detect, nrow(em), round(100*n_detect/nrow(em))),
                sprintf("KW p = %.3f", kw)))

  # Spearman rho (raw Pearson in parens) for MAPT / TREM2 / TFRC
  for (tg in targets[targets %in% names(P)]) {
    rho <- spearman_r(em[[P[[su]]]], em[[P[[tg]]]])$rho
    raw <- pearson_r (em[[P[[su]]]], em[[P[[tg]]]])$r
    rows <- rbind(rows, data.frame(Subunit = su,
      Feature = sprintf("%s — Spearman \u03c1 (raw)", tg), Value = fmt(rho, raw)))
  }
  # vs immunoassay (Elecsys) Tau
  if (has_elec) {
    rho_e <- spearman_r(em[[P[[su]]]], em$ElecsysTau)$rho
    rows <- rbind(rows, data.frame(Subunit = su,
      Feature = "vs immunoassay (Elecsys) Tau", Value = sprintf("\u03c1 = %+.2f", rho_e)))
  }
  # partial rho | MAPT for LAMP2, CTSB
  for (tg in c("LAMP2","CTSB")) {
    if (!tg %in% names(P)) next
    pr <- spearman_partial(em[[P[[su]]]], em[[P[[tg]]]], em[, P$MAPT, drop = FALSE])$r
    rows <- rbind(rows, data.frame(Subunit = su,
      Feature = sprintf("%s partial \u03c1 (| MAPT)", tg), Value = sprintf("%+.2f", pr)))
  }
  partA <- rbind(partA, rows)
}

# Inter-subunit Spearman (raw Pearson in parens)
if (all(c("V1A","V1E1") %in% names(P))) {
  rho_s <- spearman_r(em[[P$V1A]], em[[P$V1E1]])$rho
  raw_s <- pearson_r (em[[P$V1A]], em[[P$V1E1]])$r
  partA <- rbind(partA, data.frame(
    Subunit = "V1A <-> V1E1", Feature = "Inter-subunit Spearman \u03c1 (raw)",
    Value = fmt(rho_s, raw_s)))
}
partA <- rbind(partA, data.frame(
  Subunit = "Interpretation", Feature = "Subunit-specific axis",
  Value = "No robust subunit-specific axis (Pearson differences not preserved under rank analysis)"))

print(partA, row.names = FALSE)

# ── PART B: Diagnostic group changes (ANCOVA, Type II SS; DX adjusted for AGE, SEX) ─
cat("\n=== Table 2 Part B: Diagnostic Group Changes (ANCOVA, age/sex adjusted) ===\n")
sex_col <- if ("SEX" %in% names(em)) "SEX" else "PTGENDER"
proteins_B <- c("APP","LCN2","TFRC","MAPT","V1A")
partB <- do.call(rbind, lapply(proteins_B[proteins_B %in% names(P)], function(prot) {
  vals <- em[[P[[prot]]]]
  ok   <- !is.na(vals) & !is.na(em$DX) & !is.na(em$AGE) & !is.na(em[[sex_col]])
  sub  <- em[ok, ]
  kw   <- kruskal.test(vals[ok] ~ em$DX[ok])$p.value
  fit  <- lm(as.formula(sprintf("`%s` ~ DX + AGE + %s", P[[prot]], sex_col)), data = sub)
  # ANCOVA Type II SS: DX effect adjusted for AGE and SEX (covariate-adjusted test;
  # method stated in Methods).
  aov_res <- car::Anova(fit, type = 2)               # Type II SS
  dx_p  <- aov_res["DX", "Pr(>F)"]
  eta2  <- aov_res["DX", "Sum Sq"] / sum(aov_res[, "Sum Sq"])
  robust <- if (dx_p < 0.05 && kw < 0.05) "YES"
            else if (dx_p < 0.05 && kw >= 0.05) "Gained"
            else if (dx_p >= 0.05 && kw < 0.05) "Lost"
            else "n.s."
  data.frame(Protein = prot,
             KW_p        = formatC(kw,   format = "f", digits = 3),
             ANCOVA_DX_p = formatC(dx_p, format = "f", digits = 3),
             eta_sq      = round(eta2, 3),
             Robust      = robust, stringsAsFactors = FALSE)
}))
print(partB, row.names = FALSE)

write.csv(partA, "output/tables/Table2_PartA_Subunits.csv",          row.names = FALSE)
write.csv(partB, "output/tables/Table2_PartB_DiagnosticGroups.csv",  row.names = FALSE)
cat("\nSaved:\n  output/tables/Table2_PartA_Subunits.csv\n",
    "  output/tables/Table2_PartB_DiagnosticGroups.csv\n", sep = "")
