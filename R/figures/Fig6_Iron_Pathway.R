# =============================================================================
# Fig6_Iron_Pathway.R   (replaces the former Fig6_Iron_Suppression.R)
#
# Figure 6 (a-b): Iron pathway at the protein level
#
#   a — CSF TFRC by diagnostic group: group-level decline after age/sex
#       adjustment (ANCOVA p = 0.036; -4.7% DEM vs CN)
#   b — TFRC-Tau concordance is NOT consistent across platforms or against
#       immunoassay Tau: Elecsys-Tau vs TMT-TFRC rho = +0.08;
#       vs SomaScan-TFRC rho = -0.44  --> individual-level TFRC-Tau coupling
#       is not claimed.
#
# NOTE (revision): the previous version reported a TFRC<->MAPT "suppression
# effect" in which the correlation strengthened after controlling for V1A
# (+0.378 -> +0.602). That analysis used untransformed Pearson correlations
# and an individual-level V1A-Tau structure that did not survive
# distribution-robust analysis; it is no longer claimed. The iron-Tau claim is
# limited to the transcriptomic and group-level proteomic evidence. See CHANGES.md.
#
# Input (requires ADNI data access — not included in repository):
#   Emory TMT-MS CSF proteomics + DXSUM.rda + ADSL.rda
#   Roche Elecsys immunoassay (UPENNBIOMK_ROCHE_ELECSYS_*.csv)
#   SomaScan 7k post-QC matrix (independent aptamer platform)
#
# Output: output/figures/Fig6_Iron_Pathway.png
# =============================================================================
source("R/figures/utils.R")

# ── Set your local paths here ─────────────────────────────────────────────────
EMORY_PATH     <- "D:/work/emory_results/"
ADNIMERGE_PATH <- "D:/work/ADNIMERGE2/ADNIMERGE2/data/"
ADNI_AUX_PATH  <- "D:/work/adni_aux/"   # folder with Elecsys + SomaScan files
# ─────────────────────────────────────────────────────────────────────────────

em   <- load_adni_proteomics(EMORY_PATH, ADNIMERGE_PATH)
elec <- tryCatch(load_elecsys_tau(ADNI_AUX_PATH), error = function(e){message(e$message); NULL})
soma <- tryCatch(load_somascan(ADNI_AUX_PATH),    error = function(e){message(e$message); NULL})

P <- list(TFRC = "TFRC_P02786", MAPT = "MAPT_P10636")
P <- P[sapply(P, function(x) x %in% names(em))]
if (!"TFRC" %in% names(P)) stop("TFRC column not found in Emory dataset.")

# ── Panel a: TFRC by DX (group-level decline) ─────────────────────────────────
# ANCOVA on raw TFRC with age + sex, Type I (sequential) SS with DX entered first
sex_col <- if ("SEX" %in% names(em)) "SEX" else "PTGENDER"
sub <- em[!is.na(em[[P$TFRC]]) & !is.na(em$AGE) & !is.na(em[[sex_col]]), ]
fit_tfrc <- lm(as.formula(sprintf("`%s` ~ DX + AGE + %s", P$TFRC, sex_col)), data = sub)
ancova_p <- anova(fit_tfrc)["DX", "Pr(>F)"]               # Type I, DX first
pct_chg  <- 100 * (median(em[[P$TFRC]][em$DX=="DEM"], na.rm=TRUE) /
                   median(em[[P$TFRC]][em$DX=="CN"],  na.rm=TRUE) - 1)
tfrc_q <- quantile(em[[P$TFRC]], 0.995, na.rm = TRUE)

fig6a <- ggplot(em[!is.na(em[[P$TFRC]]),], aes(x = DX, y = .data[[P$TFRC]], fill = DX)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.75) +
  scale_fill_manual(values = dx_colors) +
  coord_cartesian(ylim = c(NA, tfrc_q * 1.05)) +
  annotate("text", x = 2.55, y = tfrc_q * 0.98,
           label = sprintf("ANCOVA p = %.3f\n%+.1f%% DEM vs CN", ancova_p, pct_chg),
           size = 3.3, fontface = "italic") +
  labs(x = "Diagnosis", y = "TFRC (transferrin receptor)", title = "a") +
  theme_paper + theme(legend.position = "none")

# ── Panel b: TFRC-Tau concordance across platforms (vs immunoassay Tau) ────────
r_tmt <- NA; r_soma <- NA
if (!is.null(elec)) {
  em2   <- merge(em, elec, by = "RID", all.x = TRUE)
  r_tmt <- spearman_r(em2[[P$TFRC]], em2$ElecsysTau)$rho
  if (!is.null(soma)) {
    sm <- merge(elec, soma, by = "RID")
    if ("TFRC_s" %in% names(sm)) r_soma <- spearman_r(sm$TFRC_s, sm$ElecsysTau)$rho
  }
}
bdat <- data.frame(
  comp = factor(c("Elecsys vs\nTMT-TFRC","Elecsys vs\nSomaScan-TFRC"),
                levels = c("Elecsys vs\nTMT-TFRC","Elecsys vs\nSomaScan-TFRC")),
  r = c(r_tmt, r_soma),
  platform = c("TMT-MS","SomaScan"))
fig6b <- ggplot(bdat, aes(x = comp, y = r, fill = platform)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%+.2f", r))),
            vjust = ifelse(bdat$r >= 0, -0.4, 1.3), size = 4, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = c("TMT-MS" = "#2c6fbf", "SomaScan" = "#d68910"), guide = "none") +
  ylim(-0.55, 0.25) +
  labs(x = NULL, y = "\u03c1 (TFRC vs immunoassay Tau)", title = "b") +
  theme_paper

cat(sprintf("TFRC group ANCOVA p = %.3f (%+.1f%% DEM vs CN)\n", ancova_p, pct_chg))
cat(sprintf("TFRC-Tau vs Elecsys: TMT rho = %+.2f ; SomaScan rho = %+.2f\n", r_tmt, r_soma))

fig6 <- fig6a + fig6b
save_fig(fig6, "Fig6_Iron_Pathway.png", width = 11, height = 4.6)
