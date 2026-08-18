# =============================================================================
# Fig5_CSF_Orthogonal_Test.R
#
# Figure 5 (a-h): an orthogonal test of the energy-starved lysosome at the
#                 (distribution-robust)
#
#   a  -  V1A protein abundance by diagnostic group (Kruskal-Wallis p = 0.999;
#       DEM/CN fold-change = 1.02)  -  group-level preservation
#   b  -  Energy/Demand ratio (mean z[HK1,LDHA,PKM] - z[V1A]) by group (KW p = 0.285)
#   c  -  V1A vs MAPT scatter (untransformed): high but LEVERAGE-DRIVEN Pearson
#       r = +0.86; distribution-robust Spearman rho = +0.25
#   d  -  V1A-Tau robustness ladder:
#       Pearson +0.86 -> trim top-5% +0.23 -> Spearman +0.25
#       -> rho|GFAP +0.05 -> rho|GFAP+TREM2+NfL -0.06 -> vs Elecsys Tau +0.03
#   e  -  HK1 vs MAPT scatter (Spearman rho = +0.52)
#   f  -  HK1-Tau robustness ladder: Spearman +0.52 -> rho|GFAP +0.42
#       -> rho|GFAP+V1A +0.42
#   g  -  Validation against immunoassay (Elecsys) Tau on two proteomic platforms
#       (TMT-MS and SomaScan): HK1 concordant on both; V1A and the SomaScan
#       MAPT / TFRC aptamers are not
#   h  -  V-ATPase subunit comparison (V1A vs V1E1) under rank-based analysis:
#       no robust subunit-specific axis
#
# Input (requires ADNI data access  -  not included in repository):
#   Emory TMT-MS CSF proteomics + DXSUM.rda + ADSL.rda
#   Roche Elecsys immunoassay (UPENNBIOMK_ROCHE_ELECSYS_*.csv)
#   SomaScan 7k post-QC matrix (independent aptamer platform)
#   See loaders in R/figures/utils.R.
#
# Output: output/figures/Fig5_CSF_Orthogonal_Test.png
# =============================================================================
source("R/figures/utils.R")

# -- Set your local paths here -------------------------------------------------
EMORY_PATH     <- "path/to/emory_results/"
ADNIMERGE_PATH <- "path/to/ADNIMERGE2/data/"
ADNI_AUX_PATH  <- "path/to/adni_aux/"   # folder with Elecsys + SomaScan files
# -----------------------------------------------------------------------------

em    <- load_adni_proteomics(EMORY_PATH, ADNIMERGE_PATH)
elec  <- tryCatch(load_elecsys_tau(ADNI_AUX_PATH), error = function(e){message(e$message); NULL})
soma  <- tryCatch(load_somascan(ADNI_AUX_PATH),    error = function(e){message(e$message); NULL})

# -- Protein column mapping ----------------------------------------------------
P <- list(
  V1A  = "ATP6V1A_P38606",   V1E1 = "ATP6V1E1_P36543",
  MAPT = "MAPT_P10636",      GFAP = "GFAP_P14136",
  TREM2= "TREM2_Q9NZC2",     TFRC = "TFRC_P02786",
  LAMP2= "LAMP2_P13473",     CTSB = "CTSB_P07858",
  HK1  = "HK1_P19367",       LDHA = "LDHA_P00338",
  PKM  = "PKM_P14618",       NEFL = "NEFL_P07196"
)
P <- P[sapply(P, function(x) x %in% names(em))]

# Merge immunoassay Tau (for panels d, g) if available
if (!is.null(elec)) em <- merge(em, elec, by = "RID", all.x = TRUE)
has_elec <- !is.null(elec) && "ElecsysTau" %in% names(em)

# -- Panel a: V1A by DX (group-level preservation) -----------------------------
v1a_q <- quantile(em[[P$V1A]], 0.995, na.rm = TRUE)
kw_v1a <- kruskal.test(em[[P$V1A]] ~ em$DX)$p.value
fig5a <- ggplot(em[!is.na(em[[P$V1A]]),], aes(x = DX, y = .data[[P$V1A]], fill = DX)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.75) +
  scale_fill_manual(values = dx_colors) +
  coord_cartesian(ylim = c(NA, v1a_q * 1.05)) +
  annotate("text", x = 2.6, y = v1a_q * 0.98,
           label = sprintf("KW p = %.3f", kw_v1a), size = 3.5, fontface = "italic") +
  labs(x = "Diagnosis", y = "ATP6V1A (V1A)", title = "A") +
  theme_paper + theme(legend.position = "none")

# -- Panel b: Energy/Demand ratio by DX ----------------------------------------
glyc_cols <- unlist(P[intersect(c("HK1","LDHA","PKM"), names(P))])
em$energy_z <- scale(rowMeans(scale(em[, glyc_cols, drop = FALSE]), na.rm = TRUE))[,1]
em$demand_z <- scale(em[[P$V1A]])[,1]
em$ED_diff  <- em$energy_z - em$demand_z
ed_q <- quantile(em$ED_diff, c(0.005, 0.995), na.rm = TRUE)
kw_ed <- kruskal.test(em$ED_diff ~ em$DX)$p.value
fig5b <- ggplot(em[!is.na(em$ED_diff),], aes(x = DX, y = ED_diff, fill = DX)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.75) +
  scale_fill_manual(values = dx_colors) +
  coord_cartesian(ylim = c(ed_q[1]*1.1, ed_q[2]*1.1)) +
  annotate("text", x = 2.6, y = ed_q[2]*0.95,
           label = sprintf("KW p = %.3f", kw_ed), size = 3.5, fontface = "italic") +
  labs(x = "Diagnosis", y = "Energy \u2212 Demand (z-score)", title = "B") +
  theme_paper + theme(legend.position = "none")

# -- Panel c: V1A vs MAPT scatter (untransformed; leverage-driven) -------------
scat   <- em[!is.na(em[[P$V1A]]) & !is.na(em[[P$MAPT]]), ]
r_c    <- pearson_r(scat[[P$V1A]], scat[[P$MAPT]])$r
rho_c  <- spearman_r(scat[[P$V1A]], scat[[P$MAPT]])$rho
xclip  <- quantile(scat[[P$V1A]], 0.99, na.rm = TRUE)
yclip  <- quantile(scat[[P$MAPT]], 0.99, na.rm = TRUE)
fig5c <- ggplot(scat, aes(x = .data[[P$V1A]], y = .data[[P$MAPT]], color = DX)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = dx_colors) +
  coord_cartesian(xlim = c(NA, xclip), ylim = c(NA, yclip)) +
  annotate("text", x = min(scat[[P$V1A]], na.rm = TRUE), y = yclip*0.95,
           label = sprintf("Pearson r = %+.2f (leverage-driven)\nSpearman \u03c1 = %+.2f", r_c, rho_c),
           size = 3.2, hjust = 0, color = "#1f2a52") +
  annotate("text", x = xclip, y = 0, label = "axes capped at 99th pct",
           size = 2.4, hjust = 1, color = "grey60") +
  labs(x = "ATP6V1A (V1A)", y = "MAPT (Tau)", title = "C") +
  theme_paper + theme(legend.position = c(0.85, 0.30))

# -- Panel d: V1A-Tau robustness ladder ----------------------------------------
v1a_trim <- scat[scat[[P$V1A]] < quantile(scat[[P$V1A]], 0.95, na.rm = TRUE), ]
v1a_elec <- if (has_elec) spearman_r(em[[P$V1A]], em$ElecsysTau)$rho else NA
ladder_v1a <- data.frame(
  step = factor(c("Pearson\nraw","trim\ntop-5%","Spearman","\u03c1|GFAP",
                  "\u03c1|GFAP+\nTREM2+NfL","vs Elecsys\nTau"),
                levels = c("Pearson\nraw","trim\ntop-5%","Spearman","\u03c1|GFAP",
                           "\u03c1|GFAP+\nTREM2+NfL","vs Elecsys\nTau")),
  r = c(r_c,
        pearson_r(v1a_trim[[P$V1A]], v1a_trim[[P$MAPT]])$r,
        rho_c,
        spearman_partial(em[[P$V1A]], em[[P$MAPT]], em[, P$GFAP, drop = FALSE])$r,
        spearman_partial(em[[P$V1A]], em[[P$MAPT]],
                         em[, c(P$GFAP, P$TREM2, P$NEFL)])$r,
        v1a_elec),
  kind = c("raw","trim","spearman","partial","partial","elecsys")
)
fig5d <- ggplot(ladder_v1a, aes(x = step, y = r, fill = kind)) +
  geom_col(color = "black", width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", r)),
            vjust = ifelse(ladder_v1a$r >= 0, -0.4, 1.3), size = 3, fontface = "bold") +
  scale_fill_manual(values = c(raw = "#c0392b", trim = "#d68910", spearman = "#d68910",
                               partial = "#9aa0b0", elecsys = "#1f2a52")) +
  ylim(-0.2, 1.0) +
  labs(x = "Robustness control", y = "r (V1A \u2194 MAPT)", title = "D") +
  theme_paper + theme(legend.position = "none",
                      axis.text.x = element_text(size = 7))

# -- Panel e: HK1 vs MAPT scatter ----------------------------------------------
hk_dat <- em[!is.na(em[[P$HK1]]) & !is.na(em[[P$MAPT]]), ]
rho_hk <- spearman_r(hk_dat[[P$HK1]], hk_dat[[P$MAPT]])$rho
hclip  <- quantile(hk_dat[[P$HK1]], 0.99, na.rm = TRUE)
yclip2 <- quantile(hk_dat[[P$MAPT]], 0.99, na.rm = TRUE)
fig5e <- ggplot(hk_dat, aes(x = .data[[P$HK1]], y = .data[[P$MAPT]], color = DX)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = dx_colors) +
  coord_cartesian(xlim = c(NA, hclip), ylim = c(NA, yclip2)) +
  annotate("text", x = min(hk_dat[[P$HK1]], na.rm = TRUE), y = yclip2*0.95,
           label = sprintf("Spearman \u03c1 = %+.2f", rho_hk),
           size = 3.5, hjust = 0, fontface = "bold", color = "#2e8b57") +
  labs(x = "HK1 (hexokinase-1)", y = "MAPT (Tau)", title = "E") +
  theme_paper + theme(legend.position = c(0.85, 0.30))

# -- Panel f: HK1-Tau robustness ladder ----------------------------------------
ladder_hk <- data.frame(
  step = factor(c("Spearman","\u03c1|GFAP","\u03c1|GFAP+V1A"),
                levels = c("Spearman","\u03c1|GFAP","\u03c1|GFAP+V1A")),
  r = c(rho_hk,
        spearman_partial(em[[P$HK1]], em[[P$MAPT]], em[, P$GFAP, drop = FALSE])$r,
        spearman_partial(em[[P$HK1]], em[[P$MAPT]], em[, c(P$GFAP, P$V1A)])$r)
)
fig5f <- ggplot(ladder_hk, aes(x = step, y = r)) +
  geom_col(fill = "#2e8b57", color = "black", width = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", r)), vjust = -0.4, size = 4.8, fontface = "bold") +
  ylim(0, 0.7) +
  labs(x = "Robustness control", y = "r (HK1 \u2194 MAPT)", title = "F") +
  theme_paper + theme(legend.position = "none")

# -- Panel g: concordance with immunoassay (Elecsys) Tau across platforms -----
if (has_elec) {
  prots <- c("V1A","HK1","MAPT","TFRC")
  tmt_r <- sapply(prots, function(p) spearman_r(em[[P[[p]]]], em$ElecsysTau)$rho)
  soma_r <- rep(NA, length(prots)); names(soma_r) <- prots
  if (!is.null(soma)) {
    sm <- merge(em[, c("RID","ElecsysTau")], soma, by = "RID")
    if ("HK1_s"  %in% names(sm)) soma_r["HK1"]  <- spearman_r(sm$HK1_s,  sm$ElecsysTau)$rho
    if ("MAPT_s" %in% names(sm)) soma_r["MAPT"] <- spearman_r(sm$MAPT_s, sm$ElecsysTau)$rho
    if ("TFRC_s" %in% names(sm)) soma_r["TFRC"] <- spearman_r(sm$TFRC_s, sm$ElecsysTau)$rho
  }
  gdat <- rbind(
    data.frame(Protein = prots, Platform = "TMT-MS",   r = as.numeric(tmt_r)),
    data.frame(Protein = prots, Platform = "SomaScan", r = as.numeric(soma_r)))
  gdat$Protein <- factor(gdat$Protein, levels = prots)
  fig5g <- ggplot(gdat, aes(x = Protein, y = r, fill = Platform)) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(r), "", sprintf("%+.2f", r))),
              position = position_dodge(0.8), vjust = ifelse(gdat$r >= 0, -0.4, 1.2),
              size = 2.6, na.rm = TRUE) +
    scale_fill_manual(values = c("TMT-MS" = "#2c6fbf", "SomaScan" = "#d68910")) +
    labs(x = "Validation vs immunoassay Tau", y = "\u03c1 vs Elecsys Tau", title = "G") +
    theme_paper + theme(legend.position = c(0.80, 0.85))
} else {
  fig5g <- ggplot() + annotate("text", x = 0, y = 0,
            label = "Elecsys immunoassay data required for panel g") +
            labs(title = "G") + theme_paper
}

# -- Panel h: subunit comparison V1A vs V1E1 (rank-based) ----------------------
feat   <- c("MAPT","TREM2","TFRC")
hdat <- do.call(rbind, lapply(c("V1A","V1E1"), function(su) {
  zo <- sapply(feat, function(tg) spearman_r(em[[P[[su]]]], em[[P[[tg]]]])$rho)
  pl <- c(LAMP2 = spearman_partial(em[[P[[su]]]], em[[P$LAMP2]], em[, P$MAPT, drop = FALSE])$r,
          CTSB  = spearman_partial(em[[P[[su]]]], em[[P$CTSB]],  em[, P$MAPT, drop = FALSE])$r)
  data.frame(Subunit = su,
             Feature = c("Tau","TREM2","TFRC","LAMP2|Tau","CTSB|Tau"),
             r = c(as.numeric(zo), as.numeric(pl)))
}))
hdat$Feature <- factor(hdat$Feature, levels = c("Tau","TREM2","TFRC","LAMP2|Tau","CTSB|Tau"))
fig5h <- ggplot(hdat, aes(x = Feature, y = r, fill = Subunit)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("V1A" = "#1f2a52", "V1E1" = "#9aa0b0")) +
  labs(x = "Subunit comparison (rank-based)", y = "Spearman \u03c1", title = "H") +
  theme_paper + theme(legend.position = c(0.85, 0.85),
                      axis.text.x = element_text(size = 12))

# -- Combine & save (4 rows x 2 cols, portrait) --------------------------------
fig5 <- (fig5a + fig5b) / (fig5c + fig5d) / (fig5e + fig5f) / (fig5g + fig5h)
save_fig(fig5, "Fig5_CSF_Orthogonal_Test.png", width = 11, height = 15)
