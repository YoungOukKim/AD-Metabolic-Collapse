# =============================================================================
# Fig5_CSF_Validation.R
#
# Figure 5 (a-h): CSF proteomic validation of the energy-starved lysosome
#
#   a — V1A protein abundance by diagnostic group (KW p = 0.999)
#   b — Energy/Demand ratio by diagnostic group (KW p = 0.285)
#   c — V1A vs MAPT scatter colored by diagnosis (r = +0.884)
#   d — Progressive partial correlation bars (GFAP / TREM2 / NfL control)
#   e — HK1 vs MAPT scatter (r = +0.799)
#   f — DX-stratified V1A-MAPT correlations (MCI peak: r = +0.924)
#   g — V1A vs V1E1 correlation heatmap
#   h — Partial r with lysosomal markers after Tau control
#
# Input (requires ADNI data access — not included in repository):
#   Emory TMT-MS CSF proteomics + DXSUM.rda + ADSL.rda
#   Download from https://adni.loni.usc.edu
#   See load_adni_proteomics() in R/figures/utils.R for details.
#
# Output: output/figures/Fig5_CSF_Validation.png
# =============================================================================
source("R/figures/utils.R")

# ── Set your local paths here ─────────────────────────────────────────────────
EMORY_PATH     <- "D:/work/emory_results/"
ADNIMERGE_PATH <- "D:/work/ADNIMERGE2/ADNIMERGE2/data/"
# ─────────────────────────────────────────────────────────────────────────────

em <- load_adni_proteomics(EMORY_PATH, ADNIMERGE_PATH)

# ── Protein column mapping ────────────────────────────────────────────────────
# Column names follow GeneName_UniProtID convention (e.g., "ATP6V1A_P38606")
P <- list(
  V1A  = "ATP6V1A_P38606",   V1E1 = "ATP6V1E1_P36543",
  MAPT = "MAPT_P10636",      GFAP = "GFAP_P14136",
  APP  = "APP_P05067",       TREM2= "TREM2_Q9NZC2",
  LAMP1= "LAMP1_P11279",     LAMP2= "LAMP2_P13473",
  CTSB = "CTSB_P07858",      CTSD = "CTSD_P07339",
  HK1  = "HK1_P19367",       LDHA = "LDHA_P00338",
  LDHB = "LDHB_P07195",      PKM  = "PKM_P14618",
  TFRC = "TFRC_P02786",      LCN2 = "LCN2_P80188",
  NEFL = "NEFL_P07196"
)
P_found   <- P[sapply(P, function(x) x %in% names(em))]
P_missing <- names(P)[!names(P) %in% names(P_found)]
if (length(P_missing) > 0)
  message("Proteins not found in dataset: ", paste(P_missing, collapse=", "))
P <- P_found

# ── Panel a: V1A by DX ────────────────────────────────────────────────────────
v1a_q <- quantile(em[[P$V1A]], 0.995, na.rm=TRUE)
fig5a <- ggplot(em[!is.na(em[[P$V1A]]),],
                aes(x=DX, y=.data[[P$V1A]], fill=DX)) +
  geom_boxplot(outlier.size=0.5, alpha=0.7) +
  scale_fill_manual(values=dx_colors) +
  coord_cartesian(ylim=c(NA, v1a_q*1.05)) +
  annotate("text", x=2, y=v1a_q*0.95,
           label="KW p = 0.999", size=3.5, fontface="italic") +
  labs(x="Diagnosis", y="ATP6V1A (V1A)", title="a") +
  theme_paper + theme(legend.position="none")

# ── Panel b: Energy/Demand ratio ──────────────────────────────────────────────
glyc_cols <- unlist(P[intersect(c("HK1","LDHA","PKM"), names(P))])
em$energy_z <- scale(rowMeans(em[, glyc_cols, drop=FALSE], na.rm=TRUE))[,1]
em$demand_z <- scale(em[[P$V1A]])[,1]
em$ED_diff  <- em$energy_z - em$demand_z
ed_q <- quantile(em$ED_diff, c(0.005,0.995), na.rm=TRUE)

fig5b <- ggplot(em[!is.na(em$ED_diff),], aes(x=DX, y=ED_diff, fill=DX)) +
  geom_boxplot(outlier.size=0.5, alpha=0.7) +
  scale_fill_manual(values=dx_colors) +
  coord_cartesian(ylim=c(ed_q[1]*1.1, ed_q[2]*1.1)) +
  annotate("text", x=2, y=ed_q[2]*0.9,
           label="KW p = 0.285", size=3.5, fontface="italic") +
  labs(x="Diagnosis", y="Energy \u2212 Demand (z-score)", title="b") +
  theme_paper + theme(legend.position="none")

# ── Panel c: V1A vs MAPT scatter ─────────────────────────────────────────────
scat      <- em[!is.na(em[[P$V1A]]) & !is.na(em[[P$MAPT]]),]
r_c       <- cor(scat[[P$V1A]], scat[[P$MAPT]], use="complete.obs")
v1a_clip  <- quantile(scat[[P$V1A]], 0.995, na.rm=TRUE)
mapt_clip <- quantile(scat[[P$MAPT]], 0.995, na.rm=TRUE)

fig5c <- ggplot(scat, aes(x=.data[[P$V1A]], y=.data[[P$MAPT]], color=DX)) +
  geom_point(alpha=0.35, size=1) +
  geom_smooth(aes(group=1), method="lm", se=TRUE, color="black", linewidth=0.6) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.4, linetype="dashed") +
  scale_color_manual(values=dx_colors) +
  coord_cartesian(xlim=c(NA,v1a_clip), ylim=c(NA,mapt_clip)) +
  annotate("text", x=min(scat[[P$V1A]],na.rm=TRUE), y=mapt_clip*0.90,
           label=sprintf("r = +%.3f\nn = %d", r_c, nrow(scat)),
           size=3.5, fontface="bold", hjust=0) +
  labs(x="ATP6V1A", y="MAPT (Tau)", title="c") +
  theme_paper + theme(legend.position=c(0.15,0.75))

# ── Panel d: Partial correlation bars ─────────────────────────────────────────
pc_gfap <- partial_cor(em[[P$V1A]], em[[P$MAPT]], em[[P$GFAP]])
pc_gt   <- partial_cor_double(em[[P$V1A]], em[[P$MAPT]], em[[P$GFAP]], em[[P$TREM2]])
pc_gn   <- partial_cor_double(em[[P$V1A]], em[[P$MAPT]], em[[P$GFAP]], em[[P$NEFL]])
pc_gtn  <- partial_cor_triple(em[[P$V1A]], em[[P$MAPT]],
                               em[[P$GFAP]], em[[P$TREM2]], em[[P$NEFL]])

pc_bar <- data.frame(
  Control = factor(
    c("None","| GFAP","| GFAP+TREM2","| GFAP+NfL","| GFAP+TREM2+NfL"),
    levels = c("None","| GFAP","| GFAP+TREM2","| GFAP+NfL","| GFAP+TREM2+NfL")),
  r = c(r_c, pc_gfap$r, pc_gt$r, pc_gn$r, pc_gtn$r)
)
fig5d <- ggplot(pc_bar, aes(x=Control, y=r, fill=Control)) +
  geom_col(color="black", width=0.6) +
  scale_fill_manual(values=c("None"="#2c3e50","| GFAP"="#2980b9",
                             "| GFAP+TREM2"="#8e44ad","| GFAP+NfL"="#16a085",
                             "| GFAP+TREM2+NfL"="#c0392b")) +
  geom_text(aes(label=sprintf("%.3f",r)), vjust=-0.3, size=3.2, fontface="bold") +
  labs(x="Partial correlation control", y="r (V1A \u2194 MAPT)", title="d") +
  ylim(0,1.0) + theme_paper +
  theme(legend.position="none", axis.text.x=element_text(angle=20, hjust=1))

# ── Panel e: HK1 vs MAPT ─────────────────────────────────────────────────────
hk_dat    <- em[!is.na(em[[P$HK1]]) & !is.na(em[[P$MAPT]]),]
r_hk      <- cor(hk_dat[[P$HK1]], hk_dat[[P$MAPT]], use="complete.obs")
hk_clip   <- quantile(hk_dat[[P$HK1]], 0.995, na.rm=TRUE)
mapt_clip2 <- quantile(hk_dat[[P$MAPT]], 0.995, na.rm=TRUE)

fig5e <- ggplot(hk_dat, aes(x=.data[[P$HK1]], y=.data[[P$MAPT]], color=DX)) +
  geom_point(alpha=0.35, size=1) +
  geom_smooth(aes(group=1), method="lm", se=TRUE, color="black", linewidth=0.6) +
  scale_color_manual(values=dx_colors) +
  coord_cartesian(xlim=c(NA,hk_clip), ylim=c(NA,mapt_clip2)) +
  annotate("text", x=min(hk_dat[[P$HK1]],na.rm=TRUE), y=mapt_clip2*0.90,
           label=sprintf("r = %+.3f\nn = %d", r_hk, nrow(hk_dat)),
           size=3.5, fontface="bold", hjust=0) +
  labs(x="HK1 (Hexokinase 1)", y="MAPT (Tau)", title="e") +
  theme_paper + theme(legend.position=c(0.15,0.65))

# ── Panel f: DX-stratified correlations ──────────────────────────────────────
dx_r <- do.call(rbind, lapply(c("CN","MCI","DEM"), function(dx) {
  sub <- scat[scat$DX == dx,]
  ct  <- cor.test(sub[[P$V1A]], sub[[P$MAPT]])
  data.frame(DX=dx, r=ct$estimate, n=nrow(sub))
}))
dx_r$DX <- factor(dx_r$DX, levels=c("CN","MCI","DEM"))

fig5f <- ggplot(dx_r, aes(x=DX, y=r, fill=DX)) +
  geom_col(color="black", width=0.6) +
  geom_text(aes(label=sprintf("r=%.3f\nn=%d",r,n)),
            vjust=-0.3, size=3, fontface="bold") +
  scale_fill_manual(values=dx_colors) + ylim(0,1.1) +
  labs(x="Diagnosis", y="r (V1A \u2194 MAPT)", title="f") +
  theme_paper + theme(legend.position="none")

# ── Panel g: V1A vs V1E1 heatmap ─────────────────────────────────────────────
targets <- c("MAPT","GFAP","TREM2","TFRC","LAMP2","CTSB")
hm_dat  <- do.call(rbind, lapply(targets[targets %in% names(P)], function(tg) {
  r1 <- cor(em[[P$V1A]],  em[[P[[tg]]]], use="complete.obs")
  r2 <- cor(em[[P$V1E1]], em[[P[[tg]]]], use="complete.obs")
  rbind(data.frame(Target=tg, Subunit="V1A",  r=r1),
        data.frame(Target=tg, Subunit="V1E1", r=r2))
}))
hm_dat$Target <- factor(hm_dat$Target, levels=rev(targets))

fig5g <- ggplot(hm_dat, aes(x=Subunit, y=Target, fill=r)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=sprintf("%+.2f",r)), size=3.5, fontface="bold") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0, limits=c(-0.5,1)) +
  labs(x=NULL, y=NULL, title="g", fill="r") +
  theme_paper + theme(legend.position="right",
                      axis.text.x=element_text(face="bold", size=11))

# ── Panel h: Partial r with lysosomal markers | MAPT ─────────────────────────
lyso_pc <- do.call(rbind, lapply(c("V1A","V1E1"), function(sub) {
  do.call(rbind, lapply(c("LAMP2","CTSB"), function(tg) {
    if (!all(c(sub,tg) %in% names(P))) return(NULL)
    pr <- partial_cor(em[[P[[sub]]]], em[[P[[tg]]]], em[[P$MAPT]])
    data.frame(Subunit=sub, Target=tg, partial_r=pr$r)
  }))
}))
lyso_pc$label <- paste0(lyso_pc$Subunit, "\u2194", lyso_pc$Target)

fig5h <- ggplot(lyso_pc, aes(x=label, y=partial_r, fill=Subunit)) +
  geom_col(color="black", width=0.6) +
  geom_text(aes(label=sprintf("%+.3f",partial_r)),
            vjust=ifelse(lyso_pc$partial_r>0,-0.3,1.3),
            size=3, fontface="bold") +
  geom_hline(yintercept=0, linewidth=0.5) +
  scale_fill_manual(values=c("V1A"="#e74c3c","V1E1"="#3498db")) +
  labs(x=NULL, y="Partial r (| MAPT)", title="h") +
  theme_paper + theme(legend.position="none",
                      axis.text.x=element_text(angle=30, hjust=1))

# ── Combine & save ────────────────────────────────────────────────────────────
fig5 <- (fig5a+fig5b)/(fig5c+fig5d)/(fig5e+fig5f)/(fig5g+fig5h)
save_fig(fig5, "Fig5_CSF_Validation.png", width=14, height=18)
