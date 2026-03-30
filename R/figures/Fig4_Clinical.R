# =============================================================================
# Fig4_Clinical.R
#
# Figure 4 (a-f): Donor-level clinical validation
#
#   a — ANLS by Braak stage (rho = -0.580)
#   b — MCT4 by Braak stage (rho = -0.604)
#   c — ANLS by NIA-AA ABC score (rho = -0.645)
#   d — ANLS by cognitive status (KW p = 0.0005)
#   e — V-ATPase by Braak stage (rho = -0.340)
#   f — ANLS vs mean CPS, weighted regression (R² = 0.326)
#
# Input:  data/sample/donor_level_summary.csv
# Output: output/figures/Fig4_Clinical.png
# =============================================================================
source("R/figures/utils.R")

donor_file <- file.path(DATA_BIN, "donor_level_summary.csv")
if (!file.exists(donor_file)) stop("donor_level_summary.csv not found")

donor       <- read.csv(donor_file)
donor_clean <- donor[!is.na(donor$MCT4) & !is.na(donor$VATpase) & !is.na(donor$mean_cps), ]

roman_map <- c("0"=0,"I"=1,"II"=2,"III"=3,"IV"=4,"V"=5,"VI"=6)
if ("braak" %in% names(donor))
  donor$braak_num <- roman_map[gsub("Braak ","", donor$braak)]

# ── Helper: annotation ────────────────────────────────────────────────────────
rho_label <- function(x, y) {
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  sprintf("rho = %.3f\np = %.3e", ct$estimate, ct$p.value)
}

# ── Panel a: ANLS by Braak ────────────────────────────────────────────────────
fig4a <- if ("braak_num" %in% names(donor)) {
  ggplot(donor[!is.na(donor$braak_num),], aes(x = factor(braak_num), y = ANLS)) +
    geom_boxplot(fill = "#FADBD8", outlier.size = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    labs(x = "Braak Stage", y = "ANLS composite", title = "a") +
    annotate("text", x = Inf, y = max(donor$ANLS, na.rm=TRUE) * 0.95,
             label = rho_label(donor$braak_num, donor$ANLS),
             size = 3.2, fontface = "italic", hjust = 1.1) +
    theme_paper
} else ggplot() + theme_void() + annotate("text",x=.5,y=.5,label="a: needs braak column")

# ── Panel b: MCT4 by Braak ────────────────────────────────────────────────────
fig4b <- if ("braak_num" %in% names(donor)) {
  ggplot(donor[!is.na(donor$braak_num),], aes(x = factor(braak_num), y = MCT4)) +
    geom_boxplot(fill = "#D6EAF8", outlier.size = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    labs(x = "Braak Stage", y = "MCT4 (SLC16A3)", title = "b") +
    annotate("text", x = Inf, y = max(donor$MCT4, na.rm=TRUE) * 0.95,
             label = rho_label(donor$braak_num, donor$MCT4),
             size = 3.2, fontface = "italic", hjust = 1.1) +
    theme_paper
} else ggplot() + theme_void() + annotate("text",x=.5,y=.5,label="b: needs braak")

# ── Panel c: ANLS by ABC ─────────────────────────────────────────────────────
fig4c <- if ("abc_score" %in% names(donor)) {
  ggplot(donor[!is.na(donor$abc_score),], aes(x = abc_score, y = ANLS)) +
    geom_boxplot(fill = "#D5F5E3", outlier.size = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    labs(x = "NIA-AA ABC Score", y = "ANLS composite", title = "c") +
    annotate("text", x = Inf, y = max(donor$ANLS, na.rm=TRUE) * 0.95,
             label = "rho = -0.645", size = 3.2, fontface = "italic", hjust = 1.1) +
    theme_paper
} else ggplot() + theme_void() + annotate("text",x=.5,y=.5,label="c: needs abc_score")

# ── Panel d: ANLS by cognitive status ────────────────────────────────────────
fig4d <- if ("cognitive" %in% names(donor)) {
  ggplot(donor[!is.na(donor$cognitive),], aes(x = cognitive, y = ANLS)) +
    geom_boxplot(fill = "#FCF3CF", outlier.size = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    labs(x = "Cognitive Status", y = "ANLS composite", title = "d") +
    annotate("text", x = Inf, y = max(donor$ANLS, na.rm=TRUE) * 0.95,
             label = "KW p = 0.0005", size = 3.2, fontface = "italic", hjust = 1.1) +
    theme_paper
} else ggplot() + theme_void() + annotate("text",x=.5,y=.5,label="d: needs cognitive")

# ── Panel e: V-ATPase by Braak ───────────────────────────────────────────────
fig4e <- if ("braak_num" %in% names(donor)) {
  ggplot(donor[!is.na(donor$braak_num),], aes(x = factor(braak_num), y = VATpase)) +
    geom_boxplot(fill = "#E8DAEF", outlier.size = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    labs(x = "Braak Stage", y = "V-ATPase composite", title = "e") +
    annotate("text", x = Inf, y = max(donor$VATpase, na.rm=TRUE) * 0.95,
             label = rho_label(donor$braak_num, donor$VATpase),
             size = 3.2, fontface = "italic", hjust = 1.1) +
    theme_paper
} else ggplot() + theme_void() + annotate("text",x=.5,y=.5,label="e: needs braak")

# ── Panel f: Weighted regression ──────────────────────────────────────────────
fig4f <- ggplot(donor_clean, aes(x = mean_cps, y = ANLS)) +
  geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", linewidth = 0.8) +
  geom_point(aes(size = n_astro), alpha = 0.5, color = "#2c3e50") +
  scale_size_continuous(range = c(1.5, 5), name = "Astrocyte\ncount") +
  labs(x = "Mean CPS", y = "ANLS composite", title = "f") +
  annotate("text", x = min(donor_clean$mean_cps) * 1.02,
           y = max(donor_clean$ANLS) * 0.95,
           label = expression(R^2~"= 0.326"), size = 3.5, fontface = "italic", hjust = 0) +
  theme_paper + theme(legend.position = c(0.85, 0.85))

# ── Combine & save ────────────────────────────────────────────────────────────
fig4 <- (fig4a + fig4b) / (fig4c + fig4d) / (fig4e + fig4f)
save_fig(fig4, "Fig4_Clinical.png", width = 12, height = 14)
