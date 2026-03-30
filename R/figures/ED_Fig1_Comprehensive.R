# =============================================================================
# ED_Fig1_Comprehensive.R
#
# Extended Data Figure 1: Comprehensive gene expression changes
#
# Bar chart of % change (Bins 0.2-0.4 -> 0.6-0.8) for all analyzed genes,
# grouped by functional category.
#
# Input:  data/sample/astro_bin_means.csv
# Output: output/figures/ED_Fig1_Comprehensive.png
# =============================================================================
source("R/figures/utils.R")

astro <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
astro <- add_composites_astro(astro)

all_genes <- c("SLC16A3","HK2","LDHA","PDK1","SLC16A1","SLC2A1","PKM","PFKFB3",
               "TFRC","FTH1","FTL","CP","SLC40A1",
               "LAMP1","LAMP2","CTSB","CTSD","LIPA",
               "ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0D1","ATP6V1C1",
               "SLC9A6","MTOR","SOX9","STAT3","LCN2","PTGDS")

cat_map <- c(SLC16A3="ANLS",HK2="ANLS",LDHA="ANLS",PDK1="ANLS",
             SLC16A1="ANLS",SLC2A1="ANLS",PKM="ANLS",PFKFB3="ANLS",
             TFRC="Iron",FTH1="Iron",FTL="Iron",CP="Iron",SLC40A1="Iron",
             LAMP1="Lysosomal",LAMP2="Lysosomal",CTSB="Lysosomal",
             CTSD="Lysosomal",LIPA="Lysosomal",
             ATP6V1A="V-ATPase",ATP6V1B2="V-ATPase",ATP6V0A1="V-ATPase",
             ATP6V0D1="V-ATPase",ATP6V1C1="V-ATPase",
             SLC9A6="Signaling",MTOR="Signaling",SOX9="Signaling",
             STAT3="Signaling",LCN2="Signaling",PTGDS="Signaling")

available <- intersect(all_genes, names(astro))
ed1_df <- do.call(rbind, lapply(available, function(g) {
  early <- mean(astro[[g]][astro$bin %in% c(0.2,0.3,0.4)], na.rm=TRUE)
  late  <- mean(astro[[g]][astro$bin %in% c(0.6,0.7,0.8)], na.rm=TRUE)
  data.frame(Gene = g, pct = (late-early)/early*100, Category = cat_map[g])
}))
ed1_df <- ed1_df[order(ed1_df$pct),]
ed1_df$Gene <- factor(ed1_df$Gene, levels = ed1_df$Gene)

ed_fig1 <- ggplot(ed1_df, aes(x = Gene, y = pct, fill = Category)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) + coord_flip() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("ANLS"="e74c3c","Iron"="#e67e22",
                               "Lysosomal"="#95a5a6","Signaling"="#8e44ad",
                               "V-ATPase"="#3498db")) +
  labs(x = NULL, y = "% Change (Bins 0.2\u20130.4 \u2192 0.6\u20130.8)") +
  theme_paper + theme(legend.position = "right")

save_fig(ed_fig1, "ED_Fig1_Comprehensive.png", width = 8, height = 10)
