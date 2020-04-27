###############################################################################
############################# FIGURE 1 ########################################
###############################################################################

# Load packages
library(viridis)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggmap)
library(ggthemes)
library(readxl)
library(tidyverse)

# Load ggplots
date <- Sys.Date()
tsne_male_pbmc <- readRDS("current/1-PBMC/results/R_objects/ggplots/tsne_time_points_PBMC_Seuratv3.rds")
tsne_male_pbmc <- tsne_male_pbmc$male
umap_cll <- readRDS("current/2-CLL/results/R_objects/ggplots/umap_RT_CLL.rds")
pc1_vs_time <- readRDS("current/1-PBMC/results/R_objects/ggplots/pc1_time_Seuratv3_gg.rds")
dotplot <- readRDS("current/1-PBMC/results/R_objects/ggplots/dotplot_top_genes_signature.rds")
ma_plot <- readRDS("current/6-Reviews/results/R_objects/ma_plot_all_types_pbmc.rds")
vln_cold_score <- readRDS("current/1-PBMC/results/R_objects/ggplots/violin_cold_shock_score.rds")
roc <- readRDS("current/1-PBMC/results/R_objects/ggplots/roc_curve_pbmc.rds")

# Create Violin plot scATAC-seq results (Giovanni's results)
atac_df <- read_excel(path = "current/1-PBMC/data/misc/violin_plot_giovanni.xlsx")
colnames(atac_df) <- c("CLL-DN", "CLL-UP", "Monocyte-DN", "Monocyte-UP", "T cell-DN", "T cell-UP")
atac_df <- atac_df %>% 
  gather(key = "id", value = "scaled_p_value") %>% 
  separate(col = "id", into = c("cell_type", "is_up"), sep = "-") %>% 
  mutate(is_up = factor(is_up), cell_type = factor(cell_type, levels = c("T cell", "Monocyte", "CLL")))
levels(atac_df$is_up) <- c("down", "up")
violin_atac <- atac_df %>%
  ggplot(aes(cell_type, scaled_p_value, fill = is_up, color = is_up)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c("T-cell", "Monocyte", "CLL")) +
  scale_fill_manual("", values = c("dodgerblue3", "firebrick3")) +
  scale_color_manual("", values = c("dodgerblue3", "firebrick3")) +
  labs(x = "", y = "z-score", fill = "", color = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black")) +
  guides(fill = guide_legend(reverse = TRUE))
saveRDS(atac_df, "current/1-PBMC/results/R_objects/atac_df_enhancers.rds")

# Create barplot TFBM enrichment analysis
df <- read_excel("current/doc/tables/supplementary_table1.xlsx", col_names = TRUE, sheet = "Integrative ATAC RNA")
df <- df %>% 
  filter(names %in% c("JUNB", "FOSL1", "Stat6", "IRF9")) %>% 
  dplyr::select("names", "DN peaks", "UP peaks") %>% 
  gather(key = "is_up", value = "log10_p", - "names") %>% 
  dplyr::mutate(is_up = ifelse(is_up == "UP peaks", "up", "down"),
         is_up = factor(is_up, levels = c("up", "down")))
df$names[df$names == "Stat6"] <- "STAT6"
df$names <- factor(df$names, levels = c("JUNB", "FOSL1", "STAT6", "IRF9"))
barplot_atac <- df %>% 
  ggplot(aes(names, log10_p, color = is_up, fill = is_up)) +
  geom_col(position = "dodge") +
  scale_color_manual("", values = c("firebrick3", "dodgerblue3")) +
  scale_fill_manual("", values = c("firebrick3", "dodgerblue3")) +
  labs(x = "", y = "-log10 (p-value)", color = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black", face = "italic"), 
        legend.position = "top",
        plot.margin = unit(c(0, 0.1, 0, 0.5), "cm"))


# Save legends, modify axis, etc
tsne_male_pbmc <- tsne_male_pbmc + 
  theme(legend.position = "top", 
        legend.spacing.x = unit(0.05, 'cm'))
umap_cll <- umap_cll + 
  scale_color_manual(values = c("#999999", "#92e8df", "yellow2", "limegreen", "#632c63", "#e4624e")) +
  theme(legend.position = "top", legend.spacing.y = unit(0.005, 'cm'), legend.spacing.x = unit(0.05, 'cm')) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))
strip_labs <- c("T-cell", "NK", "Monocyte", "B-cell")
names(strip_labs) <- c("T", "NK", "Monocyte", "B")
pc1_vs_time <- pc1_vs_time +
  geom_boxplot(outlier.size = 0.1, fatten = 1) +
  facet_grid(. ~ cell_type, labeller = labeller(cell_type = strip_labs)) +
  scale_x_discrete("Time (h)", labels = as.character(c(0, 2, 8, 24, 48))) +
  theme(strip.text = element_text(size = 10), 
        axis.title.x = element_text(size = 11, face = "plain"),
        axis.title.y = element_text(size = 11, face = "plain"), 
        axis.text.x = element_text(size = 9.5))
ma_plot <- ma_plot + 
  ylab(expression("log"[2]~" (biased / unbiased)")) +
  scale_color_manual("", values = c("gray78", "green4"), labels = c("sig.", "no sig.")) +
  theme(legend.position = "top", legend.spacing.x = unit(0.05, 'cm')) +
  guides(colour = guide_legend(override.aes = list(size=2.5)))
dotplot <- dotplot +
  theme(legend.title = element_blank(), 
        legend.text = element_blank(), 
        legend.spacing.x = unit(0.20, 'cm'))
roc <- roc + 
  geom_line(size = 1) +
  scale_color_manual("", values = c("limegreen", "darkgray"), labels = c("time score", "random"))
  theme(legend.position = "right")
x_titl <- expression("-log"[10]*"(p-value)")
plotlist <- list(tsne_male_pbmc = tsne_male_pbmc, umap_cll = umap_cll, violin_atac = violin_atac, ma_plot = ma_plot, dotplot = dotplot, barplot_atac = barplot_atac, roc = roc)
walk(names(plotlist), function(plt) {
  leg <- as_ggplot(get_legend(plotlist[[plt]]))
  print(plt)
  ggsave(
    filename = str_c("current/doc/figures/legends/", date, "_", plt, ".pdf"), 
    plot = leg, 
    width = 9, 
    height = 5,
    units = "cm"
  )
})
vln_cold_score <- vln_cold_score +
  ylab("Time Score") +
  theme(axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 11, hjust = 0.5),
        axis.text.y = element_text(size = 9.5, color = "grey30"),
        plot.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0.35), "cm"))
tsne_male_pbmc <- tsne_male_pbmc + theme_nothing()
umap_cll <- umap_cll + theme_nothing()
violin_atac <- violin_atac + theme(legend.position = "none")
ma_plot <- ma_plot + theme(legend.position = "none", 
                           axis.title.y = element_text(size = 10), 
                           axis.title.x = element_text(size = 10))
dotplot <- dotplot + 
  scale_x_discrete("Time (h)", labels = as.character(c(0, 2, 8, 24, 48))) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 9, face = "italic", color = "black"),
        axis.title = element_text(size = 12, face = "plain"), 
        axis.text.x = element_text(size = 10),
        plot.margin = unit(c(0.1, 0.2, 0.3, 0), "cm"))
barplot_atac <- barplot_atac + theme(legend.position = "none")
roc <- roc + 
  labs(x = "FPR", y = "TPR") +
  theme(legend.position = "none", 
        axis.text = element_text(size = 9.5))

# Arrange
## Row 1
tsne_umap <- plot_grid(tsne_male_pbmc, NULL, umap_cll, nrow = 1, ncol = 3, rel_widths = c(1, 0.05, 1))
tsne_umap_boxplot <- plot_grid(tsne_umap, NULL, pc1_vs_time, nrow = 3, ncol = 1, rel_heights = c(1, 0.2, 1))
null_violin <- plot_grid(NULL, violin_atac, nrow = 2, ncol = 1, rel_heights = c(1.2, 1))
row1 <- plot_grid(tsne_umap_boxplot, NULL, null_violin, nrow = 1, ncol = 3, rel_widths = c(0.67, 0.025, 0.33))

## Row 2
ma_vln <- ggarrange(plotlist = list(ma_plot, vln_cold_score), nrow = 2, ncol = 1, align = "v")
barplot_roc <- ggarrange(plotlist = list(barplot_atac, roc), nrow = 2, ncol = 1, align = "v")
ma_vln_barplot_roc <- ggarrange(plotlist = list(ma_vln, barplot_roc), ncol = 2, align = "h")
row2 <- plot_grid(dotplot, NULL, ma_vln_barplot_roc, nrow = 1, ncol = 3, rel_widths = c(0.33, 0.005, 0.67))

## Figure 1
fig1 <- plot_grid(row1, NULL, row2, nrow = 3, ncol = 1, rel_heights = c(1.1, 0.065, 1))

# Save
ggsave(plot = fig1, filename = "current/doc/figures/R/figure1.pdf", height = 23, width = 18, units = "cm")

###############################################################################

atac_df <- as.data.frame(atac_df)

# T-cell
atac_t_df <- atac_df[atac_df$cell_type == "T cell", ]
wilcox.test(scaled_p_value ~ is_up, data = atac_t_df, alternative = "two.sided")

# Monocyte
atac_mono_df <- atac_df[atac_df$cell_type == "Monocyte", ]
wilcox.test(scaled_p_value ~ is_up, data = atac_mono_df, alternative = "two.sided")
p_val <- c(0.008765, 0.03602)

# CLL
atac_cll_df <- atac_df[atac_df$cell_type == "CLL", ]
wilcox.test(scaled_p_value ~ is_up, data = atac_cll_df, alternative = "two.sided")


p_val <- c(2.827e-06, 8.821e-08, 0.005228)
names(p_val) <- c("T-cell", "Monocyte", "CLL")
p_val

