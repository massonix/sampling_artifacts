#################################################################################
################ Supplementary Figure 8: Annotation #############################
#################################################################################

# Load packages
library(viridis)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggmap)
library(ggthemes)
library(readxl)
library(tidyverse)

# Load R objects
date <- Sys.Date()
pbmc <- readRDS("current/1-PBMC/results/R_objects/pbmc_Seurat_male.rds")
cll <- readRDS("current/2-CLL/results/R_objects/cll_seurat_fig1.rds")
markers <- readRDS("current/2-CLL/results/R_objects/markers_cll_clusters.rds")
t_act <- readRDS("current/3-T_cell_activation/results/R_objects/t_act_seurat_2M.rds")

# PBMC
Idents(pbmc) <- "cell_type"
levels(Idents(pbmc)) <- c("T-cell", "NK", "Monocyte", "B-cell")
palette2 <- c("#c20a35", "#aa2edc", "#71bdd0", "#bbaa2a")
tsne_cell_type <- DimPlot(pbmc, reduction = "tsne", cols = palette2)
legend1 <- as_ggplot(get_legend(tsne_cell_type))
ggsave(
  filename = str_c("current/doc/figures/legends/", date, "_", "annotation_pbmc1", ".pdf"), 
  plot = legend1, 
  width = 9, 
  height = 5,
  units = "cm"
)
legend2 <- as_ggplot(get_legend(FeaturePlot(pbmc, features = "IL7R")))
ggsave(
  filename = str_c("current/doc/figures/legends/", date, "_", "annotation_pbmc2", ".pdf"), 
  plot = legend2, 
  width = 9, 
  height = 5,
  units = "cm"
)
tsne_cell_type <- tsne_cell_type + theme_nothing()
pbmc_markers <- c("IL7R", "GNLY", "LYZ", "MS4A1")
feature_plots <- purrr::map(pbmc_markers, function(x){
  p <- FeaturePlot(pbmc, features = x, reduction = "tsne")
  p +
    theme(axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
})
feature_plots_arr <- ggarrange(plotlist = feature_plots, nrow = 2, ncol = 2)
row_pbmc <- ggarrange(plotlist = list(tsne_cell_type, feature_plots_arr), ncol = 2, widths = c(1, 1))

# CLL
Idents(cll) <- "donor"
umap_donor <- DimPlot(cll, reduction = "umap")
legend1 <- as_ggplot(get_legend(umap_donor))
ggsave(
  filename = str_c("current/doc/figures/legends/", date, "_", "annotation_cll", ".pdf"), 
  plot = legend1, 
  width = 9, 
  height = 5,
  units = "cm"
)
umap_donor <- umap_donor + theme_nothing()
cll_markers <- c("IGLC2", "IGHG1", "IGHA1", "IL7R")
feature_plots <- purrr::map(cll_markers, function(x){
  p <- FeaturePlot(cll, features = x, reduction = "umap")
  p +
    theme(axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
})
feature_plots_arr <- ggarrange(plotlist = feature_plots, nrow = 2, ncol = 2)
row_cll <- ggarrange(plotlist = list(umap_donor, feature_plots_arr), ncol = 2, widths = c(1, 1))

# T-cell activation
Idents(t_act) <- "cell_type"
levels(Idents(t_act)) <- c("CD4 T-cell", "Cytotoxic", "Cycling", "B-cell")
tsne_cell_type <- DimPlot(t_act, reduction = "tsne", cols = c("#c20a35", "#aa2edc", "chartreuse3", "#bbaa2a"))
legend1 <- as_ggplot(get_legend(tsne_cell_type))
ggsave(
  filename = str_c("current/doc/figures/legends/", date, "_", "annotation_t_act", ".pdf"), 
  plot = legend1, 
  width = 9, 
  height = 5,
  units = "cm"
)
tsne_cell_type <- tsne_cell_type + theme_nothing()
t_act_markers <- c("IL7R", "NKG7", "TOP2A", "MS4A1")
feature_plots <- purrr::map(t_act_markers, function(x){
  p <- FeaturePlot(t_act, features = x, reduction = "tsne")
  p +
    theme(axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
})
feature_plots_arr <- ggarrange(plotlist = feature_plots, nrow = 2, ncol = 2)
row_t_act <- ggarrange(plotlist = list(tsne_cell_type, feature_plots_arr), ncol = 2, widths = c(1, 1))

# Arrange
supp8 <- plot_grid(row_pbmc, NULL, row_cll, NULL, row_t_act, nrow = 5, ncol = 1, rel_heights = c(1, 0.1, 1, 0.1, 1))

# Save
ggsave(filename = "current/doc/figures/R/supp8.pdf", height = 27, width = 18, units = "cm")


