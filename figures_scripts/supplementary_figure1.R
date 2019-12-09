#################################################################################
##Supplementary Figure 1: Generalizability of bias associated with sampling time#
#################################################################################

# Load packages
library(viridis)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggmap)
library(tidyverse)

# Load ggplots
date <- Sys.Date()
tsne_female_pbmc <- readRDS("current/1-PBMC/results/R_objects/ggplots/tsne_time_points_PBMC_Seuratv3.rds")
tsne_female_pbmc <- tsne_female_pbmc$female
tsne_smart <- readRDS(file = "current/1-PBMC/results/R_objects/ggplots/tsne_time_points_smartseq_gg.rds")
density_smart <- readRDS("current/1-PBMC/results/R_objects/ggplots/explained_variability_smartseq.rds")

# Change themes
tsne_female_pbmc <- tsne_female_pbmc + 
  theme_nothing() +
  theme(plot.margin = unit(c(0, 0, 0, 0.75), "cm"))
tsne_smart <- tsne_smart + 
  theme_nothing() +
  theme(plot.margin = unit(c(0, 0, 0, 0.75), "cm"))
density_smart <- density_smart +
  scale_color_manual("", values = c("indianred2", "steelblue2", "springgreen3"),
                     labels = c("time", "# genes", "library size")) +
  theme(legend.position = c(0.05, 0.95),
        axis.title.x = element_text(size = 11, face = "plain"),
        axis.title.y = element_text(size = 11, face = "plain"),
        axis.text.x = element_text(size = 8),
        legend.text = element_text(size = 9), 
        plot.margin = unit(c(0, 0.1, 0, 0.75), "cm"))
        

# Arrange
plot_list <- list(tsne_female_pbmc, tsne_smart, density_smart)
supp1 <- ggarrange(plotlist = plot_list, ncol = 3, labels = c("a", "b", "c"))


# Save
ggsave(plot = supp1, filename = "current/doc/figures/R/supp1.pdf", width = 18, height = 8, units = "cm")




