#################################################################################
########################### Supplementary Figure 1 ##############################
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

# Change themes
tsne_female_pbmc <- tsne_female_pbmc + 
  theme_nothing() +
  theme(plot.margin = unit(c(0, 0, 0, 0.75), "cm"))
tsne_smart <- tsne_smart + 
  theme_nothing() +
  theme(plot.margin = unit(c(0, 0, 0, 0.75), "cm"))

        
# Arrange
plot_list <- list(tsne_female_pbmc, tsne_smart)
supp1 <- ggarrange(plotlist = plot_list, ncol = 2, labels = c("a", "b"))


# Save
ggsave(plot = supp1, filename = "current/doc/figures/R/supp1.pdf", width = 14, height = 8, units = "cm")




