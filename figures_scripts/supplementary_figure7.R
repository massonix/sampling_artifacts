#################################################################################
######################### Supplementary Figure 7 ################################
#################################################################################

# Load packages
library(viridis)
library(Seurat)
library(ggpubr)
library(ggthemes)
library(cowplot)
library(ggmap)
library(tidyverse)

# Load ggplots
date <- Sys.Date()
heatmap <- readRDS("current/1-PBMC/results/R_objects/ggplots/heatmap_deg_pbmc.rds")
ridgeplot <- readRDS(file = "current/1-PBMC/results/R_objects/ggplots/ridge_plot_scores.rds")
lollipop <- readRDS("current/1-PBMC/results/R_objects/ggplots/lollipop_go_enrichment_by_celltype_pbmc.rds")

# Change legend position and save
lollipop <- lollipop + 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9.5),
        legend.spacing.x = unit(0.05, 'cm'))
leg <- as_ggplot(get_legend(lollipop))
ggsave(
  filename = "current/doc/figures/legends/lollipop_legend.pdf",
  plot = leg,
  width = 9,
  height = 5,
  units = "cm"
)
lollipop <- lollipop + 
  theme_tufte(ticks = FALSE, base_family = "sans") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,0,1,0), "cm"),
        axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 9))
lollipop

# Arrange
ridge_lollipop <- plot_grid(ridgeplot, lollipop, nrow = 2, ncol = 1, rel_heights = c(0.55, 0.45))
supp7 <- plot_grid(NULL, heatmap[[4]], ridge_lollipop, ncol = 3, nrow = 1, rel_widths = c(0.02, 0.25, 0.75))

# Save
ggsave(plot = supp7, filename = "current/doc/figures/R/supp7.pdf", width = 18, height = 12, units = "cm")


