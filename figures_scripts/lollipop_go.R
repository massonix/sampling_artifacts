###############################################################################
#######FIGURE 1: Sampling time is a major bias ################################
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


lollipop_go <- readRDS("current/1-PBMC/results/R_objects/ggplots/lollipop_go_enrichment_pbmc.rds")
lollipop_go <- lollipop_go +
  labs(y = x_titl, x = "") +
  scale_x_discrete(labels = labels_lol) +
  theme_tufte(ticks = FALSE, base_family = "sans") +
  theme(axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 9), 
        legend.position = "top",
        legend.text = element_text(size = 9.5),
        legend.spacing.x = unit(0.05, 'cm'))