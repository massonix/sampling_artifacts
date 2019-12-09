#################################################################################
##############Supplementary Figure 4: Computational correction ##################
#################################################################################

# Load packages
library(viridis)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggmap)
library(tidyverse)

# Load ggplots
umaps_cll <- readRDS("current/2-CLL/results/R_objects/ggplots/umaps_original_corrected.rds")
bootstrap <- readRDS("current/1-PBMC/results/R_objects/bootstrap_sil_width_pct_cells.rds")

# Modify themes and save legends
umaps_cll <- map(umaps_cll, ~.x + theme_nothing())
bootstrap <- bootstrap +
  scale_fill_manual("", values = c("azure3", "chartreuse2"), labels = c("original", "corrected")) +
  scale_color_manual("", values = c("azure3", "chartreuse2"), labels = c("original", "corrected")) +
  xlab("% biased cells") +
  theme_classic() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.position = c(0.9, 0.9), 
        plot.margin = unit(c(0, 1, 0, 1), "cm"))

# Arrange
supp4 <- ggarrange(
  plotlist = list(umaps_cll$original, umaps_cll$regressed, bootstrap), 
  ncol = 3, 
  labels = c("a", "", "b"),
  widths = c(0.75, 0.75, 1))
supp4

# Save
ggsave(plot = supp4, filename = "current/doc/figures/R/supp4.pdf", width = 18, height = 8, units = "cm")



