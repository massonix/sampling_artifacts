#################################################################################
###### Supplementary Figure 5: Cold shock score in T cultured cells #############
#################################################################################

# Load packages
library(Seurat)
library(tidyverse)

# Load ggplot object
vln_cold_t <- readRDS("current/3-T_cell_activation/results/R_objects/ggplots/boxplot_cold_shock_t_act.rds")

# Modify themes
vln_cold_t <- vln_cold_t +
  theme(axis.title.y = element_text(size = 11), 
        axis.text.x = element_text(size = 11), 
        legend.text = element_text(size = 11))

# Save
ggsave(filename = "current/doc/figures/final/suppX.pdf", plot = vln_cold_t, width = 12, height = 9, units = "cm")
