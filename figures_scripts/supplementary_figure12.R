###############################################################################
############################# Supplementary 12 ################################
###############################################################################

# Load packages
library(tidyverse)
library(ggpubr)
library(cowplot)

# Load ggplot objects
date <- Sys.Date()
tsne <- readRDS("current/1-PBMC/results/R_objects/ggplots/tsne_Smart-seq2_clustered_k2.rds")
barplot <- readRDS("current/1-PBMC/results/R_objects/ggplots/barplot_cluster_distr_Smart-seq2.rds")
volcano <- readRDS("current/1-PBMC/results/R_objects/ggplots/volcano_4ºC_pbmc.rds")
lineplot <- readRDS("current/6-Reviews/results/R_objects/number_deg_4C_vs_21_CLL.rds")

# Move legends
tsne <- tsne +
  labs(x = "tSNE1", y = "tSNE2", color = "") +
  theme_classic() +
  theme(plot.title = element_blank(),
        legend.position = c(0.85, 0.885), 
        legend.background = element_blank(),
        plot.margin = unit(c(0, 1, 0, 0), "cm"))
barplot <- barplot +
  ylab("Number of cells") +
  theme(legend.position = c(0.2, 0.9), 
        axis.text.x = element_text(size = 10),
        legend.background = element_blank())
volcano <- volcano +
  theme(legend.text = element_text(size = 11),
        legend.position = c(0.85, 0.15),
        legend.background = element_blank(),
        plot.margin = unit(c(0, 1, 0, 0), "cm"))
lineplot <- lineplot +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")
lineplot_legend <- as_ggplot(get_legend(lineplot))
ggsave(
  filename = str_c("current/doc/figures/legends/", date, "_", "lineplot_legend_4C", ".pdf"), 
  plot = lineplot_legend, 
  width = 9, 
  height = 5,
  units = "cm"
)
lineplot <- lineplot + theme(legend.position = "none")

# Arrange
supp12 <- ggarrange(plotlist = list(tsne, barplot, volcano, lineplot), nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))

# Save
ggsave(
  plot = supp12, 
  filename = "current/doc/figures/R/supp12.pdf", 
  width = 18.5,
  height = 18, 
  units = "cm"
)
