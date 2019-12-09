###############################################################################
##################Supplementary 5: DEA 4ºC vs fresh ###########################
###############################################################################

# Load packages
library(tidyverse)
library(ggpubr)
library(cowplot)

# Load ggplot objects
tsne <- readRDS("current/1-PBMC/results/R_objects/ggplots/tsne_Smart-seq2_clustered_k2.rds")
barplot <- readRDS("current/1-PBMC/results/R_objects/ggplots/barplot_cluster_distr_Smart-seq2.rds")
volcano <- readRDS("current/1-PBMC/results/R_objects/ggplots/volcano_4ºC_pbmc.rds")
lollipop <- readRDS("current/2-CLL/results/R_objects/ggplots/lollipop_plot_4ºC_cll.rds")

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
lollipop <- lollipop +
  theme(legend.background = element_blank(),
        legend.position = c(0.7, 0.4))

# Arrange
supp5 <- ggarrange(plotlist = list(tsne, barplot, volcano, lollipop), nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))

# Save
ggsave(
  plot = supp5, 
  filename = "current/doc/figures/final/supp5.pdf", 
  width = 18.5,
  height = 18, 
  units = "cm"
)