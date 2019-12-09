#################################################################################
###### Supplementary Figure X: T-cell activation and culturing ##################
#################################################################################

# Load packages
library(Seurat)
library(ggpubr)
library(cowplot)
library(tidyverse)

# Load ggplot object
tsnes_f2_f3 <- readRDS("current/3-T_cell_activation/results/R_objects/ggplots/tsnes_t_activation_rep2.rds")
box_m1 <- readRDS("current/3-T_cell_activation/results/R_objects/ggplots/boxplot_cold_shock_t_act.rds")
box_f2_f3 <- readRDS("current/3-T_cell_activation/results/R_objects/ggplots/boxplot_cold_shock_t_act_rep2.rds")

# Arrange tsnes
# Female2
tsnes_f2 <- tsnes_f2_f3[c("0_rep2_F2", "1_rep2_F2")]
tsnes_f2 <- purrr::map2(tsnes_f2, c("Original", "Cultured"), function(p, x) {
  p +
    ggtitle(x) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
})
tsnes_f2_arranged <- ggarrange(plotlist = tsnes_f2, nrow = 1, ncol = 2)

# Female3
tsnes_f3 <- tsnes_f2_f3[c("0_rep2_F3", "1_rep2_F3")]
tsnes_f3 <- purrr::map2(tsnes_f3, c("Original", "Cultured"), function(p, x) {
  p +
    ggtitle(x) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
})
tsnes_f3_arranged <- ggarrange(plotlist = tsnes_f3, nrow = 1, ncol = 2)
tsnes <- ggarrange(
  plotlist = list(tsnes_f2_arranged, NULL, tsnes_f3_arranged), 
  ncol = 3, 
  nrow = 1, 
  widths = c(1, 0.1, 1)
)

# Arrange boxplots
boxplots <- list(
  male1 = box_m1,
  female2 = box_f2_f3$female2,
  female3 = box_f2_f3$female3
)
legend <- as_ggplot(get_legend(boxplots$male1 + theme(legend.position = "top")))
ggsave(plot = legend, filename = "current/doc/figures/legends/boxplot_t_activation.pdf", height = 9, width = 9, units = "cm")
boxplots <- purrr::map2(boxplots, c("Donor 1", "Donor 2", "Donor 3"), function(p, donor) {
  p +
    ggtitle(donor) +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          legend.text = element_text(size = 10), 
          axis.text.x = element_text(size = 10), 
          axis.title.y = element_text(size = 11),
          plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"),
          legend.position = "none")
})
boxplots_arranged <- ggarrange(plotlist = boxplots, nrow = 1, ncol = 3)

# Final figure
suppX <- ggarrange(
  plotlist = list(tsnes, boxplots_arranged), 
  nrow = 2,
  heights = c(1, 0.75)
)
ggsave(
  filename = "current/doc/figures/R/suppX.pdf", 
  height = 14, 
  width = 18.5, 
  units = "cm"
)
