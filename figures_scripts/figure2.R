###############################################################################
################################ FIGURE 2 #####################################
###############################################################################

# Load packages
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggmap)
library(tidyverse)

# Load ggplots
## Computational correction
tsne_comp <- readRDS("current/1-PBMC/results/R_objects/ggplots/tsnes_original_correction_pbmc.rds")
kbet_comp <- readRDS("current/1-PBMC/results/R_objects/ggplots/barplot_acceptance_rate_correction_pbmc.rds")

## PBMC culturing 
tsne_t_act <- readRDS("current/3-T_cell_activation/results/R_objects/ggplots/tsnes_t_activation.rds")
tsne_t_act <- tsne_t_act[c("0M", "2M")]
kbet_t_act <- readRDS("current/3-T_cell_activation/results/R_objects/ggplots/barplot_kBET_t_act.rds")

## 4ºC storage
date <- Sys.Date()
tsne_4C_pbmc <- readRDS("current/1-PBMC/results/R_objects/ggplots/tsne_4C_pbmc.rds")
tsne_4C_cll <- readRDS("current/2-CLL/results/R_objects/ggplots/tsne_4C_CLL.rds")
kBET_4C_pbmc <- readRDS("current/1-PBMC/results/R_objects/ggplots/barplot_acceptance_rate_storage_pbmc.rds")
kBET_4C_cll <- readRDS("current/2-CLL/results/R_objects/ggplots/barplot_acceptance_rate_storage_cll.rds")

# Generate plot kBET 4ºC
kbet_4C_df <- list(PBMC = kBET_4C_pbmc$data, CLL = kBET_4C_cll$data)
kbet_4C_df <- bind_rows(kbet_4C_df, .id = "dataset")
kbet_4C <- kbet_4C_df %>% 
  mutate(dataset = factor(dataset, levels = c("PBMC", "CLL"))) %>% 
  ggplot(aes(dataset, acceptance_rate, fill = temperature)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_manual(values = c("darkorange1", "#a5cded"),
                    labels = c("21ºC", "4ºC")) +
  scale_y_continuous(limits = c(0, 20)) +
  theme_classic() +
  labs(x = "", y = "kBET (acceptance rate)", fill = "") +
  theme_classic()
  
# Save Legends
tsne_4C_pbmc <- tsne_4C_pbmc +
  scale_color_manual(values = c("#999999", "darkorange1", "#a5cded"),
                     labels = c("Fresh", "21ºC", "4ºC"))
kbet_comp <- kbet_comp +
  scale_y_continuous(limits = c(0, 80)) +
  scale_fill_manual(values = c("gray57", "limegreen"), labels = c("Original", "Corrected"))
  
plotlist <- list(tsne_4C_pbmc, tsne_t_act$`0M`, tsne_comp$original, kbet_4C, kbet_comp)
names(plotlist) <- c("tsne_4C_pbmc", "tsne_t_act", "tsne_comp_corr", "kbet_4C", "kbet_comp")
walk(names(plotlist), function(plt) {
  leg <- as_ggplot(get_legend(plotlist[[plt]]))
  print(plt)
  ggsave(
    filename = str_c("current/doc/figures/legends/", date, "_", plt, ".pdf"), 
    plot = leg, 
    width = 9, 
    height = 5,
    units = "cm"
  )
})

# Remove legends and axis, adjust fontsizes...
tsne_4C_pbmc <- tsne_4C_pbmc + theme_nothing()
tsne_4C_cll <- tsne_4C_cll + theme_nothing()
kbet_4C <- kbet_4C +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        plot.margin = unit(c(0, 0, 0, 1),"cm"))
tsne_t_act <- map(tsne_t_act, function(p) {
  p +
    scale_color_manual("", values = c("#999999", "#632c63", "#e4624e")) +
    theme_nothing()
})
kbet_t_act <- kbet_t_act + 
  theme(axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        strip.text = element_text(size = 11),
        plot.margin = unit(c(0, 0, 0, 1),"cm"))
tsne_comp <- map(tsne_comp, ~ .x + theme_nothing())
kbet_comp <- kbet_comp +
  geom_col(position = "dodge", color = "black") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        plot.margin = unit(c(0, 0, 0, 1),"cm"))

# Arrange
row1 <- plot_grid(tsne_comp$original, tsne_comp$corrected, kbet_comp, ncol = 3, nrow = 1, 
                  rel_widths = c(0.75, 0.75, 1), align = "h", axis = "b")
row2 <- plot_grid(tsne_t_act$`0M`, tsne_t_act$`2M`, kbet_t_act, ncol = 3, nrow = 1, 
                  rel_widths = c(0.75, 0.75, 1), align = "h", axis = "b")
row3 <- plot_grid(tsne_4C_pbmc, tsne_4C_cll, kbet_4C, ncol = 3, nrow = 1, 
                  rel_widths = c(0.75, 0.75, 1), align = "h", axis = "b")
fig2 <- plot_grid(row1, NULL, row2, row3, nrow = 4, ncol = 1, rel_heights = c(1, 0.05, 1, 1))

# Save
ggsave(plot = fig2, filename = "current/doc/figures/R/figure2.pdf", width = 18.5, height = 25, units = "cm")

