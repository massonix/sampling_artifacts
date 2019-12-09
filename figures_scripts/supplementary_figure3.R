###############################################################################
#######Supplementary: # detected genes across processing types#################
###############################################################################

# Load packages
library(viridis)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggmap)
library(ggthemes)
library(tidyverse)

# Load R objects
ma_plot <- readRDS("current/2-CLL/results/R_objects/ggplots/ma_plot_all_cell_types_cll.rds")
df_pbmc <- readRDS("current/1-PBMC/results/R_objects/ggplots/dataframe_ngenesVSprocessing_pbmc.rds")
df_cll <- readRDS("current/2-CLL/results/R_objects/ggplots/dataframe_ngenesVSprocessing_cll.rds")
lollipop_go_pbmc <- readRDS("current/1-PBMC/results/R_objects/ggplots/lollipop_go_enrichment_pbmc.rds")
lollipop_go_cll <- readRDS("current/2-CLL/results/R_objects/ggplots/lollipop_go_enrichment_pbmc.rds")

# Create plot 
df <- list(PBMC = df_pbmc, CLL = df_cll) %>% 
  bind_rows(.id = "dataset") %>% 
  dplyr::mutate(processing = factor(processing, levels = c("fresh", "local", "central")),
         dataset = factor(dataset, levels = c("PBMC", "CLL")))
box_plot <- ggplot(df, aes(processing, nFeature_RNA, fill = processing)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_grid(. ~ dataset) +
  scale_fill_manual("", values = c("#999999", "#632c63", "#e4624e")) +
  labs(x = "", y = "# detected genes") +
  scale_y_continuous(limits = c(0, 2250)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(size = 12)) 

# Adjust axis,  font sizes and legends
ma_plot <- ma_plot +
  theme(axis.title = element_text(size = 12),
        legend.position = "none")
labels_lol_pbmc <- rev(c(
  "Arp2/3 \ncomplex-med. \nactin nucleation", 
  "cytolysis", 
  "neg. reg. of \nendopeptidase \n activity", 
  "neg. reg. of \ntranslation"
))
x_titl <- expression("-log"[10]*"(p-value)")
lollipop_go_pbmc <- lollipop_go_pbmc +
  labs(y = x_titl, x = "") +
  scale_x_discrete(labels = labels_lol_pbmc) +
  theme_tufte(ticks = FALSE, base_family = "sans") +
  theme(axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 9), 
        legend.position = c(0.75, 0.25),
        legend.text = element_text(size = 9.5),
        legend.spacing.x = unit(0.05, 'cm'))
labels_lol_cll <- rev(c(
  "Fc receptor\n signaling pathway", 
  "activation of\n immune response", 
  "cellular response to\n oxygen levels", 
  "negative regulation of\n transcription, DNA-templated"
))
lollipop_go_cll <- lollipop_go_cll +
  labs(y = x_titl, x = "") +
  scale_x_discrete(labels = labels_lol_cll) +
  theme_tufte(ticks = FALSE, base_family = "sans") +
  theme(axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 9), 
        legend.position = c(0.75, 0.25),
        legend.text = element_text(size = 9.5),
        legend.spacing.x = unit(0.05, 'cm'))
# Arrange
supp3 <- ggarrange(
  plotlist = list(ma_plot, box_plot, lollipop_go_pbmc, lollipop_go_cll), 
  ncol = 2, 
  nrow = 2,
  # widths = c(1, 1, 0.75),
  labels = c("a", "b", "c", "d")
)

# Save
ggsave(
  plot = supp3, 
  filename = "current/doc/figures/R/supp_figure3.pdf", 
  width = 18.5,
  height = 16, 
  units = "cm"
)

# T tests
# PBMC
pbmc_fVSl <- df_pbmc[df_pbmc$processing %in% c("fresh", "local"), ]
t.test(formula = nFeature_RNA ~ processing, data = pbmc_fVSl, alternative = "two.sided")
pbmc_fVSc <- df_pbmc[df_pbmc$processing %in% c("fresh", "central"), ]
t.test(formula = nFeature_RNA ~ processing, data = pbmc_fVSc, alternative = "two.sided")

# CLL
cll_fVSl <- df_cll[df_cll$processing %in% c("fresh", "local"), ]
t.test(formula = nFeature_RNA ~ processing, data = cll_fVSl, alternative = "two.sided")
cll_fVSc <- df_cll[df_cll$processing %in% c("fresh", "central"), ]
t.test(formula = nFeature_RNA ~ processing, data = cll_fVSc, alternative = "two.sided")















