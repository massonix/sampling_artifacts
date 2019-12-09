###############################################################################
##################Supplementary 7: TFBS promoters #############################
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

# Create barplot TFBM enrichment analysis
atac_df <- read_excel(path = "current/1-PBMC/data/misc/violin_plot_giovanni.xlsx", sheet = " Ramon_promoters")
colnames(atac_df) <- c("CLL-DN", "CLL-UP", "Monocyte-DN", "Monocyte-UP", "T cell-DN", "T cell-UP")
atac_df <- atac_df %>% 
  gather(key = "id", value = "scaled_p_value") %>% 
  separate(col = "id", into = c("cell_type", "is_up"), sep = "-") %>% 
  mutate(is_up = factor(is_up), cell_type = factor(cell_type, levels = c("T cell", "Monocyte", "CLL")))
levels(atac_df$is_up) <- c("down", "up")
violin_atac <- atac_df %>%
  ggplot(aes(cell_type, scaled_p_value, fill = is_up, color = is_up)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c("T-cell", "Monocyte", "CLL")) +
  scale_fill_manual("", values = c("dodgerblue3", "firebrick3")) +
  scale_color_manual("", values = c("dodgerblue3", "firebrick3")) +
  labs(x = "", y = "Scaled p-value", fill = "", color = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black")) +
  guides(fill = guide_legend(reverse = TRUE), color = FALSE)  

# Save
ggsave(filename = "current/doc/figures/final/supp7.pdf", plot = violin_atac, width = 12, height = 9, units = "cm")
