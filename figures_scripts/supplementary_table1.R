###############################################################################
#######Supplementary table: DEA cold-shock response ###########################
###############################################################################

# Load packages
library(openxlsx)

# Load data
dea_pbmc <- readRDS("current/1-PBMC/results/R_objects/dea_results_pbmc.rds")
dea_pbmc_cell_types <- readRDS("current/1-PBMC/results/R_objects/dea_results_per_cell_type.rds")
dea_cll <- readRDS("current/2-CLL/results/R_objects/dea_results_cll.rds")

# Save as excel
dea_pbmc_cell_types <- purrr::map(dea_pbmc_cell_types, function(df) {
  df <- df[, c("gene", "avg_expr", "avg_logFC", "p_val_adj", "is_significant")]
  names(df) <- c("gene", "average_expression", "log_fc", "p_val_adj", "is_significant")
  df
})
dea <- list(
  PBMC = dea_pbmc, 
  CLL = dea_cll, 
  "T-cell" = dea_pbmc_cell_types$T, 
  NK = dea_pbmc_cell_types$NK, 
  Monocyte = dea_pbmc_cell_types$Monocyte, 
  "B-cell" = dea_pbmc_cell_types$B
)
openxlsx::write.xlsx(dea, "current/doc/tables/supplementary_table1.xlsx")
