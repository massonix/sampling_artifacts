# This file contains all the function definition used throughout this project
# They are sorted in the exact same order they are used in the analysis notebooks

normalize_hash <- function(sce) {
  # Normalizes each hashtag oligonucleotide (HTO) UMI by dividing by the  
  # geometric mean and log-transforming it.
  # 
  # Args:
  #   sce: A SingleCellExperiment object filtered to only contain the expression 
  #        of HTO.
  #      
  # Returns:
  #   The SCE object with the normalized UMI as a "logcounts" assay.
  
  # Add prior counts to avoid dividing by 0
  logcounts(sce) <- (counts(sce) + 1) %>% 
    apply(1, function(x) log(x / psych::geometric.mean(x + 0.1))) %>% 
    t()
  sce
}

cluster_cells_by_hash <- function(sce, k, dist = "mrw") {
  # Cluster cells by k-medoids based on the UMI of hashtag oligos (HTO).  
  # 
  # Args:
  #   sce: A SingleCellExperiment object normalized and filtered to only contain 
  #        the expression of HTO.
  #   k: Number of desired cluters.
  #   dist: distance metric. Default is Manhattan weighted by rank ("mrw").
  #      
  # Returns:
  #   The SCE object with a new "cluster" variable in the colData, wich contains
  #   the cluster for each cell.
  matr_norm <- logcounts(sce)
  dist_num <- distNumeric(t(matr_norm), t(matr_norm), method = dist)
  clusters <- fastkmed(dist_num, ncluster = k, iterate = 100)
  sce$cluster <- clusters$cluster
  sce
}

plot_heatmap <- function(sce, assay = "logcounts", title, order = FALSE, annotation = TRUE, legend = TRUE) {
  # Plots heatmap to visualize the predicted cell clusters based on hashtag
  # oligonucleotides (HTO) expression
  # 
  # Args:
  #   sce: A SingleCellExperiment object filtered to only contain the expression 
  #        of HTO, and with a "cluster" variable in the colData slot.
  #   assay: the matrix to use for the heatmap (logcounts or classification_mat)
  #   title: character string with the title for the heatmap 
  #   legend: logical indicating if the legend should be plotted. 
  #      
  # Returns:
  #   A pheatmap object
  
  cols <- colorRampPalette(c("#b344bb", "black", "yellow"))(30)
  # Order rows in mat_norm so cells in the same cluster appear together heatmap
  if (assay == "logcounts") {
    mat_norm <- logcounts(sce)
  } else {
    mat_norm <- metadata(sce)$classification_matrix * 1
    sce$cluster <- sce$cell_identity
  }
  barcodes <- colData(sce) %>% 
    as.data.frame() %>% 
    rownames_to_column("barcode")  %>% 
    dplyr::arrange(cluster)
  
  if (order) {
    barcodes <- barcodes %>% 
      dplyr::mutate(cluster = factor(cluster, 
                                     levels = c(rownames(mat_norm), "unassigned"), 
                                     ordered = TRUE)) %>% 
      dplyr::arrange(cluster) 
  }
  mat_norm <- mat_norm[, barcodes$barcode]
  
  # Change cluster labels so the heatmap shows a strong diagonal
  if (order) {
    sce$cluster <- factor(sce$cluster, levels = rownames(mat_norm), ordered = TRUE) 
  }
  
  # Create a cell annotation df to visualize the cluster id in the heatmap
  cell_annot <- data.frame(cluster = barcodes$cluster)          
  rownames(cell_annot) <- barcodes$barcode
  
  # Plot and save the heatmap
  if (annotation) {
    hto_clusters_heat <- pheatmap(
      mat_norm, 
      cluster_rows = FALSE, 
      cluster_cols = FALSE,
      scale = "none",
      color = cols,
      annotation_col = cell_annot,
      labels_col = "",
      legend = legend,
      main = title
    )
  } else {
    hto_clusters_heat <- pheatmap(
      mat_norm, 
      cluster_rows = FALSE, 
      cluster_cols = FALSE,
      scale = "none",
      color = cols,
      labels_col = "",
      legend = legend,
      main = title
    )
  }
  hto_clusters_heat
}

plot_heatmap2 <- function(sce, title) {
  # Plots heatmap to visualize the predicted cell clusters based on hashtag
  # oligonucleotides (HTO) expression
  # 
  # Args:
  #   sce: A SingleCellExperiment object filtered to only contain the expression 
  #        of HTO, and with a "cluster" variable in the colData slot.
  #   assay: the matrix to use for the heatmap (logcounts or classification_mat)
  #      
  # Returns:
  #   A pheatmap object

  cols <- colorRampPalette(c("black", "yellow"))(20)
  
  # Order rows in mat_norm so cells in the same cluster appear together heatmap
  mat_norm <- metadata(sce)$classification_matrix * 1
  sce$cluster <- sce$cell_identity
  barcodes <- colData(sce) %>% 
    as.data.frame() %>% 
    rownames_to_column("barcode")  %>% 
    dplyr::arrange(cluster)
  barcodes <- barcodes %>% 
    dplyr::mutate(cluster = factor(cluster, 
      levels = c(rownames(mat_norm), "unassigned"), 
      ordered = TRUE)) %>% 
    dplyr::arrange(cluster) 
  mat_norm <- mat_norm[, barcodes$barcode]
  
  # Change cluster labels so the heatmap shows a strong diagonal
  sce$cluster <- factor(sce$cluster, levels = rownames(mat_norm), ordered = TRUE) 
  
  # Change row and col labels, add gaps betwen clusters
  cell_annot <- data.frame(cluster = barcodes$cluster)          
  rownames(cell_annot) <- barcodes$barcode
  gap_cols <- match(unique(cell_annot$cluster), cell_annot$cluster)
  row_lab <- str_replace(rownames(mat_norm), pattern = "_", replacement = " ")
  row_lab <- str_replace(row_lab, pattern = "4C", replacement = "4ºC")
  pos_col_lab <- c()
  tab_row <- table(sce$cell_identity)[c(rownames(mat_norm), "unassigned")]
  cell_count <- 0
  for (i in 1:length(tab_row)) {
    pos_col_lab[i] <- cell_count + (tab_row[i] / 2)
    cell_count <- cell_count + tab_row[i]
  }
  pos_col_lab <- round(pos_col_lab)
  col_lab <- rep("", ncol(mat_norm))
  col_lab[pos_col_lab] <- c(row_lab, "unassigned")
  
  # Separate mutiplets and negative assignments
  mat_unassign <- mat_norm[, colnames(sce)[sce$cell_identity == "unassigned"]]
  bar_unassign <- colnames(mat_unassign)
  mutiplet_bar <- c()
  for (lab in rownames(mat_unassign)) {
    sel <- bar_unassign[mat_unassign[lab, ] == 1]
    sel <- sel[!(sel %in% mutiplet_bar)]
    mutiplet_bar <- c(mutiplet_bar, sel)
  }
  negative_bar <- bar_unassign[!(bar_unassign %in% mutiplet_bar)]
  column_selection <- c(
    colnames(mat_norm)[1:(length(colnames(mat_norm)) - length(bar_unassign))], 
    mutiplet_bar, 
    negative_bar
  )
  mat_norm <- mat_norm[, column_selection]
  gap_cols <- c(gap_cols, ncol(mat_norm) - length(negative_bar))
  sum_assign <- sum(sce$cell_identity != "unassigned")
  col_lab[col_lab == "unassigned"]  <- ""
  col_lab[sum_assign + (length(mutiplet_bar) / 2)] <- "Multiplet"
  col_lab[sum_assign + length(mutiplet_bar) + (length(negative_bar) / 2)] <- "Negative"
  
  # Plot and save the heatmap
  hto_clusters_heat <- pheatmap(
    mat_norm, 
    # gaps_col = gap_cols[2:length(gap_cols)],
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    scale = "none",
    color = cols,
    labels_row = row_lab,
    labels_col = col_lab,
    angle_col = 45,
    main = title,
    legend = FALSE
  )
  hto_clusters_heat
}

classify_cells_by_hash  <- function(sce, conditions) {
  # Classify cells by condition based on its hashtag oligo (HTO) expression 
  # 
  # Args:
  #   sce: A SingleCellExperiment object normalized and filtered to only contain 
  #        the expression of HTO. The colData slot needs to contain the cluster
  #        assigned by k-medoids
  #   conditions: character vector that specifies the conditions encoded by HTO
  #      
  # Returns:
  #   The SCE object with new variables in the colData slot:
  #    condition: the condition assigned. In the case of multiplets, they are
  #               separated by a slash (/).
  #    is_assigned: logical variable that indicates whether a cell is assigned
  #                 to one and only one condition.
  #   Vector of thresholds as metadata(sce). 
  mat_umi <- logcounts(sce)
  classification_mat <- matrix(, nrow = 0, ncol = ncol(mat_umi))
  threshs <- c()
  
  for (hto in conditions) {
    # Exclude cells from the cluster in which the HTO is highly expressed
    barcode_sub <- sce[, sce$cluster != hto]
    background <- mat_umi[hto, colnames(barcode_sub)] 
    
    # Exclude highest 0.5% values as potential outliers
    outlier_cutoff <- quantile(background, 0.95)
    back_no_outliers <- background[background < outlier_cutoff]
    
    # Transform back_no_outliers to have positive integers 
    # (required to fit a negative binomial distribution)
    transformed_back <- 1000 * back_no_outliers
    transformed_back2 <- round(transformed_back + abs(min(transformed_back)) + 1)
    
    # Fit a negative binomial distribution 
    fit <- fitdist(transformed_back2, "nbinom", method = "mle") 
    
    # Calculate q=0.99 quantile of the fitted distribution
    hto_cutoff_tr <- quantile(fit, 0.99)$quantiles[1, 1]
    hto_cutoff <- (hto_cutoff_tr - abs(min(transformed_back)) - 1) / 1000
    threshs[hto] <- hto_cutoff
    
    # Use the threshold in (4) to assign each cell (barcode) to an HTO
    hto_classifier <- mat_umi[hto, ] > hto_cutoff
    classification_mat <- rbind(classification_mat, hto_classifier)
  }
  rownames(classification_mat) <- conditions
  metadata(sce)$classification_thresholds <- threshs
  metadata(sce)$classification_matrix <- classification_mat
  
  # Create the cell_identity vector
  cells <- colnames(classification_mat)
  cell_identity <- map_chr(cells, function(cell) {  
    hto_vect <- classification_mat[, cell]
    class <- ifelse(
      sum(hto_vect) == 1, 
      names(hto_vect)[hto_vect], 
      "unassigned"
    )}
  )
  sce$cell_identity <- cell_identity
  sce
}

plot_thresholds <- function(gg, threshs, df, bat, don) {
  # Plot classification thresholds to already created ggplots 
  # 
  # Args:
  #   gg: ggplot with the hashtag oligos (HTO) distributions.
  #   threshs: double vector indicating the thresholds to plot as vertical lines
  #            in the histograms.
  #   df: data frame with the data from which the original distribution was 
  #       created
  #      
  # Returns:
  #   ggplot object with the thresholds plotted.
  df <- filter(df, batch == bat)
  names(threshs) <- str_remove(names(threshs), "_RT")
  names(threshs) <- str_replace(names(threshs), "_4C", " 4ºC")
  if (any(str_detect(names(threshs), "Bioabank"))) {
    names(threshs)[str_detect(names(threshs), "Bioabank")] <-  "24h Bioabank"
  }
  
  for (i in 1:length(threshs)) {
    gg <- gg +
      geom_vline(
        data = df %>% 
          filter(donor == don, hashtag == names(threshs)[i]) %>% 
          mutate(cutoff = threshs[i]),
        aes(xintercept = cutoff),
        color = "darkblue",
        linetype = "dashed"
      )
  }
  gg
}

plot_tsne <- function(sce, 
                      exprs_values, 
                      color_by, 
                      colors,
                      point_size,
                      point_alpha,
                      title = FALSE, 
                      subtitle = FALSE) {
  # Plots a tSNE for a given SingleCellExperiment object
  # 
  # Args:
  #   sce: A SingleCellExperiment object.
  #   exprs_values: either "counts" or "logcounts"
  #   color_by: variable in the colData SCE slot used to color each point.
  #   colors: character vector specifying the colors for the points.
  #   title: string specifying the tSNE title.
  #   subtitle: if specified, adds a subtitle to the tSNE.
  # 
  # Returns:
  #   ggplot object with the resulting tSNE.
  set.seed(1)
  sce <- runTSNE(object = sce, exprs_values = exprs_values)
  tsne <- reducedDim(sce, "TSNE") %>%
    as.data.frame() %>% 
    set_names(c("TSNE1", "TSNE2")) %>% 
    mutate(new_var = colData(sce)[[color_by]]) %>% 
    ggplot(aes(TSNE1, TSNE2, color = new_var)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual("", values = colors) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          legend.text = element_text(size = 11)) 
  
  if (is.character(title)) {
    tsne <- tsne + ggtitle(title)
  } 
  if (is.character(subtitle)) {
    tsne <- tsne + labs(subtitle = subtitle)
  } 

  tsne
}

find_gene_signature <- function(sce, n_genes = 100, random = FALSE) {
  # Find top 100 differentially expressed genes between two conditions.
  # 
  # Args: 
  #   sce: A SingleCellExperiment object with a "label" binary variable
  #        which labels cells as "affected" or "unaffected".
  #   n_genes: numeric indicating how many genes the gene signature should have
  #   random: logical indicating whether or not to find a random singature.
  # 
  # Returns:
  #   A dataframe with 2 variables: "gene" (vector of the top 100 DEG ranked by
  #   p-value) and "sign" (-1 if downregulated, +1 if upregulated). If random = TRUE,
  #   100 extra rows are added with random gene symbols from rownames(sce), and a 
  #   3rd variable "is_random" is added.
  seurat <- Convert(from = sce, to = "seurat")
  seurat <- SetAllIdent(object = seurat, id = "label")
  seurat <- ScaleData(seurat)
  dea_output <- FindMarkers(
    seurat, 
    ident.1 = "affected", 
    ident.2 = "unaffected", 
    test.use = "MAST",
    logfc.threshold = 0
  )
  output_df <- data.frame(
    gene = rownames(dea_output)[1:n_genes],
    sign = sign(dea_output$avg_logFC[1:n_genes]),
    p_val_adj = dea_output$p_val_adj[1:n_genes],
    log_fc = dea_output$avg_logFC[1:n_genes]
  )
  if (random == TRUE) {
    random_signature <- sample(rownames(sce), size = n_genes, replace = FALSE)
    sign_rand <- sign(dea_output[random_signature, "avg_logFC"])
    sign_rand[is.na(sign_rand)] <- 0 
    log_fc_rand <- dea_output[random_signature, "avg_logFC"]
    log_fc_rand[is.na(log_fc_rand)] <- 0
    p_val_adj_rand <- dea_output[random_signature, "p_val_adj"]
    p_val_adj_rand[is.na(p_val_adj_rand)] <- 1
    random_df <- data.frame(
      gene = random_signature, 
      sign = sign_rand, 
      log_fc = log_fc_rand, 
      p_val_adj = p_val_adj_rand
    )
    output_df <- rbind(output_df, random_df)
    output_df$is_random <- c(rep(FALSE, n_genes), rep(TRUE, n_genes))
  } 
  output_df
}

calc_time_score <- function(sce, signature_df, random = FALSE) {
  # Calculate time score for every cell in sce.
  # 
  # Args: 
  #   sce: A SingleCellExperiment object.
  #   signature_df: a dataframe with the variables "gene" (chr vector with the 
  #                 gene symbols of the signature ordered by importance) and 
  #                 "sign" (-1 if downregulated, +1 if upregulated). It can 
  #                 contain a "is_random" to distinguish genes from the random signature.
  #   random: logical indicating whether or not to calculate random time score.
  # 
  # Returns:
  #   Original sce with a new variable "time_score", which is computed by adding together
  #   the weighted scaled counts of the genes in the signature. The weights correspond to 
  #   the inverse ranking in the signature, signed with the direction of the DE. If 
  #   random = TRUE it adds an extra "time_score_random" variable.
  
  row_selection <- rownames(sce) %in% signature_df$gene
  sce_sub <- sce[row_selection, ]
  seurat <- Convert(from = sce_sub, to = "seurat")
  seurat <- ScaleData(seurat)
  
  if ("is_random" %in% colnames(signature_df)) {
    signature_df_sub <- signature_df[!(signature_df$is_random), ]
  } else {
    signature_df_sub <- signature_df
  }
  
  signature_df_sub$weight <- rev(1:nrow(signature_df_sub)) * signature_df_sub$sign
  time_score <- map_dbl(colnames(sce), function(cell) {
    norm_counts <- seurat@scale.data[as.character(signature_df_sub$gene), cell]
    sum(norm_counts * signature_df_sub$weight)
  })
  sce$time_score <- time_score
  
  if (random == TRUE) {
    signature_df_sub <- signature_df[signature_df$is_random, ]
    signature_df_sub$weight <- rev(1:nrow(signature_df_sub)) * signature_df_sub$sign
    time_score_random <- map_dbl(colnames(sce), function(cell) {
      norm_counts <- seurat@scale.data[as.character(signature_df_sub$gene), cell]
      sum(norm_counts * signature_df_sub$weight)
    })
    sce$time_score_random <- time_score_random
  }
  sce
}

calc_time_score_v3 <- function(sce, signature_df, random = FALSE) {
  # Calculate time score for every cell in sce. Works with Seurat 3
  # 
  # Args: 
  #   sce: A SingleCellExperiment object.
  #   signature_df: a dataframe with the variables "gene" (chr vector with the 
  #                 gene symbols of the signature ordered by importance) and 
  #                 "sign" (-1 if downregulated, +1 if upregulated). It can 
  #                 contain a "is_random" to distinguish genes from the random signature.
  #   random: logical indicating whether or not to calculate random time score.
  # 
  # Returns:
  #   Original sce with a new variable "time_score", which is computed by adding together
  #   the weighted scaled counts of the genes in the signature. The weights correspond to 
  #   the inverse ranking in the signature, signed with the direction of the DE. If 
  #   random = TRUE it adds an extra "time_score_random" variable.
  
  row_selection <- rownames(sce) %in% signature_df$gene
  sce_sub <- sce[row_selection, ]
  seurat <- as.Seurat(sce_sub)
  seurat <- ScaleData(seurat)
  
  if ("is_random" %in% colnames(signature_df)) {
    signature_df_sub <- signature_df[!(signature_df$is_random), ]
  } else {
    signature_df_sub <- signature_df
  }
  
  signature_df_sub$weight <- rev(1:nrow(signature_df_sub)) * signature_df_sub$sign
  time_score <- map_dbl(colnames(sce), function(cell) {
    indexes <- match(as.character(signature_df_sub$gene), rownames(seurat))
    indexes <- indexes[!is.na(indexes)]
    norm_counts <- seurat@assays$RNA@scale.data[indexes, cell]
    signature_df_sub_sub <- filter(signature_df_sub, gene %in% names(norm_counts))
    sum(norm_counts * signature_df_sub_sub$weight)
  })
  sce$time_score <- time_score
  
  if (random == TRUE) {
    signature_df_sub <- signature_df[signature_df$is_random, ]
    signature_df_sub$weight <- rev(1:nrow(signature_df_sub)) * signature_df_sub$sign
    time_score_random <- map_dbl(colnames(sce), function(cell) {
      indexes <- match(as.character(signature_df_sub$gene), rownames(seurat))
      indexes <- indexes[!is.na(indexes)]
      norm_counts <- seurat@assays$RNA@scale.data[indexes, cell]
      signature_df_sub_sub <- filter(signature_df_sub, gene %in% names(norm_counts))
      sum(norm_counts * signature_df_sub_sub$weight)
    })
    sce$time_score_random <- time_score_random
  }
  sce
}

test_time_score <- function(sce, random = FALSE, return_ROC = FALSE) {
  # Tests the discriminative power of "time_score" with different accuracies.
  # 
  # Args: 
  #   sce: A SingleCellExperiment object with a "label" variable that contains the observed
  #        classes (affected/unaffected), a time_score and potentially a time_score_random 
  #        variables.
  #   random: logical indicating whether the accuracies should be calculated for the real
  #           time score or for the random one.
  #   return_ROC: logical indicating whether a data frame containing the TPR and FPR for 
  #               several time_score cutoffs should be returned or not.
  # Returns: a data frame with the variables "accuracy_measure" (sensitivity, etc.),
  #          "value". If return_ROC = TRUE, returns another data frame
  #          to plot a ROC curve.
  
  # Define "prediction" and "performance" objects from the ROCR package
  if (random) {
    time_score <- sce$time_score_random
  } else {
    time_score <- sce$time_score
  }
  labs <- factor(sce$label, c("unaffected", "affected"), ordered = TRUE)
  pred <- prediction(predictions = time_score, labels = labs)
  perf <- performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
  
  # Find index optimal cutoff as the one that maximizes specificity and sensitivity
  find_opt_cut <- function(perf, pred){
    cut.ind <- mapply(FUN = function(x, y, p){
      d <- (x - 0)^2 + (y - 1)^2
      ind <- which(d == min(d))
      ind
    }, perf@x.values, perf@y.values, pred@cutoffs)
  }
  ind_opt_cut <- find_opt_cut(perf, pred)
  
  # Compute accuracy metrices
  accuracies <- c("sens", "spec", "acc", "prec")
  values <- c()
  for (metric in accuracies) {
    curr_perf <- performance(pred, metric)
    value <- curr_perf@y.values[[1]][ind_opt_cut]
    values <- c(values, value)
  }
  accuracies_df <- data.frame(accuracies, values)
  
  # Return
  if (return_ROC == TRUE) {
    roc_df <- data.frame(fpr = unlist(perf@x.values), tpr = unlist(perf@y.values))
    list(accuracies = accuracies_df, roc = roc_df)
  } else {
    accuracies_df
  }
}

find_var_genes <- function(sce) {
  # Takes a SingleCellExperiment object, finds and filters for HVG. Return SCE
  fit_var <- trendVar(sce, use.spikes = FALSE) 
  decomp_var <- decomposeVar(sce, fit_var)
  top_hvgs <- order(decomp_var$bio, decreasing = TRUE)
  top_20_pct_hvgs <- top_hvgs[1:(0.2 * length(top_hvgs))]
  sce <- sce[top_20_pct_hvgs, ]
  sce
}

get_GOenrichment <- function(target, universe) {
  # Performs a Gene Ontology enrichment analysis for a given target gene set and
  # universe.
  #
  # Args:
  #   target: character vector with the mgi symbols corresponding to the target
  #           set.
  #   universe: integer vector with the entrez symbols corresponding to the gene
  #             universe. 
  #
  # Returns:
  #   A data.frame with the GO id, the p-value, the Odds score and the 
  #   description of every enriched GO term.
  
  params <- new("GOHyperGParams", geneIds = target, 
                universeGeneIds = universe, annotation = "org.Hs.eg.db",
                ontology = "BP", pvalueCutoff = 1, 
                conditional = TRUE, testDirection = "over")
  hgOver_df <- summary(hyperGTest(params))
  go_results <- hgOver_df %>% 
    arrange(desc(OddsRatio))
  go_results
}

plot_histogram <- function(sce, qc_metric, title, log_scale = FALSE) {
  histogram <- colData(sce) %>% 
    as.data.frame() %>% 
    ggplot(aes_string(qc_metric)) + 
    geom_histogram(bins = 100, alpha = 0.75) +
    xlab(title) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw()
  if (log_scale) {
    histogram + scale_x_log10()
  } else {
    histogram
  }
}

plot_scatter <- function(sce, qc_metrics, alpha, title_x, title_y, log_x = FALSE, log_y = FALSE) {
  scatter <- colData(sce) %>% 
    as.data.frame() %>% 
    ggplot(aes_string(qc_metrics[1], qc_metrics[2])) +
    geom_point(size = 0.8, alpha = alpha) +
    labs(x = title_x, y = title_y) +
    theme_classic2()
  if (log_x & log_y) {
    scatter + scale_x_log10() + scale_y_log10()
  } else if (log_x) {
    scatter + scale_x_log10()
  } else if (log_y) {
    scatter + scale_y_log10()
  } else {
    scatter
  }
}


########################################################################
###############################DEPRECATED###############################
########################################################################

calculate_stress_metric3 <- function(sce, stress_df, cell) {
  # Calculates a stress metric for each cell in sce.
  # 
  # Args:
  #   sce: A SingleCellExperiment object.
  #   stress_df: a data frame that, for each time-affected gene contains the
  #              ensembl gene id ("ensembl"), the gene symbol ("symbol") and
  #              the log fold-change between time-affected and -unaffected cells
  #              (log_FC). Arranged by increasing adjusted p-value.
  #   cell: chr indicating the cell barcode for which to compute the stress metric.
  # 
  # Returns:
  #   the stress metric for cell (as double).
  abs_fold_changes <- stress_df$log_FC[1:100]
  norm_fc <- abs_fold_changes / max(abs_fold_changes) * 10
  logcounts_stress <- assays(sce)$logcounts_gene[, cell]
  stress_metric <- sum(abs_log_fc * logcounts_stress)
  stress_metric
}

perform_dea <- function(sce) {
  # Carries out a differential expression analysis between 2 groups of cells.
  # 
  # Args:
  #   sce: A SingleCellExperiment object, with a "label" binary variable
  # 
  # Returns:
  #   A data frame with the significant differentially expressed genes 
  #   (adj p < 0.01), ranked by log fold-change.
  
  train_sce <- sce[, !(colnames(sce) %in% test_sets[[1]])]
  train_sca <- SceToSingleCellAssay(train_sce)
  zlm_out <- zlm(~ label, train_sca)
  ref <- colnames(attr(zlm_out, "coefC"))[2]
  summary_zlm <- summary(zlm_out, doLRT = ref)
  fit <- summary_zlm$datatable
  p_values <- fit[contrast == ref & component == "H", 
                  .(primerid, `Pr(>Chisq)`)]
  log_FC <- fit[contrast == ref & component == "logFC", 
                .(primerid, coef)]
  fit <- merge(p_values, log_FC, by = "primerid")
  fit <- fit %>% 
    set_colnames(c("gene", "p_value", "log_FC")) %>% 
    filter(!is.na(log_FC)) %>% 
    mutate(adj_p_value = p.adjust(p_value, method = "fdr"),
           log10_adj_p_value = -1 * log10(adj_p_value),
           significance = ifelse(adj_p_value < 0.01, TRUE, FALSE)) %>% 
    arrange(abs(log_FC))
  fit
}

get_stress_signature <- function(sce, test_sce) {
  # Gets gene signature for sce, test it on test_set
  # 
  # Args:
  #   sce: A SingleCellExperiment object
  #   test_set: a SingleCellExperiment object to act as test set.
  # 
  # Returns:
  #   A list, with the genes that composed the gene signature as the first element
  #   and a data frame with the accuracy metrics as the second element.
  
  # Define training set
  train_sce <- sce[, !(colnames(sce) %in% colnames(test_sce))]
  
  # Perform differential expression analysis (DEA) on training set
  train_seurat <- Convert(from = train_sce, to = "seurat")
  train_seurat <- SetAllIdent(object = train_seurat, id = "label")
  train_seurat <- ScaleData(train_seurat, vars.to.regress = "nGene")
  train_de_genes <- FindMarkers(
    train_seurat, 
    ident.1 = "affected", 
    ident.2 = "unaffected", 
    test.use = "MAST",
    logfc.threshold = 0
  )
  signature <- rownames(train_de_genes)[1:100]
  log_fcs <- list(real = train_de_genes$avg_logFC[1:100])
  
  # Get random signature
  random_signature <- sample(rownames(train_sce), size = 100, replace = FALSE)
  log_fcs[["rand"]] <- train_de_genes[random_signature, "avg_logFC"]
  log_fcs[["rand"]][is.na(log_fcs[["rand"]])] <- 0 
  
  # Normalize gene-wise so they are comparable among them (z-score)
  test_sce <- test_sce[union(signature, random_signature), ]
  assays(test_sce)$logcounts_gene <- t(apply(
    logcounts(test_sce), 
    1, 
    function(x) (x - mean(x)) / (sd(x) + 1)
  ))
  
  # Calculate stress metric for each cell. Find optimal cutoff
  out_real_rand <- list(
    real = list(signature = signature), 
    rand = list(signature = random_signature)
  )
  for (sign in c("real", "rand")) {
    curr_sce <- test_sce[out_real_rand[[sign]]$signature, ]
    
    # Stress metric: for each cell, multiply the z-score of the counts of a gene
    # in the signature by the logFC, sum all values.
    stress_metric <- map_dbl(colnames(curr_sce), function(cell) {
      logcounts_stress <- assays(curr_sce)$logcounts_gene[, cell]
      stress_metric <- sum(log_fcs[[sign]] * logcounts_stress)
      stress_metric
    })
    curr_sce$stress_metric <- stress_metric
    pred <- prediction(
      curr_sce$stress_metric, 
      factor(curr_sce$label, ordered = TRUE)
    )
    cost <- performance(pred, "cost")
    cutoff <- pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
    curr_sce$prediction <- ifelse(
      curr_sce$stress_metric > cutoff, 
      "affected", 
      "unaffected"
    )
    out_real_rand[[sign]] <- c(out_real_rand[[sign]], list(
      log_fc = log_fcs,
      cutoff = cutoff,
      stress_metric = curr_sce$stress_metric,
      predictions = curr_sce$prediction
    ))
  }
  
  # Create output list
  output <- list(
    barcodes = colnames(test_sce), 
    labels = test_sce$label, 
    real = out_real_rand$real,
    random = out_real_rand$rand
  )
  output
}

get_stress_signature_seurat <- function(sce, test_sce) {
  # Gets gene signature for sce, test it on test_set using the AddModuleScore
  # function from Seurat.
  # 
  # Args:
  #   sce: A SingleCellExperiment object
  #   test_set: a SingleCellExperiment object to act as test set.
  # 
  # Returns:
  #   A list, with the genes that composed the gene signature as the first element
  #   and a data frame with the accuracy metrics as the second element.
  
  # Define training set
  train_sce <- sce[, !(colnames(sce) %in% colnames(test_sce))]
  
  # Perform differential expression analysis (DEA) on training set
  train_seurat <- Convert(from = train_sce, to = "seurat")
  train_seurat <- SetAllIdent(object = train_seurat, id = "label")
  train_seurat <- ScaleData(train_seurat, vars.to.regress = "nGene")
  train_de_genes <- FindMarkers(
    train_seurat, 
    ident.1 = "affected", 
    ident.2 = "unaffected", 
    test.use = "MAST",
    logfc.threshold = 0
  )
  signature <- rownames(train_de_genes)[1:100]
  
  # Get random signature
  random_signature <- sample(rownames(train_sce), size = 100, replace = FALSE)
  
  # Scale gene-wise so they are comparable among them
  test_seurat <- Convert(from = test_sce, to = "seurat")
  test_seurat <- SetAllIdent(object = test_seurat, id = "label")
  test_seurat <- ScaleData(train_seurat, vars.to.regress = c("nGene", "batch"))
  
  # Compute scores and predictions
  time_random_l <- map(c("time_score", "random_score"), function(score) {
    test_seurat <- AddModuleScore(
      test_seurat, 
      genes.list = signature, 
      ctrl.size = 5, 
      enrich.name = score
    )
    pred <- prediction(
      test_seurat@meta.data[[score]], 
      factor(test_seurat@ident, ordered = TRUE)
    )
    cost <- performance(pred, "cost")
    cutoff <- pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
    prediction <- ifelse(
      test_seurat@meta.data[[score]] > cutoff, 
      "affected", 
      "unaffected"
    )
    list(score = test_seurat@meta.data[[score]], cutoff, prediction)
  })
  names(time_random_l) <- c("time", "random")
  time_random_l$time$signature <- signature
  time_random_l$random$signature <- random_signature
  
  # Create and return output list
  output <- c(
    list(barcodes = rownames(test_seurat@meta.data, labels = test_seurat@ident)),
    time_random_l
  )
  output
}

calc_time_score2 <- function(sce, signature_df, random = FALSE) {
  # Calculate time score for every cell in sce.
  # 
  # Args: 
  #   sce: A SingleCellExperiment object.
  #   signature_df: a dataframe with the variables "gene" (chr vector with the 
  #                 gene symbols of the singature ordered by importance) and 
  #                 "sign" (-1 if downregulated, +1 if upregulated). It can 
  #                 contain a "is_random" to distinguish genes from the random signature.
  #   vars_to_regress: chr vector with the variables to regress out.
  #   random: logical indicating whether or not to calculate random time score.
  # 
  # Returns:
  #   Original sce with a new variable "time_score", which is computed by adding together
  #   the weighted scaled counts of the genes in the signature. The weights correspond to 
  #   the inverse ranking in the signature, signed with the direction of the DE. If 
  #   random = TRUE it adds an extra "time_score_random" variable.
  
  row_selection <- rownames(sce) %in% signature_df$gene
  sce_sub <- sce[row_selection, ]
  assays(sce_sub)$logcounts_gene <- t(apply(
    logcounts(sce_sub), 
    1, 
    function(x) (x - mean(x)) / (sd(x) + 1)
  ))
  
  time_score <- map_dbl(colnames(sce), function(cell) {
    logcounts_stress <- assays(sce)$logcounts_gene[, cell]
    sum(log_fcs[[sign]] * logcounts_stress)
  })
  
  map_dbl(colnames(sce), function(cell) {
    norm_counts <- seurat@scale.data[as.character(signature_df_sub$gene), cell]
    sum(norm_counts * signature_df_sub$log_fc)
  })
  seurat <- Convert(from = sce_sub, to = "seurat")
  seurat <- ScaleData(seurat, vars.to.regress = vars_to_regress)
  
  if ("is_random" %in% colnames(signature_df)) {
    signature_df_sub <- signature_df[!(signature_df$is_random), ]
  } else {
    signature_df_sub <- signature_df
  }
  # Stress metric: for each cell, multiply the z-score of the counts of a gene
  # in the signature by the logFC, sum all values.
  stress_metric <- map_dbl(colnames(curr_sce), function(cell) {
    logcounts_stress <- assays(curr_sce)$logcounts_gene[, cell]
    stress_metric <- sum(log_fcs[[sign]] * logcounts_stress)
    stress_metric
  })
  time_score <- map_dbl(colnames(sce), function(cell) {
    norm_counts <- seurat@scale.data[as.character(signature_df_sub$gene), cell]
    sum(norm_counts * signature_df_sub$log_fc)
  })
  sce$time_score <- time_score
  
  if (random == TRUE) {
    signature_df_sub <- signature_df[signature_df$is_random, ]
    time_score_random <- map_dbl(colnames(sce), function(cell) {
      norm_counts <- seurat@scale.data[as.character(signature_df_sub$gene), cell]
      sum(norm_counts * signature_df_sub$log_fc)
    })
    sce$time_score_random <- time_score_random
  }
  sce
}

find_gene_signature2 <- function(sce, random = FALSE) {
  # Find top 100 differentially expressed genes between two conditions.
  # 
  # Args: 
  #   sce: A SingleCellExperiment object with a "label" binary variable
  #        which labels cells as "affected" or "unaffected".
  #   vars_to_regress: chr vector with the variables to regress out.
  #   random: logical indicating whether or not to find a random singature.
  # 
  # Returns:
  #   A dataframe with 2 variables: "gene" (vector of the top 100 DEG ranked by
  #   p-value) and "sign" (-1 if downregulated, +1 if upregulated). If random = TRUE,
  #   100 extra rows are added with random gene symbols from rownames(sce), and a 
  #   3rd variable "is_random" is added.
  seurat <- Convert(from = sce, to = "seurat")
  seurat <- SetAllIdent(object = seurat, id = "label")
  seurat <- ScaleData(seurat)
  dea_output <- FindMarkers(
    seurat, 
    ident.1 = "affected", 
    ident.2 = "unaffected", 
    test.use = "MAST",
    logfc.threshold = 0
  )
  output_df <- data.frame(
    gene = rownames(dea_output)[1:100],
    log_fc = dea_output$avg_logFC[1:100]
  )
  if (random == TRUE) {
    random_signature <- sample(rownames(sce), size = 100, replace = FALSE)
    log_fc_rand <- dea_output[random_signature, "avg_logFC"]
    log_fc_rand[is.na(log_fc_rand)] <- 0 
    random_df <- data.frame(gene = random_signature, log_fc = log_fc_rand)
    output_df <- rbind(output_df, random_df)
    output_df$is_random <- c(rep(FALSE, 100), rep(TRUE, 100))
  } 
  output_df
}

plot_tsne2 <- function(sce, 
                       color_by, 
                       colors,
                       point_size) {
  # Plots a tSNE for a given SingleCellExperiment object
  # 
  # Args:
  #   sce: A SingleCellExperiment object.
  #   exprs_values: either "counts" or "logcounts"
  #   color_by: variable in the colData SCE slot used to color each point.
  #   colors: character vector specifying the colors for the points.
  #   title: string specifying the tSNE title.
  #   subtitle: if specified, adds a subtitle to the tSNE.
  # 
  # Returns:
  #   ggplot object with the resulting tSNE.
  set.seed(1)
  tsne <- reducedDim(sce, "TSNE") %>%
    as.data.frame() %>% 
    set_names(c("TSNE1", "TSNE2")) %>% 
    mutate(new_var = colData(sce)[[color_by]]) %>% 
    ggplot(aes(TSNE1, TSNE2, color = new_var)) +
    geom_point(size = point_size) +
    scale_color_manual("", values = colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.text = element_text(size = 11)) 
}

pre_process_seurat <- function(object, vars_to_regress = NULL) {
  # Finds HVG, scales the data, runs PCA, tSNE, and UMAP on a seurat object
  object %>%   
    FindVariableFeatures() %>% 
    ScaleData(vars.to.regress = vars_to_regress) %>% 
    RunPCA() %>% 
    RunTSNE(reduction = "pca", dims = 1:15) %>% 
    RunUMAP(reduction = "pca", dims = 1:15) 
}