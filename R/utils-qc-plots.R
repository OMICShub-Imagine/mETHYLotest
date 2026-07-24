#' (Internal) Plots PCA and t-SNE (Simplified 2D version)
#'
#' @description
#' Internal utility function. Does not compute anything, formats
#' ggplot graphics from pre-calculated results.
#'
#' @param pca_result The pre-calculated `prcomp` object.
#' @param tsne_result A pre-calculated data.frame containing TSNE1 and TSNE2.
#' @param beta Beta matrix (only for colnames).
#' @param pd Phenotype data frame.
#' @param couleurs_origines Vector of colors for 'Origin'.
#' @param doPCA (Logical) Generate the PCA plot?
#' @param doTSNE (Logical) Generate the t-SNE plot?
#'
#' @return A list containing `pca_plot` (ggplot) and `tsne_plot` (ggplot).
#'
#' @noRd
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual guides guide_legend
#' @importFrom ggplot2 theme_minimal ggtitle labs facet_wrap theme
#' @importFrom ggrepel geom_text_repel
#'
mf_run_pca_tsne <- function(pca_result,
                            tsne_result = NULL, # <-- MODIFIED: Receives the result
                            beta,
                            pd,
                            couleurs_origines = NULL,
                            doPCA = TRUE,
                            doTSNE = TRUE) {

  p_champ <- NULL
  myLoad_t_SNE <- NULL

  # --- Common data preparation ---
  pd_data <- data.frame(
    Sample_Status = pd$Sample_Status,
    SampleName = colnames(beta),
    Tissue = pd$Tissue,
    Plate = pd$Sample_Plate,
    Origin = pd$Origin
  )

  # --- 1. PCA plot (Unchanged, was already correct) ---
  if (doPCA && !is.null(pca_result)) {
    pca_df <- cbind(
      as.data.frame(pca_result$x[, c("PC1", "PC2")]),
      pd_data
    )
    pc1_var <- round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2)
    pc2_var <- round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2)

    p_champ <- ggplot(
      pca_df,
      aes(x = PC1, y = PC2, color = Origin, shape = Sample_Status)
    ) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(
        aes(label = SampleName),
        size = 2.5, color = "black", box.padding = 0.4,
        point.padding = 0.3, max.overlaps = Inf, show.legend = FALSE
      ) +
      facet_wrap(~ Plate) +
      theme_minimal() +
      ggtitle("[Import/Filtering] PCA by Origin, Status (and Plate)") +
      labs(
        x = paste0("PC1 (", pc1_var, "%)"),
        y = paste0("PC2 (", pc2_var, "%)"),
        color = "Sample Origin",
        shape = "Sample Status"
      ) +
      theme(legend.position = "bottom")

    if (!is.null(couleurs_origines)) {
      p_champ <- p_champ +
        scale_color_manual(values = couleurs_origines, name = "Sample Origin")
    }
  }

  # --- 2. t-SNE plot (MODIFIED) ---
  # The Rtsne::Rtsne computation section has been REMOVED

  if (doTSNE && !is.null(tsne_result)) { # Checks if data was received

    message("Generating t-SNE plot from pre-calculated data.")

    # Combines received t-SNE results with metadata
    df_tsne_combined <- cbind(
      tsne_result, # Uses object provided as argument
      pd_data
    )

    myLoad_t_SNE <- ggplot(
      df_tsne_combined,
      aes(x = TSNE1, y = TSNE2, color = Origin, shape = Sample_Status)
    ) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(
        aes(label = SampleName),
        size = 2.5, color = "black", box.padding = 0.4,
        point.padding = 0.3, max.overlaps = Inf, show.legend = FALSE
      ) +
      facet_wrap(~ Plate) +
      theme_minimal() +
      ggtitle("[Import/Filtering] t-SNE by Origin, Status (and Plate)") +
      labs(
        x = "t-SNE 1", y = "t-SNE 2",
        color = "Sample Origin", shape = "Sample Status"
      ) +
      theme(legend.position = "bottom")

    if (!is.null(couleurs_origines)) {
      myLoad_t_SNE <- myLoad_t_SNE +
        scale_color_manual(values = couleurs_origines, name = "Sample Origin")
    }

  } else if (doTSNE) {
    message("WARNING: doTSNE=TRUE but tsne_result was NULL. t-SNE plot skipped.")
  }

  # --- Return ---
  return(list(
    pca_plot = p_champ,
    tsne_plot = myLoad_t_SNE
  ))
}

#' (Internal) Plot the SVD correlation heatmap
#'
#' @description
#' Uses PCAtools to generate the correlation heatmap from a
#' pre-calculated PCA object.
#'
#' @param p_obj The PRE-CALCULATED `pca` object from PCAtools.
#' @param pheno_columns Columns of `pd` to correlate.
#'
#' @return An `eigencorplot` graphic object.
#'
#' @noRd
#' @importFrom PCAtools eigencorplot
#'
mf_run_svd_cor <- function(p_obj, pheno_columns) { # <-- Signature MISE À JOUR

  # La section PCAtools::pca(...) a été SUPPRIMÉE

  if (is.null(p_obj)) {
    message("AVERTISSEMENT SVD: Objet PCA non fourni. Heatmap sautée.")
    return(NULL)
  }

  # 1. Identifier les métadonnées qui ont de la variance
  vars_qui_varient <- sapply(p_obj$metadata[, pheno_columns, drop = FALSE], function(col) {
    length(unique(col)) > 1
  })
  metavars_valides <- pheno_columns[vars_qui_varient]

  if(length(metavars_valides) == 0) {
    message("AVERTISSEMENT SVD: Aucune des pheno_columns fournies n'a de variance. Heatmap sautée.")
    return(NULL)
  }

  # 2. Identifier les PCs qui ont de la variance
  components_valides <- names(which(p_obj$variance > 0.001))
  if(length(components_valides) == 0) {
    message("AVERTISSEMENT SVD: Aucune composante PCA n'a de variance. Heatmap sautée.")
    return(NULL)
  }

  # 3. Générer la heatmap de corrélation
  svd_plot <- PCAtools::eigencorplot(
    p_obj,
    metavars = metavars_valides,
    components = components_valides,
    main = "SVD Correlation Heatmap (Phenotypes vs PCs)"
  )

  return(svd_plot)
}
