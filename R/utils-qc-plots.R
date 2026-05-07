#' (Interne) Trace la PCA et le t-SNE (Version simplifiée 2D)
#'
#' @description
#' Fonction utilitaire interne. Ne calcule rien, met en forme les
#' graphiques ggplot à partir des résultats pré-calculés.
#'
#' @param pca_result L'objet `prcomp` pré-calculé.
#' @param tsne_result Un data.frame pré-calculé contenant TSNE1 et TSNE2.
#' @param beta Matrice Beta (uniquement pour les colnames).
#' @param pd Data frame phénotypique.
#' @param couleurs_origines Vecteur de couleurs pour 'Origin'.
#' @param doPCA (Logique) Générer le graphique PCA ?
#' @param doTSNE (Logique) Générer le graphique t-SNE ?
#'
#' @return Une liste contenant `pca_plot` (ggplot) et `tsne_plot` (ggplot).
#'
#' @noRd
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual guides guide_legend
#' @importFrom ggplot2 theme_minimal ggtitle labs facet_wrap theme
#' @importFrom ggrepel geom_text_repel
#'
mf_run_pca_tsne <- function(pca_result,
                            tsne_result = NULL, # <-- MODIFIÉ : Reçoit le résultat
                            beta,
                            pd,
                            couleurs_origines = NULL,
                            doPCA = TRUE,
                            doTSNE = TRUE) {

  p_champ <- NULL
  myLoad_t_SNE <- NULL

  # --- Préparation des données communes ---
  pd_data <- data.frame(
    Sample_Status = pd$Sample_Status,
    SampleName = colnames(beta),
    Tissue = pd$Tissue,
    Plate = pd$Sample_Plate,
    Origin = pd$Origin
  )

  # --- 1. Graphique PCA (Inchangé, c'était déjà correct) ---
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
      ggtitle("[Import/Filtering] PCA par Origine, Statut (et Plaque)") +
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

  # --- 2. Graphique t-SNE (MODIFIÉ) ---
  # La section de calcul Rtsne::Rtsne a été SUPPRIMÉE

  if (doTSNE && !is.null(tsne_result)) { # Vérifie si on a reçu les données

    message("Génération du graphique t-SNE à partir des données pré-calculées.")

    # Combine les résultats t-SNE reçus avec les métadonnées
    df_tsne_combined <- cbind(
      tsne_result, # Utilise l'objet fourni en argument
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
      ggtitle("[Import/Filtering] t-SNE par Origine, Statut (et Plaque)") +
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
    message("AVERTISSEMENT: doTSNE=TRUE mais tsne_result était NULL. Graphique t-SNE sauté.")
  }

  # --- Retour ---
  return(list(
    pca_plot = p_champ,
    tsne_plot = myLoad_t_SNE
  ))
}

#' (Interne) Trace la heatmap de corrélation SVD
#'
#' @description
#' Utilise PCAtools pour générer la heatmap de corrélation à partir d'un
#' objet PCA pré-calculé.
#'
#' @param p_obj L'objet `pca` PRÉ-CALCULÉ de PCAtools.
#' @param pheno_columns Colonnes de `pd` à corréler.
#'
#' @return Un objet graphique `eigencorplot`.
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
