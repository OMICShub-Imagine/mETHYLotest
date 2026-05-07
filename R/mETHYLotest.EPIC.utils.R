# ==============================================================================
# EPIC_utils.R — QC plot helpers
# ==============================================================================


#' SVD Covariate Correlation Heatmap
#'
#' @description
#' Computes and plots the correlation between principal components (from a
#' \code{PCAtools::pca} object) and phenotype variables. Uses
#' \code{PCAtools::eigencorplot} internally.
#'
#' @param p_obj         A \code{PCAtools::pca} object (output of
#'   \code{PCAtools::pca()}).
#' @param pheno_columns Character vector of column names from \code{p_obj$metadata}
#'   to correlate with PCs. Only numeric-convertible or factor columns are kept.
#' @param n_pcs         Integer. Maximum number of PCs to display (default 10).
#' @param corFUN        Character. Correlation method passed to
#'   \code{eigencorplot}: \code{"pearson"}, \code{"spearman"} or
#'   \code{"kendall"} (default \code{"pearson"}).
#'
#' @return A \code{ggplot} / \code{recordedplot} object, or \code{NULL} if the
#'   plot cannot be generated.
#'
#' @export
#' @importFrom PCAtools eigencorplot
mETHYLotest.EPIC.utils.SVDCorrelation <- function(p_obj,
                                                 pheno_columns,
                                                 n_pcs   = 10L,
                                                 corFUN  = "pearson") {

  if (is.null(p_obj))
    stop("'p_obj' must be a PCAtools::pca object.")

  if (is.null(p_obj$metadata) || nrow(p_obj$metadata) == 0L)
    stop("'p_obj$metadata' is empty. ",
         "Re-run PCAtools::pca() with a non-NULL 'metadata' argument.")

  # ── 1. Filter pheno_columns to those present in metadata ─────────────────
  available <- intersect(pheno_columns, colnames(p_obj$metadata))
  if (length(available) == 0L) {
    warning("[mETHYLotest SVD] None of the requested pheno_columns (",
            paste(pheno_columns, collapse = ", "),
            ") found in pca metadata. Available: ",
            paste(colnames(p_obj$metadata), collapse = ", "))
    return(NULL)
  }

  dropped <- setdiff(pheno_columns, available)
  if (length(dropped) > 0L)
    message("[mETHYLotest SVD] Columns not found in metadata (skipped): ",
            paste(dropped, collapse = ", "))

  # ── 2. Keep only columns that can be coerced to numeric ───────────────────
  # eigencorplot requires numeric or factor metadata
  usable <- Filter(function(col) {
    x <- p_obj$metadata[[col]]
    is.numeric(x) || is.factor(x) || is.logical(x) ||
      !all(is.na(suppressWarnings(as.numeric(as.character(x)))))
  }, available)

  # After filtering usable columns — before running eigencorplot
  if (length(usable) < 2L) {
    warning(
      "[mETHYLotest SVD] At least 2 numeric/factor variables are required for ",
      "eigencorplot. Only ", length(usable), " found after filtering: ",
      if (length(usable) == 1L) paste(usable, collapse = ", ") else "none",
      ". SVD heatmap skipped. Add biological covariates (e.g. Batch, Sex, Age) ",
      "to your metadata file."
    )
    return(NULL)
  }

  if (length(usable) == 0L) {
    warning("[mETHYLotest SVD] No numeric/factor columns remain after filtering.")
    return(NULL)
  }

  if (length(usable) < length(available))
    message("[mETHYLotest SVD] Non-numeric columns dropped: ",
            paste(setdiff(available, usable), collapse = ", "))

  # ── 3. Coerce → numeric, detect and handle large IDs ──────────────────────
  # Sentrix IDs (12-digit numbers like 207127940107) and other identifier
  # columns break eigencorplot internal matrix arithmetic.
  # Strategy: values > 1e9 OR all-unique character → convert to factor ranks.

  meta_clean <- as.data.frame(
    lapply(usable, function(col_name) {
      x       <- p_obj$metadata[[col_name]]
      n_total <- length(x[!is.na(x)])
      n_uniq  <- length(unique(x[!is.na(x)]))

      # Detect all-unique columns (likely IDs, not biological covariates)
      if (n_uniq == n_total && n_total > 2L) {
        message("[mETHYLotest SVD] '", col_name,
                "' has all-unique values (likely an ID column) — ",
                "converting to factor ranks.")
        return(as.numeric(as.factor(x)))
      }

      if (is.numeric(x) || is.logical(x)) {
        # Large numbers (Sentrix IDs etc.) → factor ranks
        if (is.numeric(x) && max(abs(x), na.rm = TRUE) > 1e9) {
          message("[mETHYLotest SVD] '", col_name,
                  "' contains large numbers (>1e9) — ",
                  "converting to factor ranks.")
          return(as.numeric(as.factor(x)))
        }
        return(as.numeric(x))
      }

      # Character / factor
      num <- suppressWarnings(as.numeric(as.character(x)))
      if (!all(is.na(num))) {
        if (max(abs(num), na.rm = TRUE) > 1e9) {
          message("[mETHYLotest SVD] '", col_name,
                  "' parses as large numbers — converting to factor ranks.")
          return(as.numeric(as.factor(x)))
        }
        return(num)
      }

      as.numeric(as.factor(x))
    }),
    stringsAsFactors = FALSE
  )
  colnames(meta_clean) <- usable
  rownames(meta_clean) <- rownames(p_obj$metadata)

  p_clean          <- p_obj
  p_clean$metadata <- meta_clean

  # ── 4. Limit to n_pcs ─────────────────────────────────────────────────────
  n_pcs_avail <- ncol(p_obj$rotated)
  n_pcs_use   <- min(n_pcs, n_pcs_avail)

  message("[mETHYLotest SVD] Running eigencorplot: ",
          n_pcs_use, " PCs x ",
          length(usable), " variable(s): ",
          paste(usable, collapse = ", "))

  # ── 5. Generate plot ───────────────────────────────────────────────────────
  tryCatch({
    plt <- PCAtools::eigencorplot(
      pcaobj        = p_clean,
      components    = PCAtools::getComponents(p_clean,
                                              seq_len(n_pcs_use)),
      metavars      = usable,
      col           = c("darkblue", "blue2", "white", "red2", "darkred"),
      cexCorval     = 0.7,
      fontCorval    = 2L,
      posLab        = "bottomleft",
      rotLabX       = 45,
      posColKey     = "top",
      cexLabColKey  = 1.2,
      scale         = TRUE,
      corFUN        = corFUN,
      corUSE        = "pairwise.complete.obs",
      corMultipleTestCorrection = "BH",
      signifSymbols = c("***", "**", "*", ""),
      signifCutpoints = c(0, 0.001, 0.01, 0.05, 1)
    )
    message("[mETHYLotest SVD] Eigencorplot generated.")
    plt
  }, error = function(e) {
    warning("[mETHYLotest SVD] eigencorplot failed: ", e$message)
    NULL
  })
}


#' Static PCA and t-SNE Plots for HTML Report
#'
#' @description
#' Generates static \code{ggplot2} scatter plots of PCA and t-SNE results
#' for inclusion in the HTML QC report.
#'
#' @param pca_result    Output of \code{stats::prcomp}. \code{NULL} skips PCA.
#' @param tsne_result   Data frame with columns \code{TSNE1}, \code{TSNE2},
#'   row-named by sample. \code{NULL} skips t-SNE.
#' @param beta          Beta matrix (probes x samples) — used to retrieve
#'   sample names.
#' @param pd            Phenotype data frame, row-named by sample.
#' @param sample_colors Named character vector of colours per group level.
#'   \code{NULL} uses ggplot2 defaults.
#' @param color_col     Character. Column of \code{pd} used for point colour
#'   (default \code{"Treatment"}, then \code{"Sample_Status"}).
#' @param doPCA         Logical. Generate PCA plot? (default \code{TRUE})
#' @param doTSNE        Logical. Generate t-SNE plot? (default \code{TRUE})
#'
#' @return A named list: \code{$pca_plot} and \code{$tsne_plot}. Each element
#'   is a \code{ggplot} object or \code{NULL}.
#'
#' @export
#' @import ggplot2
mETHYLotest.EPIC.utils.PlotPCATSNE <- function(pca_result    = NULL,
                                              tsne_result   = NULL,
                                              beta,
                                              pd,
                                              sample_colors = NULL,
                                              color_col     = NULL,
                                              doPCA         = TRUE,
                                              doTSNE        = TRUE) {

  out <- list(pca_plot = NULL, tsne_plot = NULL)

  # ── Resolve colour column ─────────────────────────────────────────────────
  if (is.null(color_col)) {
    color_col <- if ("Treatment"    %in% colnames(pd)) "Treatment"
    else if ("Sample_Status" %in% colnames(pd)) "Sample_Status"
    else colnames(pd)[1L]
  }
  if (!color_col %in% colnames(pd)) {
    warning("[mETHYLotest Plot] '", color_col,
            "' not found in pd. Colour disabled.")
    color_col <- NULL
  }

  # ── Shared ggplot2 theme ──────────────────────────────────────────────────
  base_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title  = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # ── PCA plot ──────────────────────────────────────────────────────────────
  if (doPCA && !is.null(pca_result)) {
    message("[mETHYLotest Plot] Generating static PCA plot...")

    var_exp <- summary(pca_result)$importance[2L, ] * 100
    df      <- as.data.frame(pca_result$x[, 1:2, drop = FALSE])
    df$Sample <- rownames(df)

    if (!is.null(color_col))
      df[[color_col]] <- pd[df$Sample, color_col]

    p <- ggplot2::ggplot(
      df,
      if (!is.null(color_col))
        ggplot2::aes(x = PC1, y = PC2,
                     colour = .data[[color_col]],
                     label  = Sample)
      else
        ggplot2::aes(x = PC1, y = PC2, label = Sample)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::geom_text(
        nudge_y  = diff(range(df$PC2)) * 0.03,
        size     = 3,
        check_overlap = TRUE
      ) +
      ggplot2::labs(
        title  = "PCA — PC1 vs PC2",
        x      = sprintf("PC1 (%.1f%%)", var_exp[["PC1"]]),
        y      = sprintf("PC2 (%.1f%%)", var_exp[["PC2"]]),
        colour = color_col
      ) +
      base_theme

    if (!is.null(sample_colors) && !is.null(color_col))
      p <- p + ggplot2::scale_colour_manual(values = sample_colors)

    out$pca_plot <- p
  }

  # ── t-SNE plot ────────────────────────────────────────────────────────────
  if (doTSNE && !is.null(tsne_result)) {
    message("[mETHYLotest Plot] Generating static t-SNE plot...")

    df        <- tsne_result
    df$Sample <- rownames(df)

    if (!is.null(color_col))
      df[[color_col]] <- pd[df$Sample, color_col]

    p <- ggplot2::ggplot(
      df,
      if (!is.null(color_col))
        ggplot2::aes(x = TSNE1, y = TSNE2,
                     colour = .data[[color_col]],
                     label  = Sample)
      else
        ggplot2::aes(x = TSNE1, y = TSNE2, label = Sample)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::geom_text(
        nudge_y  = diff(range(df$TSNE2)) * 0.03,
        size     = 3,
        check_overlap = TRUE
      ) +
      ggplot2::labs(
        title  = "t-SNE",
        x      = "t-SNE 1",
        y      = "t-SNE 2",
        colour = color_col
      ) +
      base_theme

    if (!is.null(sample_colors) && !is.null(color_col))
      p <- p + ggplot2::scale_colour_manual(values = sample_colors)

    out$tsne_plot <- p
  }

  out
}

# ==============================================================================
# EPIC QC Plot Utilities
# ==============================================================================


# ── Beta Value Density ─────────────────────────────────────────────────────────

#' Beta Value Distribution Plot
#'
#' @description
#' Per-sample density plot of beta values using native plotly traces.
#' The number of probes used for density estimation adapts automatically
#' to the array size and sample count:
#' \itemize{
#'   \item If \code{n_probes_sample = NULL} (default): uses all probes
#'     when n_samples <= 24, otherwise subsamples adaptively.
#'   \item Explicit integer value: forces that subsample size.
#'   \item \code{Inf}: always uses all probes.
#' }
#'
#' @param myLoad            ChAMP object with \code{$beta} and \code{$pd}.
#' @param sample_colors     Named colour vector per group level.
#' @param n_probes_sample   Integer, \code{NULL} (adaptive), or \code{Inf}
#'   (all probes). See description.
#' @param seed              Integer. Random seed for reproducible subsampling
#'   (default 42).
#'
#' @return A \code{plotly} object or \code{NULL}.
#' @export
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom stats density
mETHYLotest.EPIC.utils.PlotBetaDensity <- function(myLoad,
                                                  sample_colors   = NULL,
                                                  n_probes_sample = NULL,
                                                  seed            = 42L) {
  if (is.null(myLoad$beta)) {
    warning("[mETHYLotest Beta Density] '$beta' not found in myLoad.")
    return(NULL)
  }

  message("[mETHYLotest Beta Density] Generating beta distribution plot...")

  beta     <- myLoad$beta
  n_probes <- nrow(beta)
  n_samps  <- ncol(beta)

  # ── Adaptive subsampling logic ────────────────────────────────────────────
  #
  # Target: enough probes for stable density, few enough for fast rendering.
  #
  # Rules:
  #   <= 24 samples  → all probes (density() is fast, shapes matter)
  #   25-48 samples  → 200 000 probes max
  #   > 48 samples   → 100 000 probes max
  #
  # Override with explicit n_probes_sample or Inf for all probes.

  target <- if (!is.null(n_probes_sample)) {
    if (is.infinite(n_probes_sample)) n_probes
    else as.integer(n_probes_sample)
  } else if (n_samps <= 24L) {
    n_probes   # all probes
  } else if (n_samps <= 48L) {
    min(200000L, n_probes)
  } else {
    min(100000L, n_probes)
  }

  sampled <- target < n_probes
  if (sampled) {
    set.seed(seed)
    beta <- beta[sample(n_probes, target), ]
    message("[mETHYLotest Beta Density] Subsampled: ",
            format(target, big.mark = ","), " / ",
            format(n_probes, big.mark = ","), " probes (",
            round(target / n_probes * 100, 1), "%).")
  } else {
    message("[mETHYLotest Beta Density] Using all ",
            format(n_probes, big.mark = ","), " probes.")
  }

  # ── Group colour mapping ──────────────────────────────────────────────────
  pd_cols   <- colnames(myLoad$pd)
  color_col <- if ("Treatment"     %in% pd_cols) "Treatment"
  else if ("Sample_Status" %in% pd_cols) "Sample_Status"
  else pd_cols[1L]

  groups <- setNames(
    as.character(myLoad$pd[[color_col]]),
    rownames(myLoad$pd)
  )

  # ── Build one trace per sample ────────────────────────────────────────────
  p <- plotly::plot_ly()

  for (samp in colnames(beta)) {
    vals <- beta[, samp]
    vals <- vals[!is.na(vals) & vals >= 0 & vals <= 1]
    if (length(vals) < 2L) next

    dens <- stats::density(vals, from = 0, to = 1, n = 512L)
    grp  <- groups[[samp]]

    col <- if (!is.null(sample_colors) && grp %in% names(sample_colors))
      sample_colors[[grp]]
    else NULL

    p <- plotly::add_trace(
      p,
      x             = dens$x,
      y             = dens$y,
      type          = "scatter",
      mode          = "lines",
      name          = samp,
      legendgroup   = grp,
      line          = if (!is.null(col)) list(color = col, width = 1.5)
      else list(width = 1.5),
      text          = paste0("Sample: ", samp, "<br>Group: ", grp),
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        "Beta: %{x:.3f}<br>",
        "Density: %{y:.4f}",
        "<extra></extra>"
      )
    )
  }

  # ── Subtitle / caption ────────────────────────────────────────────────────
  caption <- if (sampled)
    paste0("Density estimated on ",
           format(target, big.mark = ","), " / ",
           format(n_probes, big.mark = ","),
           " probes (", round(target / n_probes * 100, 1), "% — random subsample, seed = ", seed, ")")
  else
    paste0("Density estimated on all ",
           format(n_probes, big.mark = ","), " probes")

  p %>%
    plotly::layout(
      title     = list(
        text = paste0("Beta Value Distribution per Sample",
                      "<br><sup>", caption, "</sup>"),
        x    = 0.5
      ),
      xaxis     = list(
        title    = "Beta Value",
        range    = c(0, 1),
        tickvals = seq(0, 1, 0.25)
      ),
      yaxis     = list(title = "Density"),
      hovermode = "x unified",
      legend    = list(title = list(text = color_col))
    )
}


# ── Sample Correlation Heatmap ─────────────────────────────────────────────────

#' Sample Correlation Heatmap (interactive)
#'
#' @description
#' Computes the Pearson (or Spearman) inter-sample correlation matrix and
#' renders it as an interactive \code{plotly} heatmap with hierarchical
#' clustering order. Outlier samples are flagged via an attribute.
#'
#' Uses \code{plotly::plot_ly(type = "heatmap")} directly to avoid
#' \code{ggplotly} bridging issues with \code{geom_tile} + \code{geom_text}.
#'
#' @param myLoad  ChAMP object with \code{$beta} and \code{$pd}.
#' @param method  \code{"pearson"} (default) or \code{"spearman"}.
#'
#' @return A \code{plotly} object with an \code{"outlier_samples"} attribute,
#'   or \code{NULL}.
#' @export
#' @importFrom plotly plot_ly layout
#' @importFrom stats hclust as.dist
mETHYLotest.EPIC.utils.PlotCorHeatmap <- function(myLoad,
                                                 method = "pearson") {
  if (is.null(myLoad$beta)) {
    warning("[mETHYLotest Cor Heatmap] '$beta' not found in myLoad.")
    return(NULL)
  }

  message("[mETHYLotest Cor Heatmap] Computing sample correlation matrix (",
          method, ")...")

  beta    <- myLoad$beta
  pd_cols <- colnames(myLoad$pd)

  # ── Correlation matrix ────────────────────────────────────────────────────
  cor_mat <- cor(beta, use = "pairwise.complete.obs", method = method)

  # ── Hierarchical clustering order ─────────────────────────────────────────
  hc           <- hclust(as.dist(1 - cor_mat), method = "complete")
  sample_order <- hc$labels[hc$order]
  cor_ordered  <- cor_mat[sample_order, sample_order]

  # ── Outlier detection ─────────────────────────────────────────────────────
  mean_cor    <- colMeans(cor_mat)
  outlier_thr <- mean(mean_cor) - 2 * stats::sd(mean_cor)
  outliers    <- names(mean_cor)[mean_cor < outlier_thr]

  if (length(outliers) > 0L)
    message("[mETHYLotest Cor Heatmap] Potential outlier(s): ",
            paste(outliers, collapse = ", "))
  else
    message("[mETHYLotest Cor Heatmap] No outliers detected.")

  # ── Colour scale ──────────────────────────────────────────────────────────
  cor_min <- max(0.7, floor(min(cor_mat) * 10) / 10)

  colorscale <- list(
    c(0,   "#2166ac"),
    c(0.5, "white"),
    c(1,   "#d6604d")
  )

  # ── Annotation bar (group colour strip) ───────────────────────────────────
  color_col <- if ("Treatment"     %in% pd_cols) "Treatment"
  else if ("Sample_Status" %in% pd_cols) "Sample_Status"
  else NULL

  # ── Hover text matrix ─────────────────────────────────────────────────────
  hover_text <- outer(
    sample_order, sample_order,
    FUN = function(x, y)
      paste0("Sample X: ", x,
             "<br>Sample Y: ", y,
             "<br>Correlation: ",
             round(cor_ordered[cbind(x, y)], 4))
  )

  # ── Build plotly heatmap ──────────────────────────────────────────────────
  # Direct plot_ly — avoids all ggplotly geom_tile / geom_text issues
  p <- plotly::plot_ly(
    z             = cor_ordered,
    x             = sample_order,
    y             = sample_order,
    type          = "heatmap",
    colorscale    = colorscale,
    zmin          = cor_min,
    zmax          = 1,
    text          = hover_text,
    hovertemplate = "%{text}<extra></extra>",
    colorbar      = list(
      title  = paste0(tools::toTitleCase(method), "<br>Correlation"),
      len    = 0.6
    )
  )

  # ── Add correlation value annotations ─────────────────────────────────────
  # Only for small datasets (annotations become unreadable at large n)
  n_samp <- ncol(beta)

  if (n_samp <= 20L) {
    annotations <- list()
    for (i in seq_along(sample_order)) {
      for (j in seq_along(sample_order)) {
        val <- cor_ordered[sample_order[i], sample_order[j]]
        annotations <- c(annotations, list(
          list(
            x         = sample_order[j],
            y         = sample_order[i],
            text      = sprintf("%.2f", val),
            showarrow = FALSE,
            font      = list(
              size  = if (n_samp <= 10L) 11 else 8,
              color = if (abs(val - 0.85) < 0.1) "black" else "black"
            )
          )
        ))
      }
    }
    p <- p %>% plotly::layout(annotations = annotations)
  }

  # ── Layout ────────────────────────────────────────────────────────────────
  subtitle <- if (length(outliers) > 0L)
    paste0("Potential outlier(s): ", paste(outliers, collapse = ", "))
  else
    "No outliers detected (criterion: mean cor &lt; mean - 2 SD)"

  p <- p %>%
    plotly::layout(
      title  = list(
        text = paste0(
          "Sample Correlation Heatmap (",
          tools::toTitleCase(method), ")",
          "<br><sup>", subtitle, "</sup>"
        ),
        x = 0.5
      ),
      xaxis  = list(
        title     = NULL,
        tickangle = 45,
        side      = "bottom"
      ),
      yaxis  = list(
        title     = NULL,
        autorange = "reversed"
      ),
      margin = list(l = 100, b = 100, t = 80, r = 80)
    )

  # Attach outliers as attribute for the UI warning box
  attr(p, "outlier_samples") <- outliers
  p
}


# ── Detection P-value Summary ──────────────────────────────────────────────────

#' Detection P-value Summary Plot
#'
#' @description
#' Computes the percentage of probes failing the detection p-value threshold
#' per sample and displays the result as a bar chart. Samples exceeding the
#' recommended 5 percent removal threshold are highlighted in red.
#'
#' @param myLoad     ChAMP object with \code{$detP} and \code{$pd}.
#' @param det_p_cut  Numeric. Detection p-value threshold (default 0.01).
#'
#' @return A \code{ggplot2} object or \code{NULL}.
#' @export
#' @import ggplot2
mETHYLotest.EPIC.utils.PlotDetPSummary <- function(myLoad,
                                                  det_p_cut = 0.01) {
  if (is.null(myLoad$detP)) {
    warning("[mETHYLotest DetP] '$detP' not found in myLoad. ",
            "Run champ.load() with method = 'ChAMP'.")
    return(NULL)
  }

  message("[mETHYLotest DetP] Computing detection p-value summary...")

  pd_cols   <- colnames(myLoad$pd)
  color_col <- if ("Treatment"     %in% pd_cols) "Treatment"
  else if ("Sample_Status" %in% pd_cols) "Sample_Status"
  else pd_cols[1L]

  n_probes   <- nrow(myLoad$detP)
  pct_failed <- colMeans(myLoad$detP > det_p_cut, na.rm = TRUE) * 100

  fail_threshold <- 5   # standard EPIC QC threshold (%)

  df_detp <- data.frame(
    Sample     = names(pct_failed),
    Pct_Failed = as.numeric(pct_failed),
    stringsAsFactors = FALSE
  )
  df_detp[[color_col]] <- myLoad$pd[df_detp$Sample, color_col]
  df_detp$Flag         <- df_detp$Pct_Failed > fail_threshold

  # Order from highest to lowest failure rate
  df_detp$Sample <- factor(
    df_detp$Sample,
    levels = df_detp$Sample[order(df_detp$Pct_Failed, decreasing = TRUE)]
  )

  n_flagged <- sum(df_detp$Flag)

  p <- ggplot2::ggplot(
    df_detp,
    ggplot2::aes(x = Sample, y = Pct_Failed, fill = Flag)
  ) +
    ggplot2::geom_bar(stat = "identity", colour = "black",
                      linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = fail_threshold,
                        colour    = "red",
                        linetype  = "dashed",
                        linewidth = 0.8) +
    ggplot2::annotate(
      "text",
      x      = Inf,
      y      = fail_threshold + max(pct_failed) * 0.03,
      label  = paste0("Removal threshold: ", fail_threshold, "%"),
      hjust  = 1.05, size = 3.5, colour = "red"
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
      labels = c("FALSE" = "Pass", "TRUE" = "Fail"),
      name   = paste0("detP > ", det_p_cut)
    ) +
    ggplot2::labs(
      title    = "Detection P-value QC per Sample",
      subtitle = if (n_flagged > 0)
        paste0(n_flagged, " sample(s) exceed the ",
               fail_threshold, "% failure threshold")
      else
        "All samples pass the detection p-value threshold",
      x        = "Sample",
      y        = paste0("% Probes with detP > ", det_p_cut),
      caption  = paste0(format(n_probes, big.mark = ","),
                        " probes evaluated")
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(
        hjust  = 0.5, size = 11,
        colour = if (n_flagged > 0) "red" else "grey40"
      ),
      plot.caption     = ggplot2::element_text(size = 9, colour = "grey50"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title     = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


# ── Sex Prediction ─────────────────────────────────────────────────────────────

#' Predict Sex from Chromosomal Methylation
#'
#' @description
#' Computes median beta values on chrX and chrY probes per sample and
#' classifies samples as Male (chrY median > 0.1) or Female. If a declared
#' sex column (\code{Sex} or \code{Gender}) is present in \code{$pd},
#' mismatches between declared and predicted sex are flagged.
#'
#' @param myLoad  ChAMP object with \code{$beta} and \code{$pd}.
#'
#' @return A named list:
#'   \itemize{
#'     \item \code{$plot}       ggplot2 object.
#'     \item \code{$data}       Data frame with per-sample predictions.
#'     \item \code{$mismatches} Character vector of mismatched sample names.
#'   }
#'   Returns \code{NULL} if annotation cannot be obtained.
#'
#' @export
#' @import ggplot2
mETHYLotest.EPIC.utils.PlotSexPrediction <- function(myLoad) {

  if (is.null(myLoad$beta)) {
    warning("[mETHYLotest Sex] '$beta' not found in myLoad.")
    return(NULL)
  }

  message("[mETHYLotest Sex] Predicting sex from chromosomal methylation...")

  beta    <- myLoad$beta
  pd_cols <- colnames(myLoad$pd)

  # Infer array type from Sample_Plate
  arraytype <- "EPICv1"
  if ("Sample_Plate" %in% pd_cols) {
    plates <- unique(myLoad$pd$Sample_Plate)
    if (any(grepl("EPICv2", plates, ignore.case = TRUE)))
      arraytype <- "EPICv2"
    else if (any(grepl("450",  plates, ignore.case = TRUE)))
      arraytype <- "450K"
  }
  message("[mETHYLotest Sex] Array type for annotation: ", arraytype)

  # Get annotation
  anno <- .mf_get_sex_chr_annotation(rownames(beta), arraytype)
  if (is.null(anno)) {
    warning("[mETHYLotest Sex] Chromosome annotation unavailable. ",
            "Sex prediction skipped.")
    return(NULL)
  }

  probes_x <- intersect(anno$probe[anno$chr %in% c("chrX", "X")],
                        rownames(beta))
  probes_y <- intersect(anno$probe[anno$chr %in% c("chrY", "Y")],
                        rownames(beta))

  message("[mETHYLotest Sex] chrX probes: ", length(probes_x),
          " | chrY probes: ", length(probes_y))

  if (length(probes_x) == 0 && length(probes_y) == 0) {
    warning("[mETHYLotest Sex] No sex-chromosome probes found.")
    return(NULL)
  }

  # Median beta per sample per sex chromosome
  df_sex <- data.frame(
    Sample   = colnames(beta),
    median_X = if (length(probes_x) > 0)
      apply(beta[probes_x, , drop = FALSE], 2,
            stats::median, na.rm = TRUE)
    else rep(NA_real_, ncol(beta)),
    median_Y = if (length(probes_y) > 0)
      apply(beta[probes_y, , drop = FALSE], 2,
            stats::median, na.rm = TRUE)
    else rep(NA_real_, ncol(beta)),
    stringsAsFactors = FALSE
  )

  # Prediction rule: chrY median > 0.1 = Male
  df_sex$Predicted_Sex <- ifelse(
    !is.na(df_sex$median_Y) & df_sex$median_Y > 0.1,
    "Male", "Female"
  )

  # Compare with declared sex if available
  sex_col <- c("Sex", "Gender", "sex", "gender", "SEX", "GENDER")
  sex_col <- sex_col[sex_col %in% pd_cols]
  sex_col <- if (length(sex_col) > 0L) sex_col[1L] else NULL

  if (!is.null(sex_col)) {
    df_sex$Declared_Sex <- myLoad$pd[df_sex$Sample, sex_col]
    declared_norm       <- tolower(trimws(df_sex$Declared_Sex))
    is_pred_male        <- df_sex$Predicted_Sex == "Male"
    is_decl_male        <- grepl("^m", declared_norm)
    is_decl_female      <- grepl("^f", declared_norm)
    df_sex$Match        <- ifelse(
      (is_pred_male  & is_decl_male) | (!is_pred_male & is_decl_female),
      "Match", "Mismatch"
    )
  }

  mismatches <- if ("Match" %in% colnames(df_sex))
    df_sex$Sample[df_sex$Match == "Mismatch"]
  else character(0)

  # Build plot depending on available sex chromosomes
  has_x <- length(probes_x) > 0L
  has_y <- length(probes_y) > 0L

  p <- if (has_x && has_y)
    .mf_sex_scatter(df_sex, probes_x, probes_y, sex_col)
  else if (has_x)
    .mf_sex_bar(df_sex, "median_X", probes_x, "chrX", "No chrY probes")
  else
    .mf_sex_bar(df_sex, "median_Y", probes_y, "chrY", "No chrX probes")

  if (length(mismatches) > 0)
    message("[mETHYLotest Sex] WARNING: ",
            length(mismatches), " mismatch(es): ",
            paste(mismatches, collapse = ", "))
  else
    message("[mETHYLotest Sex] No declared/predicted sex mismatches.")

  list(plot = p, data = df_sex, mismatches = mismatches)
}


# ── Internal helpers ────────────────────────────────────────────────────────────

#' @keywords internal
.mf_get_sex_chr_annotation <- function(probe_ids, arraytype) {

  # ── 1. Sélection du package d'annotation minfi ────────────────────────────
  pkg <- if (grepl("EPICv2", arraytype, ignore.case = TRUE))
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
  else if (grepl("EPIC", arraytype, ignore.case = TRUE))
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  else if (grepl("450", arraytype, ignore.case = TRUE))
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  else {
    message("[mETHYLotest Sex] Unknown array type '", arraytype,
            "' — trying EPICv1 annotation as fallback.")
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  }

  # ── 2. Tentative via Locations dataset (minfi annotation packages) ─────────
  # These packages expose data as named datasets: Locations, Manifest, Other…
  # NOT as a single object named after the package.
  if (requireNamespace(pkg, quietly = TRUE)) {

    local_env <- new.env(parent = emptyenv())
    anno <- tryCatch({

      utils::data("Locations", package = pkg, envir = local_env)
      locs <- local_env[["Locations"]]

      if (is.null(locs))
        stop("'Locations' dataset loaded as NULL.")

      if (!"chr" %in% colnames(locs))
        stop("'chr' column not found in Locations. ",
             "Available: ", paste(colnames(locs), collapse = ", "))

      result <- data.frame(
        probe = rownames(locs),
        chr   = as.character(locs[["chr"]]),
        stringsAsFactors = FALSE
      )
      result <- result[!is.na(result$probe) & !is.na(result$chr), ]

      n_x <- sum(result$chr %in% c("chrX", "X"))
      n_y <- sum(result$chr %in% c("chrY", "Y"))
      message("[mETHYLotest Sex] Annotation from ", pkg,
              "::Locations | probes: ", nrow(result),
              " | chrX: ", n_x, " | chrY: ", n_y)
      result

    }, error = function(e) {
      message("[mETHYLotest Sex] ", pkg, "::Locations failed: ", e$message)
      NULL
    })

    if (!is.null(anno)) return(anno)

  } else {
    message("[mETHYLotest Sex] '", pkg, "' not installed.",
            "\n  Install: BiocManager::install('", pkg, "')")
  }

  # ── 3. Fallback ChAMP probe.features (EPICv1 / 450K uniquement) ───────────
  # Not available for EPICv2
  champ_dataset <- if (grepl("EPICv1|^EPIC$", arraytype, ignore.case = TRUE))
    "probe.features.epic"
  else if (grepl("450", arraytype, ignore.case = TRUE))
    "probe.features"
  else
    NULL

  if (!is.null(champ_dataset) && requireNamespace("ChAMP", quietly = TRUE)) {

    local_env <- new.env(parent = emptyenv())
    anno <- tryCatch({
      utils::data(list    = champ_dataset,
                  package = "ChAMP",
                  envir   = local_env)
      feat    <- local_env[[ls(local_env)[1L]]]
      chr_col <- intersect(c("CHR", "chr"), colnames(feat))[1L]

      if (is.na(chr_col))
        stop("No chromosome column. Available: ",
             paste(colnames(feat), collapse = ", "))

      result <- data.frame(
        probe = rownames(feat),
        chr   = as.character(feat[[chr_col]]),
        stringsAsFactors = FALSE
      )
      result <- result[!is.na(result$probe) & !is.na(result$chr), ]
      n_x    <- sum(result$chr %in% c("chrX", "X"))
      n_y    <- sum(result$chr %in% c("chrY", "Y"))
      message("[mETHYLotest Sex] Annotation from ChAMP::", champ_dataset,
              " | probes: ", nrow(result),
              " | chrX: ", n_x, " | chrY: ", n_y)
      result

    }, error = function(e) {
      message("[mETHYLotest Sex] ChAMP::", champ_dataset,
              " failed: ", e$message)
      NULL
    })

    if (!is.null(anno)) return(anno)
  }

  # ── 4. Toutes les sources ont échoué ──────────────────────────────────────
  message("[mETHYLotest Sex] No annotation source available for '",
          arraytype, "'.",
          "\n  EPICv2 : BiocManager::install('",
          "IlluminaHumanMethylationEPICv2anno.20a1.hg38')",
          "\n  EPICv1 : BiocManager::install('",
          "IlluminaHumanMethylationEPICanno.ilm10b4.hg19')",
          "\n  450K   : BiocManager::install('",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19')")
  NULL
}


#' @keywords internal
.mf_sex_scatter <- function(df_sex, probes_x, probes_y, sex_col) {

  has_match <- "Match" %in% colnames(df_sex)

  base_aes <- if (has_match)
    ggplot2::aes(x      = median_X,
                 y      = median_Y,
                 colour = Predicted_Sex,
                 shape  = Match,
                 label  = Sample)
  else
    ggplot2::aes(x      = median_X,
                 y      = median_Y,
                 colour = Predicted_Sex,
                 label  = Sample)

  y_range <- diff(range(df_sex$median_Y, na.rm = TRUE))

  p <- ggplot2::ggplot(df_sex, base_aes) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggplot2::geom_text(
      nudge_y       = max(y_range * 0.03, 0.005),
      size          = 3,
      check_overlap = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed",
                        colour = "grey50", linewidth = 0.6) +
    ggplot2::annotate("text",
                      x = -Inf, y = 0.105,
                      label  = "Y threshold (0.1)",
                      hjust  = -0.05, size = 3, colour = "grey40") +
    ggplot2::scale_colour_manual(
      values = c("Female" = "#e74c3c", "Male" = "#3498db"),
      name   = "Predicted Sex"
    ) +
    ggplot2::labs(
      title   = "Predicted Sex from Methylation Data",
      x       = paste0("Median Beta - chrX (", length(probes_x), " probes)"),
      y       = paste0("Median Beta - chrY (", length(probes_y), " probes)"),
      caption = "Males: chrY median > 0.1  |  Females: chrY median <= 0.1"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption     = ggplot2::element_text(size = 9, colour = "grey50"),
      legend.title     = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (has_match)
    p <- p + ggplot2::scale_shape_manual(
      values = c("Match" = 16L, "Mismatch" = 4L),
      name   = "vs Declared"
    )

  # Circle mismatches in red
  if (has_match) {
    mm <- df_sex[df_sex$Match == "Mismatch", ]
    if (nrow(mm) > 0)
      p <- p + ggplot2::geom_point(data   = mm,
                                   colour = "red",
                                   size   = 7,
                                   shape  = 1,
                                   stroke = 1.5)
  }

  p
}


#' @keywords internal
.mf_sex_bar <- function(df_sex, y_col, probes, chr_label, caption_note) {
  df_sex$y_val <- df_sex[[y_col]]
  ggplot2::ggplot(df_sex,
                  ggplot2::aes(x = Sample, y = y_val,
                               fill = Predicted_Sex)) +
    ggplot2::geom_bar(stat = "identity", colour = "black",
                      linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      values = c("Female" = "#e74c3c", "Male" = "#3498db"),
      name   = "Predicted Sex"
    ) +
    ggplot2::labs(
      title   = paste0("Predicted Sex - ", chr_label, " Methylation"),
      x       = "Sample",
      y       = paste0("Median Beta - ", chr_label,
                       " (", length(probes), " probes)"),
      caption = caption_note
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, colour = "grey50"),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
    )
}


#' Cell Type Correction via ChAMP (champ.refbase)
#'
#' @description
#' Wrapper around \code{ChAMP::champ.refbase()}. Estimates cell type
#' proportions using the reference-based Houseman method and returns
#' the corrected beta matrix.
#'
#' @param beta       Numeric matrix (probes x samples), normalised.
#' @param pd         Phenotype data frame, row-named by sample.
#'   Must contain \code{Sample_Name}.
#' @param arraytype  \code{"EPIC"}, \code{"EPICv2"}, or \code{"450K"}.
#' @param save_plots Logical. Save proportion bar plot? (default TRUE)
#' @param output_dir Directory for plots and proportion table.
#'
#' @return Named list: \code{$beta_corrected}, \code{$proportions}.
#'   Returns \code{NULL} on failure.
#'
#' @export
#' @importFrom ChAMP champ.refbase
mETHYLotest.EPIC.utils.CellTypeCorrection <- function(beta,
                                                     pd,
                                                     arraytype  = "EPIC",
                                                     save_plots = TRUE,
                                                     output_dir = NULL) {

  if (!requireNamespace("ChAMP", quietly = TRUE))
    stop("[mETHYLotest CellType] Package 'ChAMP' is required.")

  message("[mETHYLotest CellType] Running champ.refbase() ",
          "(arraytype = ", arraytype, ")...")

  # ── 1. Run champ.refbase ──────────────────────────────────────────────────
  refbase_result <- tryCatch({

    ChAMP::champ.refbase(
      beta      = beta,
      pd        = pd,
      arraytype = arraytype
    )

  }, error = function(e) {
    stop("[mETHYLotest CellType] champ.refbase() failed: ", e$message)
  })

  beta_corrected <- refbase_result$beta
  proportions    <- as.data.frame(refbase_result$CellFraction)
  proportions$Sample <- rownames(proportions)

  cell_cols <- setdiff(colnames(proportions), "Sample")
  message("[mETHYLotest CellType] Cell types estimated: ",
          paste(cell_cols, collapse = ", "))

  # ── 2. Save outputs ───────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Proportion table
    prop_csv <- file.path(output_dir, "cell_type_proportions.csv")
    write.csv(proportions, prop_csv, row.names = FALSE)
    message("[mETHYLotest CellType] Proportions saved: ", prop_csv)

    # Bar plot
    if (save_plots && requireNamespace("ggplot2", quietly = TRUE)) {

      prop_long <- reshape(
        proportions,
        varying   = cell_cols,
        v.names   = "Proportion",
        timevar   = "CellType",
        times     = cell_cols,
        direction = "long"
      )

      p <- ggplot2::ggplot(
        prop_long,
        ggplot2::aes(x = Sample, y = Proportion, fill = CellType)
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_y_continuous(
          labels = scales::percent_format(accuracy = 1)
        ) +
        ggplot2::labs(
          title = "Estimated Cell Type Proportions (champ.refbase)",
          x     = NULL,
          y     = "Proportion"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
          legend.title = ggplot2::element_text(face = "bold")
        )

      plot_png <- file.path(output_dir, "cell_type_proportions.png")
      ggplot2::ggsave(plot_png, p, width = 10, height = 5, dpi = 150)
      message("[mETHYLotest CellType] Plot saved: ", plot_png)
    }
  }

  message("[mETHYLotest CellType] Correction complete.")

  list(
    beta_corrected = beta_corrected,
    proportions    = proportions
  )
}
