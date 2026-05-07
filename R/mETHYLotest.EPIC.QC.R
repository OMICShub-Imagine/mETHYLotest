#' EPIC Array QC Pipeline
#'
#' @description
#' Runs the full QC pipeline on a ChAMP \code{myLoad} object:
#' \enumerate{
#'   \item PCA (\code{stats::prcomp})
#'   \item t-SNE (\code{Rtsne})
#'   \item SVD covariate correlation (\code{PCAtools::pca})
#'   \item Beta value density distribution
#'   \item Sample correlation heatmap
#'   \item Detection p-value summary
#'   \item Sex prediction from chromosomal methylation
#'   \item Excel export of all QC metrics
#'   \item HTML report generation from Rmd template
#'   \item Interactive sample-removal UI (\code{mETHYLotest.EPIC.QC.UI})
#' }
#'
#' @param myLoad             ChAMP object with \code{$beta} and \code{$pd}.
#' @param beta_matrix        Optional replacement beta matrix — overrides
#'   \code{myLoad$beta} if provided.
#' @param outputDir          Output directory for Excel and HTML report.
#'   Defaults to \code{tempdir()/mETHYLotest_QC} if \code{NULL}.
#' @param pheno_columns      Columns of \code{$pd} used for SVD correlation.
#'   Auto-detected from extra columns if \code{NULL}.
#' @param sample_colors      Named colour vector per group label
#'   (e.g. \code{c("CTL" = "steelblue", "CASE" = "firebrick")}).
#' @param doPCA              Logical. Compute PCA? (default \code{TRUE})
#' @param doTSNE             Logical. Compute t-SNE? (default \code{TRUE})
#' @param doSVD              Logical. Compute SVD heatmap? (default \code{TRUE})
#' @param doBetaDensity      Logical. Plot beta distribution? (default \code{TRUE})
#' @param doCorHeatmap       Logical. Plot sample correlation heatmap?
#'   (default \code{TRUE})
#' @param doDetPSummary      Logical. Plot detection p-value summary?
#'   Requires \code{myLoad$detP}. (default \code{TRUE})
#' @param det_p_cut          Numeric. Detection p-value threshold used for
#'   the summary plot. (default \code{0.01})
#' @param doSexPrediction    Logical. Predict sex from chromosomal methylation?
#'   (default \code{TRUE})
#' @param run_interactive_ui Logical. Launch interactive UI? (default \code{TRUE})
#'
#' @return Invisibly returns a character vector of sample names flagged for
#'   removal. Returns \code{character(0)} if no samples are flagged or the UI
#'   is skipped.
#'
#' @seealso
#'   \code{\link{mETHYLotest.EPIC.QC.UI}},
#'   \code{\link{mETHYLotest.EPIC.utils.SVDCorrelation}},
#'   \code{\link{mETHYLotest.EPIC.utils.PlotPCATSNE}},
#'   \code{\link{mETHYLotest.EPIC.utils.PlotBetaDensity}},
#'   \code{\link{mETHYLotest.EPIC.utils.PlotCorHeatmap}},
#'   \code{\link{mETHYLotest.EPIC.utils.PlotDetPSummary}},
#'   \code{\link{mETHYLotest.EPIC.utils.PlotSexPrediction}}
#'
#' @export
#' @importFrom stats prcomp
#' @importFrom Rtsne Rtsne
#' @importFrom PCAtools pca
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom rmarkdown render
mETHYLotest.EPIC.QC <- function(myLoad,
                               beta_matrix        = NULL,
                               outputDir          = NULL,
                               pheno_columns      = NULL,
                               sample_colors      = NULL,
                               doPCA              = TRUE,
                               doTSNE             = TRUE,
                               doSVD              = TRUE,
                               doBetaDensity      = TRUE,   # fixed: was missing
                               doCorHeatmap       = TRUE,   # fixed: was missing
                               doDetPSummary      = TRUE,   # fixed: was missing
                               det_p_cut          = 0.01,   # fixed: was hardcoded
                               doSexPrediction    = TRUE,   # fixed: was missing
                               run_interactive_ui = TRUE) {

  message("[mETHYLotest EPIC QC] Starting QC pipeline.")

  # ── 1. Validation ──────────────────────────────────────────────────────────
  if (is.null(myLoad$beta) || is.null(myLoad$pd))
    stop("'myLoad' must contain '$beta' and '$pd'.")

  if (!is.null(beta_matrix)) {
    message("[mETHYLotest EPIC QC] Replacing beta matrix from argument.")
    myLoad$beta <- beta_matrix
  }

  if (!all(colnames(myLoad$beta) %in% rownames(myLoad$pd)))
    stop(
      "Column names of '$beta' do not all match row names of '$pd'. ",
      "Ensure rownames(myLoad$pd) are set to sample names."
    )

  # ── 2. Output directory ────────────────────────────────────────────────────
  if (is.null(outputDir))
    outputDir <- file.path(tempdir(), "mETHYLotest_QC")

  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
    message("[mETHYLotest EPIC QC] Output directory created: ", outputDir)
  }

  # ── 3. Auto-detect pheno_columns for SVD ──────────────────────────────────
  if (is.null(pheno_columns)) {
    # Slide and Array are Sentrix technical identifiers, not biological
    # covariates — exclude them from SVD correlation analysis
    reserved <- c("Sample_Name", "Sample_ID", "Basename",
                  "Sample_Plate", "Treatment", "Sample_Status",
                  "Slide", "Array", "Sentrix_ID", "Sentrix_Position")
    pheno_columns <- setdiff(colnames(myLoad$pd), reserved)

    if (length(pheno_columns) == 0L) {
      pheno_columns <- colnames(myLoad$pd)[1L]
      message("[mETHYLotest EPIC QC] No extra pd columns for SVD; ",
              "falling back to '", pheno_columns, "'.")
    } else {
      message("[mETHYLotest EPIC QC] SVD columns auto-detected: ",
              paste(pheno_columns, collapse = ", "))
    }
  }

  # ── 4. Initialisation ─────────────────────────────────────────────────────
  samples_to_remove <- character(0)   # always defined

  qc_results <- list(
    doPCA           = doPCA,
    doTSNE          = doTSNE,
    doSVD           = doSVD,
    doBetaDensity   = doBetaDensity,
    doCorHeatmap    = doCorHeatmap,
    doDetPSummary   = doDetPSummary,
    doSexPrediction = doSexPrediction,
    plots           = list(),
    pca_prcomp      = NULL,
    pca_tools       = NULL
  )
  excel_data_list <- list()
  tsne_df         <- NULL

  # ── 5. PCA (stats::prcomp) ─────────────────────────────────────────────────
  if (doPCA || doTSNE || run_interactive_ui) {
    message("[mETHYLotest EPIC QC] Computing PCA (stats::prcomp)...")

    qc_results$pca_prcomp <- stats::prcomp(
      t(myLoad$beta), center = TRUE, scale. = TRUE
    )

    if (doPCA)
      excel_data_list$PCA_Scores <-
      as.data.frame(qc_results$pca_prcomp$x)

    # ── 6. t-SNE ──────────────────────────────────────────────────────────
    if (doTSNE || run_interactive_ui) {
      message("[mETHYLotest EPIC QC] Computing t-SNE (Rtsne)...")

      n_pcs     <- min(30L, ncol(qc_results$pca_prcomp$x))
      pcs_input <- qc_results$pca_prcomp$x[, seq_len(n_pcs), drop = FALSE]
      n_samp    <- nrow(pcs_input)
      max_perp  <- floor((n_samp - 1) / 3)

      if (max_perp < 1L) {
        message("[mETHYLotest EPIC QC] t-SNE skipped: too few samples (n = ",
                n_samp, ").")
      } else {
        tsne_res          <- Rtsne::Rtsne(pcs_input,
                                          perplexity       = min(30L, max_perp),
                                          check_duplicates = FALSE,
                                          pca              = FALSE)
        tsne_df           <- as.data.frame(tsne_res$Y)
        colnames(tsne_df) <- c("TSNE1", "TSNE2")
        rownames(tsne_df) <- colnames(myLoad$beta)

        if (doTSNE)
          excel_data_list$tSNE_Coordinates <- tsne_df

        message("[mETHYLotest EPIC QC] t-SNE complete (perplexity = ",
                min(30L, max_perp), ").")
      }
    }
  }

  # ── 7. SVD (PCAtools::pca) ─────────────────────────────────────────────────
  if (doSVD || run_interactive_ui) {
    message("[mETHYLotest EPIC QC] Computing SVD (PCAtools::pca)...")

    pd_aligned <- myLoad$pd[colnames(myLoad$beta), , drop = FALSE]

    qc_results$pca_tools <- PCAtools::pca(
      myLoad$beta,
      metadata  = pd_aligned,
      center    = TRUE,
      removeVar = 0.1
    )

    if (doSVD) {
      message("[mETHYLotest EPIC QC] Generating SVD correlation heatmap...")
      qc_results$plots$svd_plot <- mETHYLotest.EPIC.utils.SVDCorrelation(
        p_obj         = qc_results$pca_tools,
        pheno_columns = pheno_columns
      )
      # eigencorplot returns a recordedplot/gtable, not a ggplot
      # → no $data slot: skip Excel export for SVD
    }
  }

  # ── 8. Static PCA / t-SNE plots (for HTML report) ─────────────────────────
  if (doPCA || doTSNE) {
    message("[mETHYLotest EPIC QC] Generating static PCA / t-SNE plots...")
    plot_list <- mETHYLotest.EPIC.utils.PlotPCATSNE(
      pca_result    = qc_results$pca_prcomp,
      tsne_result   = tsne_df,
      beta          = myLoad$beta,
      pd            = myLoad$pd,
      sample_colors = sample_colors,
      doPCA         = doPCA,
      doTSNE        = doTSNE
    )
    qc_results$plots$pca_plot  <- plot_list$pca_plot
    qc_results$plots$tsne_plot <- plot_list$tsne_plot
  }

  # ── 9. Beta density ────────────────────────────────────────────────────────
  # fixed: was after HTML report → now before so it appears in the report
  if (doBetaDensity) {
    message("[mETHYLotest EPIC QC] Generating beta density plot...")
    qc_results$plots$beta_density <-
      mETHYLotest.EPIC.utils.PlotBetaDensity(myLoad, sample_colors)
  }

  # ── 10. Correlation heatmap ────────────────────────────────────────────────
  if (doCorHeatmap) {
    message("[mETHYLotest EPIC QC] Generating sample correlation heatmap...")
    qc_results$plots$cor_heatmap <-
      mETHYLotest.EPIC.utils.PlotCorHeatmap(myLoad)

    # Export correlation matrix to Excel
    if (!is.null(qc_results$plots$cor_heatmap)) {
      cor_mat <- cor(myLoad$beta,
                     use    = "pairwise.complete.obs",
                     method = "pearson")
      excel_data_list$Sample_Correlation <- as.data.frame(cor_mat)
    }
  }

  # ── 11. Detection p-value summary ─────────────────────────────────────────
  if (doDetPSummary) {
    if (is.null(myLoad$detP)) {
      message("[mETHYLotest EPIC QC] '$detP' not found — ",
              "detection p-value summary skipped.")
    } else {
      message("[mETHYLotest EPIC QC] Generating detection p-value summary...")
      qc_results$plots$detp_summary <-
        mETHYLotest.EPIC.utils.PlotDetPSummary(myLoad,
                                              det_p_cut = det_p_cut)

      # Export % failed per sample to Excel
      pct_failed <- colMeans(myLoad$detP > det_p_cut,
                             na.rm = TRUE) * 100
      excel_data_list$DetP_Summary <- data.frame(
        Sample         = names(pct_failed),
        Pct_Failed     = round(as.numeric(pct_failed), 4L),
        Flag           = pct_failed > 5,
        stringsAsFactors = FALSE
      )
    }
  }

  # ── 12. Sex prediction ─────────────────────────────────────────────────────
  if (doSexPrediction) {
    message("[mETHYLotest EPIC QC] Predicting sex from methylation data...")
    qc_results$plots$sex_pred <-
      mETHYLotest.EPIC.utils.PlotSexPrediction(myLoad)

    # Export sex prediction table to Excel
    if (!is.null(qc_results$plots$sex_pred$data))
      excel_data_list$Sex_Prediction <- qc_results$plots$sex_pred$data
  }

  # ── 13. Excel export ───────────────────────────────────────────────────────
  # fixed: moved after all plots so all data is captured
  if (length(excel_data_list) > 0L) {
    message("[mETHYLotest EPIC QC] Saving QC metrics to Excel...")

    xl_path <- file.path(normalizePath(outputDir), "QC_metrics.xlsx")
    wb      <- openxlsx::createWorkbook()

    for (sheet in names(excel_data_list)) {
      openxlsx::addWorksheet(wb, sheet)
      openxlsx::writeData(wb,
                          sheet    = sheet,
                          x        = excel_data_list[[sheet]],
                          rowNames = TRUE)
    }

    openxlsx::saveWorkbook(wb, file = xl_path, overwrite = TRUE)
    message("[mETHYLotest EPIC QC] Excel saved: ", xl_path)
  }

  # ── 14. HTML report ────────────────────────────────────────────────────────
  # fixed: moved after all plots so report contains everything
  tmpl <- system.file("rmarkdown", "mETHYLotest.EPIC.QC_Report.Rmd",
                      package = "mETHYLotest")

  if (nchar(tmpl) == 0L) {
    warning("[mETHYLotest EPIC QC] Rmd template not found; ",
            "HTML report skipped.")
  } else {
    out_html <- file.path(normalizePath(outputDir),
                          "mETHYLotest_EPIC_QC_Report.html")
    tryCatch({
      rmarkdown::render(tmpl,
                        output_file = out_html,
                        params      = list(qc_data = qc_results),
                        quiet       = TRUE)
      message("[mETHYLotest EPIC QC] HTML report saved: ", out_html)
    }, error = function(e)
      warning("[mETHYLotest EPIC QC] HTML report failed: ", e$message))
  }

  # ── 15. Interactive UI ─────────────────────────────────────────────────────
  if (run_interactive_ui) {
    message("[mETHYLotest EPIC QC] Launching interactive QC UI...")

    # ── 11. Interactive UI ─────────────────────────────────────────────────────
    samples_to_remove <- mETHYLotest.EPIC.QC.UI(
      myLoad       = myLoad,
      pca_result   = qc_results$pca_prcomp,
      tsne_result  = tsne_df,
      svd_plot     = qc_results$plots$svd_plot,
      beta_density = qc_results$plots$beta_density,
      cor_heatmap  = qc_results$plots$cor_heatmap,
      detp_summary = qc_results$plots$detp_summary,
      sex_pred     = qc_results$plots$sex_pred
    )

    if (is.null(samples_to_remove))
      samples_to_remove <- character(0)

    if (length(samples_to_remove) > 0L)
      message("[mETHYLotest EPIC QC] ", length(samples_to_remove),
              " sample(s) flagged: ",
              paste(samples_to_remove, collapse = ", "))
    else
      message("[mETHYLotest EPIC QC] No samples flagged for removal.")
  }

  message("[mETHYLotest EPIC QC] Pipeline complete.")
  invisible(samples_to_remove)
}
