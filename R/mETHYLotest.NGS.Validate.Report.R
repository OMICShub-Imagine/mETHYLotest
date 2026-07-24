#' Generate HTML Validation Report for WGBS Long-Read
#'
#' @description
#' Generates a comprehensive HTML report from \code{mETHYLotest.NGS.validate()}
#' results. Mirrors the Shiny Validator UI visualisations (PCA, Heatmap,
#' SVM metrics, Variable Importance).
#'
#' @param validation_results List returned by \code{mETHYLotest.NGS.validate()}.
#' @param output_dir Directory to save the report.
#' @param filename Name of the output HTML file.
#' @param rmd_path Optional path to Rmd template override.
#'
#' @return Path to the generated HTML file (invisibly).
#' @export
mETHYLotest.NGS.Validate.Report <- function(validation_results,
                                               output_dir = getwd(),
                                               filename   = "NGS_Validation_Report.html",
                                               rmd_path   = NULL) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("[mETHYLotest] Package 'rmarkdown' is required.")
  }

  if (is.null(validation_results) || length(validation_results) == 0L) {
    stop("[mETHYLotest Validate Report] Validation results are empty.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # ── Export JSON Metrics and PNGs for Native UI ──
  tryCatch({
    message("[mETHYLotest Validate Report] Exporting native JSON metrics and PNGs...")
    json_data <- list()
    for (sig in names(validation_results)) {
      res <- validation_results[[sig]]
      safe_sig <- gsub("[^[:alnum:]_]", "_", sig)

      pca_df <- NULL
      if (!is.null(res$pca)) {
        pca_df <- as.data.frame(res$pca$x[, 1:2])
        pca_df$Class <- as.character(res$raw_data$Class)
        pca_df$Sample <- rownames(pca_df)
        
        # Save PCA PNG
        p <- ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1, y=PC2, color=Class)) +
             ggplot2::geom_point(size=3) + ggplot2::theme_minimal() +
             ggplot2::labs(title=paste("PCA -", sig))
        ggplot2::ggsave(file.path(output_dir, paste0("Validation_PCA_", safe_sig, ".png")), p, width=6, height=5, bg="white")
      }

      var_imp <- NULL
      if (!is.null(res$svm_model)) {
        imp <- tryCatch(caret::varImp(res$svm_model, scale=FALSE)$importance, error=function(e) NULL)
        if (!is.null(imp)) {
          score <- if ("Overall" %in% colnames(imp)) imp$Overall else apply(imp, 1, max)
          var_imp <- data.frame(Locus = rownames(imp), Score = score)
          var_imp <- head(var_imp[order(var_imp$Score, decreasing=TRUE), ], 20)
        }
      }

      cm <- NULL
      if (!is.null(res$confusion)) {
        cm <- as.data.frame(res$confusion)
      }

      metrics <- list(
        accuracy = res$accuracy,
        silhouette = res$silhouette,
        auc = res$auc,
        sensitivity = res$sensitivity,
        specificity = res$specificity
      )

      json_data[[sig]] <- list(
        pca = pca_df,
        varImp = var_imp,
        confusion = cm,
        metrics = metrics,
        foldAccuracies = res$fold_accuracies
      )
    }
    jsonlite::write_json(json_data, file.path(output_dir, "Validation_Metrics.json"), auto_unbox = TRUE)
  }, error = function(e) {
    warning("[mETHYLotest Validate Report] Failed to export JSON/PNG: ", e$message)
  })

  # ── Template resolution ──
  if (is.null(rmd_path)) {
    rmd_template <- system.file("rmarkdown",
                                "mETHYLotest.NGS.validate_report.Rmd",
                                package = "mETHYLotest")
    if (!nzchar(rmd_template)) {
      # Fallback locations
      candidates <- c(
        file.path(getwd(), "inst", "rmarkdown", "mETHYLotest.NGS.validate_report.Rmd"),
        file.path(getwd(), "mETHYLotest.NGS.validate_report.Rmd")
      )
      for (cc in candidates) {
        if (file.exists(cc)) { rmd_template <- cc; break }
      }
      if (!nzchar(rmd_template)) {
        stop("[mETHYLotest Validate Report] Rmd template not found. Provide 'rmd_path'.")
      }
    }
  } else {
    rmd_template <- rmd_path
  }

  if (!file.exists(rmd_template)) {
    stop("[mETHYLotest Validate Report] Template not found: ", rmd_template)
  }

  message("[mETHYLotest Validate Report] Template  : ", rmd_template)
  message("[mETHYLotest Validate Report] Signatures: ", length(validation_results))

  # ── Temp data ──
  temp_rds <- file.path(output_dir, ".temp_val_results.rds")
  saveRDS(validation_results, temp_rds)

  output_file <- file.path(output_dir, filename)

  message("[mETHYLotest Validate Report] Rendering...")

  tryCatch({
    rmarkdown::render(
      input       = rmd_template,
      output_file = filename,
      output_dir  = output_dir,
      params      = list(data_path = temp_rds),
      envir       = new.env(),
      quiet       = TRUE
    )
    message("[mETHYLotest Validate Report] Generated: ", output_file)
  }, error = function(e) {
    warning("[mETHYLotest Validate Report] Rendering failed: ", e$message)
  }, finally = {
    if (file.exists(temp_rds)) file.remove(temp_rds)
  })

  invisible(output_file)
}
