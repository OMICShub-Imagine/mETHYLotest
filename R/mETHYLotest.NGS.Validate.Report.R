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
