#' Generate HTML Summary Report
#'
#' @description
#' Compiles an R Markdown template to generate a comprehensive HTML report.
#'
#' @param project_dir Character. Absolute path to the mETHYLotest project directory.
#' @param output_file Character. Name of the output HTML file.
#' @param rmd_template Character. Path to the .Rmd template.
#' @param analyzed_samples Character vector. Sample IDs included in the analysis.
#' @param excluded_samples Character vector. Sample IDs excluded during QC.
#'
#' @return The path to the generated HTML report (invisibly).
#' @export
#' @importFrom rmarkdown render
mETHYLotest.NGS.GenerateReport <- function(project_dir,
                                           output_file      = "mETHYLotest_Analysis_Report.html",
                                           rmd_template     = NULL,
                                           analyzed_samples = NULL,
                                           excluded_samples = NULL) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required.")
  }

  if (is.null(rmd_template)) {
    rmd_template <- system.file("rmarkdown", "mETHYLotest.NGS.Report.Rmd",
                                package = "mETHYLotest")
    if (rmd_template == "") {
      candidates <- c(
        file.path(project_dir, "rmd", "mETHYLotest.NGS.Report.Rmd"),
        file.path(getwd(), "inst", "rmarkdown", "mETHYLotest.NGS.Report.Rmd"),
        file.path(getwd(), "mETHYLotest.NGS.Report.Rmd")
      )
      for (cc in candidates) {
        if (file.exists(cc)) { rmd_template <- cc; break }
      }
      if (rmd_template == "") {
        stop("Could not find RMD template.")
      }
    }
  }

  if (!file.exists(rmd_template)) {
    stop("Template not found: ", rmd_template)
  }

  output_path <- file.path(project_dir, output_file)
  message("[mETHYLotest] Generating report...")

  # ── Force local temp dir to avoid UNC/HTTP pandoc bug ──
  tmp_dir <- tempfile("report_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  tryCatch({
    rmarkdown::render(
      input             = rmd_template,
      output_file       = output_path,
      intermediates_dir = tmp_dir,
      encoding          = "UTF-8",
      params = list(
        project_dir      = project_dir,
        analyzed_samples = analyzed_samples,
        excluded_samples = excluded_samples
      ),
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
    message("[mETHYLotest] Report generated: ", output_path)
    return(invisible(output_path))
  }, error = function(e) {
    stop("[mETHYLotest] Report failed: ", e$message)
  })
}
