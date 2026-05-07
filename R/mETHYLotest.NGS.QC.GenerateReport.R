#' Generate a static HTML QC report
#'
#' @param df_meth Dataframe of methylation by chromosome
#' @param df_controls Dataframe of control genome stats (current iteration)
#' @param df_controls_all Dataframe of control genome stats for ALL samples (including excluded)
#' @param df_stats Dataframe of general position stats
#' @param df_global_meth Dataframe of global methylation per sample
#' @param min_cov The coverage threshold used
#' @param excluded_samples Vector of excluded sample IDs
#' @param project_dir Project directory (for loading saved RDS files in the report)
#' @param output_file Path to the output HTML file
#'
#' @export
mETHYLotest.NGS.GenerateQCReport <- function(df_meth,
                                             df_controls,
                                             df_controls_all = NULL,
                                             df_stats,
                                             df_global_meth = NULL,
                                             min_cov,
                                             excluded_samples,
                                             project_dir = NULL,
                                             output_file) {

  template_path <- system.file("rmarkdown", "mETHYLotest.NGS.QC_Report.Rmd",
                               package = "mETHYLotest")

  if (template_path == "") {
    stop("Template Rmd not found. Check inst/rmarkdown/mETHYLotest.NGS.QC_Report.Rmd")
  }

  message("Generating QC report: ", output_file)

  tmp_dir <- tempfile("qc_report_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  rmarkdown::render(
    input             = template_path,
    output_file       = output_file,
    intermediates_dir = tmp_dir,
    encoding          = "UTF-8",
    params = list(
      df_meth          = df_meth,
      df_controls      = df_controls,
      df_controls_all  = df_controls_all,
      df_stats         = df_stats,
      df_global_meth   = df_global_meth,
      min_cov          = min_cov,
      excluded_samples = excluded_samples,
      project_dir      = project_dir
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  message("QC report generated successfully.")
}
