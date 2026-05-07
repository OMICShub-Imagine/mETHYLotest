#' Generate mETHYLotest EPIC Final Report
#'
#' @param project_dir Character. Project root.
#' @param beta Numeric matrix (probes × samples).
#' @param pheno Data.frame, named character vector, or unnamed character vector.
#'   If a vector, it is converted internally to a data.frame with columns
#'   \code{Sample_Name} and \code{Sample_Group}.
#' @param cfg List. The \code{project_config}.
#' @param dmp_list Named list of DMP data.frames (or NULL).
#' @param signatures_dir Character or NULL. Path to Episignatures/ folder.
#' @param validation_results Named list or NULL.
#' @param top_n Integer. Top DMPs for tables/heatmaps.
#' @param open Logical. Open in browser?
#'
#' @return Path to the HTML file (invisibly).
#' @export
mETHYLotest.EPIC.Report <- function(project_dir,
                                   beta,
                                   pheno,
                                   cfg,
                                   dmp_list           = NULL,
                                   signatures_dir     = NULL,
                                   validation_results = NULL,
                                   top_n              = 50L,
                                   open               = TRUE) {

  stopifnot(dir.exists(project_dir), is.matrix(beta), is.list(cfg))

  # ── Coerce pheno to data.frame ──
  if (is.data.frame(pheno)) {
    # Already fine
  } else if (is.character(pheno) || is.factor(pheno)) {
    if (!is.null(names(pheno))) {
      # Named vector: names = Sample_Name, values = Sample_Group
      pheno <- data.frame(
        Sample_Name  = names(pheno),
        Sample_Group = as.character(pheno),
        stringsAsFactors = FALSE
      )
    } else if (length(pheno) == ncol(beta)) {
      # Unnamed vector aligned to beta columns
      pheno <- data.frame(
        Sample_Name  = colnames(beta),
        Sample_Group = as.character(pheno),
        stringsAsFactors = FALSE
      )
    } else {
      stop("pheno is a vector but its length (", length(pheno),
           ") doesn't match ncol(beta) (", ncol(beta), ").")
    }
    message("[mETHYLotest Report] pheno coerced from vector to data.frame (",
            nrow(pheno), " samples).")
  } else {
    stop("pheno must be a data.frame, named character vector, or ",
         "character vector of length ncol(beta).")
  }

  for (pkg in c("rmarkdown", "ggplot2", "pheatmap"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' required.", pkg))

  if (is.data.frame(dmp_list))
    dmp_list <- list(Main = dmp_list)

  if (is.null(signatures_dir)) {
    sd <- file.path(cfg$res_dir, "Episignatures")
    if (dir.exists(sd)) signatures_dir <- sd
  }

  rmd_path <- system.file("rmarkdown", "mETHYLotest.EPIC.Report.Rmd",
                          package = "mETHYLotest")
  if (!nzchar(rmd_path))
    stop("Rmd template not found. Is mETHYLotest installed?")

  out_file <- file.path(normalizePath(project_dir, winslash = "/"),
                        "mETHYLotest_EPIC_Report.html")

  message("\n=== Rendering mETHYLotest EPIC Report ===")

  rmarkdown::render(
    input       = rmd_path,
    output_file = basename(out_file),
    output_dir  = dirname(out_file),
    params      = list(
      project_dir        = project_dir,
      beta               = beta,
      pheno              = pheno,
      cfg                = cfg,
      dmp_list           = dmp_list,
      signatures_dir     = signatures_dir,
      validation_results = validation_results,
      top_n              = as.integer(top_n)
    ),
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )

  message("Report saved: ", out_file)
  if (isTRUE(open) && interactive()) utils::browseURL(out_file)
  invisible(out_file)
}
