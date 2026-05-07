#' Check NGS project settings and data
#'
#' @description
#' Validates the full project configuration and phenotype data:
#' \enumerate{
#'   \item Configuration file structure
#'   \item Critical paths (project dir, pheno file, results dir)
#'   \item Phenotype file contents (required columns, IDs, treatments)
#'   \item Data files existence and sizes
#'   \item Import parameters (pipeline, assembly, context)
#'   \item QC parameters
#'   \item Analysis parameters (clustering, diffmeth, annotation)
#'   \item Statistical adjustments (batch, covariates)
#'   \item System resources
#' }
#'
#' @param params_file Path to \code{project_config.R}.
#'
#' @return Invisibly returns the \code{project_config} list if all
#'   checks pass.
#' @export
#' @importFrom readxl read_excel
mETHYLotest.NGS.Check_ProjectSettings <- function(params_file = NULL) {

  null_or <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b

  .fmt_bytes <- function(bytes) {
    if (is.na(bytes) || bytes == 0L) return("0 bytes")
    units <- c("bytes", "Kb", "Mb", "Gb")
    i <- min(floor(log(bytes, 1024)), length(units) - 1L)
    paste(sprintf("%.1f", bytes / (1024 ^ i)), units[i + 1L])
  }

  errors   <- character(0)
  warnings <- character(0)

  # ========================================================================
  # 1. LOAD CONFIGURATION
  # ========================================================================

  if (is.null(params_file) || !file.exists(params_file))
    stop("[mETHYLotest Check] Config not found: ", params_file)

  config_env <- new.env()
  tryCatch(
    source(params_file, local = config_env),
    error = function(e)
      stop("[mETHYLotest Check] Error sourcing config: ", e$message))

  if (!exists("project_config", envir = config_env))
    stop("[mETHYLotest Check] Config did not create 'project_config'.")

  cfg <- config_env$project_config

  message("[mETHYLotest Check] ========================================")
  message("[mETHYLotest Check] Project: ", cfg$project_name)
  message("[mETHYLotest Check] ========================================")

  # ========================================================================
  # 2. CRITICAL PATHS
  # ========================================================================

  message("[mETHYLotest Check] --- 1. Paths ---")

  if (!is.null(cfg$project_dir) && dir.exists(cfg$project_dir))
    message("[mETHYLotest Check] OK  Project dir: ", cfg$project_dir)
  else
    errors <- c(errors, paste("Project directory not found:",
                              cfg$project_dir))

  if (!is.null(cfg$pheno_file) && file.exists(cfg$pheno_file))
    message("[mETHYLotest Check] OK  Phenotype: ", cfg$pheno_file)
  else
    errors <- c(errors, paste("Phenotype file not found:",
                              cfg$pheno_file))

  for (d in c("res_dir", "rds_dir")) {
    if (!is.null(cfg[[d]])) {
      if (!dir.exists(cfg[[d]])) {
        dir.create(cfg[[d]], recursive = TRUE)
        message("[mETHYLotest Check] OK  ", d, " created: ", cfg[[d]])
      } else {
        message("[mETHYLotest Check] OK  ", d, ": ", cfg[[d]])
      }
    }
  }

  if (length(errors) > 0L)
    stop("[mETHYLotest Check] Path errors:\n",
         paste0("  - ", errors, collapse = "\n"))

  # ========================================================================
  # 3. PHENOTYPE DATA
  # ========================================================================

  message("[mETHYLotest Check] --- 2. Phenotype data ---")

  Pheno <- tryCatch({
    ext <- tools::file_ext(cfg$pheno_file)
    if (ext %in% c("xlsx", "xls"))
      as.data.frame(readxl::read_excel(cfg$pheno_file))
    else
      read.csv(cfg$pheno_file, stringsAsFactors = FALSE)
  }, error = function(e)
    stop("[mETHYLotest Check] Cannot read phenotype: ", e$message))

  message("[mETHYLotest Check] OK  Loaded: ", nrow(Pheno), " rows, ",
          ncol(Pheno), " columns.")

  col_id   <- if (!is.null(cfg$col_sampleID))  cfg$col_sampleID
  else "Sample_ID"
  col_path <- if (!is.null(cfg$col_file_path)) cfg$col_file_path
  else "Path"
  col_tx   <- if (!is.null(cfg$col_treatment)) cfg$col_treatment
  else "Treatment"

  required_cols <- c(col_id, col_path, col_tx)
  missing_cols  <- setdiff(required_cols, colnames(Pheno))
  if (length(missing_cols) > 0L)
    errors <- c(errors, paste("Missing column(s):",
                              paste(missing_cols, collapse = ", ")))
  else
    message("[mETHYLotest Check] OK  Required columns: ",
            paste(required_cols, collapse = ", "))

  extra_cols <- setdiff(colnames(Pheno), required_cols)
  if (length(extra_cols) > 0L)
    message("[mETHYLotest Check] INFO Extra columns: ",
            paste(extra_cols, collapse = ", "))

  if (length(errors) > 0L)
    stop("[mETHYLotest Check] Phenotype errors:\n",
         paste0("  - ", errors, collapse = "\n"))

  # ── Sample IDs ─────────────────────────────────────────────────────────
  ids <- Pheno[[col_id]]
  if (any(is.na(ids) | !nzchar(trimws(ids))))
    errors <- c(errors, "Sample IDs contain NA or empty values.")
  if (length(unique(ids)) != length(ids)) {
    dup <- unique(ids[duplicated(ids)])
    errors <- c(errors, paste("Duplicate IDs:",
                              paste(dup, collapse = ", ")))
  }
  if (any(grepl("\\s", ids)))
    errors <- c(errors, "Sample IDs contain spaces.")
  if (any(grepl("[^A-Za-z0-9_.-]", ids)))
    warnings <- c(warnings,
                  "Some IDs contain special characters.")

  if (!any(grepl("Sample_ID", paste(errors, collapse = " "))))
    message("[mETHYLotest Check] OK  ", length(unique(ids)),
            " unique sample IDs.")

  # ── Treatment ──────────────────────────────────────────────────────────
  tx <- as.numeric(as.character(Pheno[[col_tx]]))
  if (any(is.na(tx)) || !all(tx %in% c(0L, 1L))) {
    errors <- c(errors, "Treatment must be 0 (Control) or 1 (Case).")
  } else {
    n0 <- sum(tx == 0L); n1 <- sum(tx == 1L)
    message("[mETHYLotest Check] OK  Treatment: ",
            n0, " control, ", n1, " case.")
    if (n0 < 2L || n1 < 2L)
      warnings <- c(warnings,
                    "A group has < 2 samples. Power may be limited.")
  }

  # ── Data files ─────────────────────────────────────────────────────────
  message("[mETHYLotest Check] --- 3. Data files ---")

  paths <- Pheno[[col_path]]
  if (any(is.na(paths))) {
    errors <- c(errors, "File paths contain NA values.")
  } else {
    missing_f <- paths[!file.exists(paths)]
    if (length(missing_f) > 0L) {
      errors <- c(errors, paste0(
        "Missing ", length(missing_f), " file(s):\n",
        paste0("    - ", utils::head(missing_f, 10L),
               collapse = "\n"),
        if (length(missing_f) > 10L)
          paste0("\n    ... and ", length(missing_f) - 10L,
                 " more.")))
    } else {
      total_sz <- sum(file.info(paths)$size, na.rm = TRUE)
      message("[mETHYLotest Check] OK  All ", length(paths),
              " files found (total: ", .fmt_bytes(total_sz), ")")
      for (i in seq_along(paths))
        message("[mETHYLotest Check]     ", ids[i], ": ",
                .fmt_bytes(file.info(paths[i])$size))
    }
  }

  # Dimension consistency
  if (length(paths) != length(ids) || length(paths) != length(tx))
    errors <- c(errors, paste0(
      "Dimension mismatch: ", length(paths), " paths, ",
      length(ids), " IDs, ", length(tx), " treatments."))

  # ========================================================================
  # 4. IMPORT PARAMETERS
  # ========================================================================

  message("[mETHYLotest Check] --- 4. Import parameters ---")

  if (!is.null(cfg$assembly) && nzchar(cfg$assembly))
    message("[mETHYLotest Check] OK  Assembly: ", cfg$assembly)
  else
    errors <- c(errors, "Assembly not defined.")

  valid_ctx <- c("CpG", "CHG", "CHH")
  if (!is.null(cfg$context) && cfg$context %in% valid_ctx)
    message("[mETHYLotest Check] OK  Context: ", cfg$context)
  else
    errors <- c(errors, paste0("Invalid context: '", cfg$context, "'"))

  if (is.list(cfg$pipeline)) {
    needed <- c("chr.col", "start.col", "coverage.col", "freqC.col")
    miss_f <- setdiff(needed, names(cfg$pipeline))
    if (length(miss_f) > 0L)
      errors <- c(errors, paste("Custom pipeline missing:",
                                paste(miss_f, collapse = ", ")))
    else
      message("[mETHYLotest Check] OK  Pipeline: custom")
    if (!is.null(cfg$coord_offset) && cfg$coord_offset != 0L)
      message("[mETHYLotest Check] OK  Offset: ", cfg$coord_offset)
  } else if (is.character(cfg$pipeline)) {
    valid_p <- c("bismarkCytosineReport", "bismarkCoverage",
                 "bismark", "amp")
    if (cfg$pipeline %in% valid_p)
      message("[mETHYLotest Check] OK  Pipeline: ", cfg$pipeline)
    else
      errors <- c(errors, paste0("Invalid pipeline: '",
                                 cfg$pipeline, "'"))
  } else {
    errors <- c(errors, "Pipeline not defined.")
  }

  if (!is.null(cfg$min_coverage) && cfg$min_coverage >= 0)
    message("[mETHYLotest Check] OK  Min coverage: ", cfg$min_coverage)
  else
    warnings <- c(warnings, "min_coverage invalid. Default: 1.")

  valid_res <- c("base", "region")
  if (!is.null(cfg$resolution) && cfg$resolution %in% valid_res)
    message("[mETHYLotest Check] OK  Resolution: ", cfg$resolution)
  else
    warnings <- c(warnings, "Invalid resolution. Default: base.")

  # ========================================================================
  # 5. QC PARAMETERS
  # ========================================================================

  message("[mETHYLotest Check] --- 5. QC parameters ---")

  if (!is.null(cfg$qc_hi_perc)) {
    if (cfg$qc_hi_perc >= 90 && cfg$qc_hi_perc <= 100)
      message("[mETHYLotest Check] OK  hi.perc: ",
              cfg$qc_hi_perc, "%")
    else
      warnings <- c(warnings, paste0(
        "qc_hi_perc=", cfg$qc_hi_perc, " outside [90,100]."))
  }

  if (!is.null(cfg$qc_lo_count) && cfg$qc_lo_count >= 0)
    message("[mETHYLotest Check] OK  lo.count: ", cfg$qc_lo_count)

  # ========================================================================
  # 6. ANALYSIS PARAMETERS
  # ========================================================================

  message("[mETHYLotest Check] --- 6. Analysis parameters ---")

  valid_d <- c("correlation", "euclidean", "maximum", "manhattan")
  valid_m <- c("ward", "complete", "average", "single")
  if (!is.null(cfg$cluster_dist) && !cfg$cluster_dist %in% valid_d)
    warnings <- c(warnings, paste0("cluster_dist '",
                                   cfg$cluster_dist, "' unusual."))
  if (!is.null(cfg$cluster_method) && !cfg$cluster_method %in% valid_m)
    warnings <- c(warnings, paste0("cluster_method '",
                                   cfg$cluster_method, "' unusual."))

  valid_od <- c("MN", "shrinkMN", "none")
  valid_tt <- c("Chisq", "F", "midPval", "fast.fisher")
  if (!is.null(cfg$diff_overdispersion) &&
      !cfg$diff_overdispersion %in% valid_od)
    errors <- c(errors, paste0("Invalid diff_overdispersion: '",
                               cfg$diff_overdispersion, "'"))
  if (!is.null(cfg$diff_test) && !cfg$diff_test %in% valid_tt)
    errors <- c(errors, paste0("Invalid diff_test: '",
                               cfg$diff_test, "'"))
  if (!is.null(cfg$diff_cutoff) &&
      (cfg$diff_cutoff < 0 || cfg$diff_cutoff > 100))
    errors <- c(errors, "diff_cutoff must be [0,100].")
  if (!is.null(cfg$diff_qvalue) &&
      (cfg$diff_qvalue <= 0 || cfg$diff_qvalue > 1))
    errors <- c(errors, "diff_qvalue must be (0,1].")

  if (!is.null(cfg$annot_diff_cutoff) &&
      (cfg$annot_diff_cutoff < 0 || cfg$annot_diff_cutoff > 100))
    warnings <- c(warnings, "annot_diff_cutoff outside [0,100].")
  if (!is.null(cfg$annot_qval_cutoff) &&
      (cfg$annot_qval_cutoff <= 0 || cfg$annot_qval_cutoff > 1))
    warnings <- c(warnings, "annot_qval_cutoff outside (0,1].")

  # ========================================================================
  # 7. STATISTICAL ADJUSTMENTS
  # ========================================================================

  message("[mETHYLotest Check] --- 7. Adjustments ---")

  if (isTRUE(cfg$perform_batch_correction)) {
    if (is.null(cfg$batch_cols) || length(cfg$batch_cols) == 0L)
      errors <- c(errors, "Batch correction ON but batch_cols empty.")
    else {
      miss_b <- setdiff(cfg$batch_cols, colnames(Pheno))
      if (length(miss_b) > 0L)
        errors <- c(errors, paste("Batch col(s) missing:",
                                  paste(miss_b, collapse = ", ")))
      else
        message("[mETHYLotest Check] OK  Batch: ",
                paste(cfg$batch_cols, collapse = ", "))
    }
  }

  if (isTRUE(cfg$perform_covariate_adj)) {
    if (is.null(cfg$covariates) || length(cfg$covariates) == 0L)
      errors <- c(errors, "Covariates ON but list empty.")
    else {
      miss_c <- setdiff(cfg$covariates, colnames(Pheno))
      if (length(miss_c) > 0L)
        errors <- c(errors, paste("Covariate col(s) missing:",
                                  paste(miss_c, collapse = ", ")))
      else
        message("[mETHYLotest Check] OK  Covariates: ",
                paste(cfg$covariates, collapse = ", "))
    }
  }

  if (isTRUE(cfg$perform_batch_correction) &&
      isTRUE(cfg$perform_covariate_adj)) {
    ov <- intersect(cfg$batch_cols, cfg$covariates)
    if (length(ov) > 0L)
      warnings <- c(warnings, paste("Batch/covariate overlap:",
                                    paste(ov, collapse = ", ")))
  }

  # ========================================================================
  # 8. SYSTEM
  # ========================================================================

  message("[mETHYLotest Check] --- 8. System ---")

  if (!is.null(cfg$num_cores) && is.numeric(cfg$num_cores) &&
      cfg$num_cores >= 1L) {
    tot <- parallel::detectCores()
    if (cfg$num_cores > tot)
      warnings <- c(warnings, paste0("num_cores(", cfg$num_cores,
                                     ") > detected(", tot, ")."))
    message("[mETHYLotest Check] OK  Cores: ",
            cfg$num_cores, "/", tot)
  } else {
    warnings <- c(warnings, "Invalid num_cores. Default: 1.")
  }

  # ========================================================================
  # 9. SUMMARY
  # ========================================================================

  message("[mETHYLotest Check] ========================================")

  if (length(warnings) > 0L) {
    message("[mETHYLotest Check] ", length(warnings), " WARNING(s):")
    for (w in warnings)
      message("[mETHYLotest Check]   ! ", w)
  }

  if (length(errors) > 0L)
    stop("[mETHYLotest Check] ", length(errors), " ERROR(s):\n",
         paste0("  - ", errors, collapse = "\n"))

  steps <- c(
    "methRead",
    if (is.list(cfg$pipeline) && !is.null(cfg$coord_offset) &&
        cfg$coord_offset != 0L)
      paste0("offset(", cfg$coord_offset, ")"),
    "QC_loop",
    paste0("unite(destrand=",
           isTRUE(cfg$unite_destrand), ")"),
    paste0("clustering(",
           null_or(cfg$cluster_dist, "correlation"), "/",
           null_or(cfg$cluster_method, "ward"), ")"),
    paste0("diffMeth(",
           null_or(cfg$diff_overdispersion, "MN"), "/",
           null_or(cfg$diff_test, "Chisq"), ")"),
    if (isTRUE(cfg$perform_batch_correction))
      paste0("batch(", paste(cfg$batch_cols, collapse = "+"), ")"),
    if (isTRUE(cfg$perform_covariate_adj))
      paste0("cov(", paste(cfg$covariates, collapse = "+"), ")"),
    "annotation", "report"
  )

  message("[mETHYLotest Check] All checks passed.")
  message("[mETHYLotest Check] Samples:  ", nrow(Pheno))
  message("[mETHYLotest Check] Pipeline: ",
          paste(steps, collapse = " -> "))
  message("[mETHYLotest Check] ========================================")

  invisible(cfg)
}
