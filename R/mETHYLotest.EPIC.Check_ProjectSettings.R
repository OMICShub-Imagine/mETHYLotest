#' Check project settings from a configuration file
#'
#' @description
#' Reads a \code{project_config.R} file and validates all settings before
#' running the pipeline:
#' \enumerate{
#'   \item Required paths exist (\code{idat_dir}, \code{pheno_file}).
#'   \item Phenotype file can be read and contains required columns
#'         (\code{Sample_Name}, \code{Sample_Plate}, \code{Sample_Group},
#'          \code{Sentrix_ID}, \code{Sentrix_Position}).
#'   \item \code{Sample_Plate} contains only valid array types.
#'   \item \code{Sentrix_Position} matches expected format.
#'   \item \code{Sample_Name} entries are unique.
#'   \item All corresponding \code{_Red.idat} and \code{_Grn.idat} files
#'         exist (derived from \code{Sentrix_ID_Sentrix_Position}).
#'   \item Configuration parameters are consistent.
#' }
#' All errors are collected and reported together.
#'
#' @param params_file Path to \code{project_config.R}.
#'
#' @return Invisibly returns the \code{project_config} list if all checks pass.
#'   Stops with a formatted summary of all errors otherwise.
#'
#' @export
#' @importFrom readxl read_excel
#' @importFrom stats na.omit
#' @importFrom utils head
mETHYLotest.EPIC.Check_ProjectSettings <- function(params_file = NULL) {

  # ── Helper: human-readable file size ───────────────────────────────────────
  .fmt_bytes <- function(bytes) {
    if (is.na(bytes) || bytes == 0L) return("0 bytes")
    units <- c("bytes", "Kb", "Mb", "Gb", "Tb")
    i     <- min(floor(log(bytes, 1024)), length(units) - 1L)
    paste(sprintf("%.1f", bytes / (1024 ^ i)), units[i + 1L])
  }

  errors   <- character(0)
  warnings <- character(0)

  # ========================================================================
  # 1. LOAD CONFIGURATION
  # ========================================================================

  if (is.null(params_file) || !file.exists(params_file))
    stop("[mETHYLotest Check] Config file must be provided and exist: ",
         params_file)

  config_env <- new.env()
  tryCatch(
    source(params_file, local = config_env),
    error = function(e)
      stop("[mETHYLotest Check] Error sourcing config: ", e$message)
  )

  if (!exists("project_config", envir = config_env))
    stop("[mETHYLotest Check] Config file did not create 'project_config'.")

  cfg <- config_env$project_config

  message("[mETHYLotest Check] ════════════════════════════════════════")
  message("[mETHYLotest Check] Project: ", cfg$project_name)
  message("[mETHYLotest Check] ════════════════════════════════════════")

  # ========================================================================
  # 2. CRITICAL PATHS
  # ========================================================================

  message("[mETHYLotest Check] --- 1. Validating paths ---")

  # Project directory
  if (!is.null(cfg$project_dir) && dir.exists(cfg$project_dir)) {
    message("[mETHYLotest Check] OK  Project dir: ", cfg$project_dir)
  } else {
    errors <- c(errors, paste("Project directory not found:",
                              cfg$project_dir))
  }

  # IDAT directory
  if (!is.null(cfg$idat_dir) && dir.exists(cfg$idat_dir)) {
    message("[mETHYLotest Check] OK  IDAT dir: ", cfg$idat_dir)
  } else {
    errors <- c(errors, paste("IDAT directory not found:", cfg$idat_dir))
  }

  # Phenotype file
  if (!is.null(cfg$pheno_file) && file.exists(cfg$pheno_file)) {
    message("[mETHYLotest Check] OK  Phenotype file: ", cfg$pheno_file)
  } else {
    errors <- c(errors, paste("Phenotype file not found:", cfg$pheno_file))
  }

  # Results directory
  if (!is.null(cfg$res_dir)) {
    if (!dir.exists(cfg$res_dir)) {
      dir.create(cfg$res_dir, recursive = TRUE)
      message("[mETHYLotest Check] OK  Results dir created: ", cfg$res_dir)
    } else {
      message("[mETHYLotest Check] OK  Results dir: ", cfg$res_dir)
    }
  }

  # RDS directory
  if (!is.null(cfg$rds_dir)) {
    if (!dir.exists(cfg$rds_dir)) {
      dir.create(cfg$rds_dir, recursive = TRUE)
      message("[mETHYLotest Check] OK  RDS dir created: ", cfg$rds_dir)
    } else {
      message("[mETHYLotest Check] OK  RDS dir: ", cfg$rds_dir)
    }
  }

  if (length(errors) > 0)
    stop("[mETHYLotest Check] Critical path errors:\n",
         paste0("  - ", errors, collapse = "\n"))

  # ========================================================================
  # 3. PHENOTYPE DATA
  # ========================================================================

  message("[mETHYLotest Check] --- 2. Validating phenotype data ---")

  Pheno <- tryCatch({
    ext <- tools::file_ext(cfg$pheno_file)
    if (ext %in% c("xlsx", "xls"))
      as.data.frame(readxl::read_excel(cfg$pheno_file))
    else
      read.csv(cfg$pheno_file, stringsAsFactors = FALSE)
  }, error = function(e)
    stop("[mETHYLotest Check] Cannot read phenotype file: ", e$message))

  message("[mETHYLotest Check] OK  Phenotype loaded: ",
          nrow(Pheno), " rows, ", ncol(Pheno), " columns.")

  # ── Required columns ─────────────────────────────────────────────────────
  required_cols <- c("Sample_Name", "Sample_Plate", "Sample_Group",
                     "Sentrix_ID", "Sentrix_Position")
  missing_cols  <- setdiff(required_cols, colnames(Pheno))

  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing required column(s):",
                              paste(missing_cols, collapse = ", ")))
    stop("[mETHYLotest Check] Missing columns:\n",
         paste0("  - ", errors, collapse = "\n"))
  }
  message("[mETHYLotest Check] OK  All required columns present.")

  # ── Sample_Name uniqueness ──────────────────────────────────────────────
  names_all <- stats::na.omit(Pheno$Sample_Name)
  names_dup <- unique(names_all[duplicated(names_all)])

  if (length(names_dup) > 0) {
    errors <- c(errors, paste0(
      "Duplicate Sample_Name(s): ",
      paste(utils::head(names_dup, 10L), collapse = ", "),
      if (length(names_dup) > 10L)
        paste0(" ... and ", length(names_dup) - 10L, " more.")))
  } else {
    message("[mETHYLotest Check] OK  All Sample_Name values unique (",
            length(names_all), " entries).")
  }

  # ── NA check on critical columns ────────────────────────────────────────
  for (col in required_cols) {
    na_count <- sum(is.na(Pheno[[col]]) | !nzchar(trimws(Pheno[[col]])))
    if (na_count > 0)
      errors <- c(errors, paste0(
        "Column '", col, "' has ", na_count, " empty/NA value(s)."))
  }

  # ── Sample_Plate validation ─────────────────────────────────────────────
  message("[mETHYLotest Check] --- 3. Validating Sample_Plate ---")

  valid_plates <- c("450K", "EPICv1", "EPICv2", "Mouse")
  plate_vals   <- unique(trimws(Pheno$Sample_Plate))
  bad_plates   <- plate_vals[!plate_vals %in% valid_plates]

  if (length(bad_plates) > 0) {
    errors <- c(errors, paste0(
      "Invalid Sample_Plate value(s): ",
      paste(bad_plates, collapse = ", "),
      ". Expected: ", paste(valid_plates, collapse = ", ")))
  } else {
    message("[mETHYLotest Check] OK  Sample_Plate: ",
            paste(plate_vals, collapse = ", "))
  }

  plate_counts <- table(Pheno$Sample_Plate)
  for (p in names(plate_counts))
    message("[mETHYLotest Check]     ", p, ": ", plate_counts[[p]],
            " samples")

  # ── Sentrix_ID validation ───────────────────────────────────────────────
  message("[mETHYLotest Check] --- 4. Validating Sentrix_ID ---")

  sid_vals  <- stats::na.omit(Pheno$Sentrix_ID)
  bad_sid   <- sid_vals[!grepl("^[0-9]+$", trimws(sid_vals))]

  if (length(bad_sid) > 0) {
    warnings <- c(warnings, paste0(
      length(bad_sid),
      " Sentrix_ID value(s) are not purely numeric: ",
      paste(utils::head(unique(bad_sid), 5L), collapse = ", ")))
  } else {
    n_slides <- length(unique(sid_vals))
    message("[mETHYLotest Check] OK  Sentrix_ID: ",
            n_slides, " unique slide(s).")
  }

  # ── Sentrix_Position validation ─────────────────────────────────────────
  message("[mETHYLotest Check] --- 5. Validating Sentrix_Position ---")

  pos_vals <- stats::na.omit(Pheno$Sentrix_Position)
  bad_pos  <- pos_vals[!grepl("^R[0-9]{2}C[0-9]{2}$", trimws(pos_vals))]

  if (length(bad_pos) > 0) {
    errors <- c(errors, paste0(
      length(bad_pos),
      " Sentrix_Position value(s) do not match 'RxxCxx' (e.g. R06C01):\n",
      paste0("    - ", utils::head(unique(bad_pos), 10L), collapse = "\n"),
      if (length(unique(bad_pos)) > 10L)
        paste0("\n    ... and ", length(unique(bad_pos)) - 10L, " more.")))
  } else {
    message("[mETHYLotest Check] OK  All Sentrix_Position values valid.")
  }

  # ── Construct basenames and check uniqueness ────────────────────────────
  message("[mETHYLotest Check] --- 6. Constructing basenames ---")

  Pheno$.basename <- paste0(
    trimws(Pheno$Sentrix_ID), "_", trimws(Pheno$Sentrix_Position))

  basenames_all  <- stats::na.omit(Pheno$.basename)
  basenames_uniq <- unique(basenames_all)
  basenames_dup  <- unique(basenames_all[duplicated(basenames_all)])

  if (length(basenames_dup) > 0) {
    warnings <- c(warnings, paste0(
      "Duplicate Sentrix_ID + Sentrix_Position combinations: ",
      paste(utils::head(basenames_dup, 5L), collapse = ", ")))
  } else {
    message("[mETHYLotest Check] OK  All basenames unique (",
            length(basenames_uniq), " entries).")
  }

  # ── Sample_Group validation ─────────────────────────────────────────────
  message("[mETHYLotest Check] --- 7. Validating Sample_Group ---")

  group_vals <- unique(trimws(Pheno$Sample_Group))
  message("[mETHYLotest Check] OK  Groups: ",
          paste(group_vals, collapse = ", "),
          " (", length(group_vals), " levels)")

  if (length(group_vals) < 2L)
    warnings <- c(warnings,
                  "Only one Sample_Group level. DMP requires >= 2 groups.")

  group_counts <- table(Pheno$Sample_Group)
  for (g in names(group_counts))
    message("[mETHYLotest Check]     ", g, ": ", group_counts[[g]],
            " samples")

  small_groups <- names(group_counts[group_counts < 3L])
  if (length(small_groups) > 0)
    warnings <- c(warnings, paste0(
      "Group(s) with < 3 samples: ",
      paste(small_groups, collapse = ", "),
      ". Statistical power may be insufficient."))

  # ========================================================================
  # 4. IDAT FILE VERIFICATION
  # ========================================================================

  message("[mETHYLotest Check] --- 8. Checking .idat files ---")

  idat_path <- cfg$idat_dir

  if (length(basenames_uniq) == 0L) {
    warnings <- c(warnings, "No basenames to check for .idat files.")
  } else {
    message("[mETHYLotest Check] Scanning ",
            length(basenames_uniq) * 2L, " .idat files...")

    missing_idat <- character(0)
    total_size   <- 0

    for (bn in basenames_uniq) {
      path_red <- file.path(idat_path, paste0(bn, "_Red.idat"))
      path_grn <- file.path(idat_path, paste0(bn, "_Grn.idat"))

      # Also check in subdirectories (some setups use Sentrix_ID as subfolder)
      sid <- sub("_.*$", "", bn)
      path_red_sub <- file.path(idat_path, sid, paste0(bn, "_Red.idat"))
      path_grn_sub <- file.path(idat_path, sid, paste0(bn, "_Grn.idat"))

      red_found <- file.exists(path_red) || file.exists(path_red_sub)
      grn_found <- file.exists(path_grn) || file.exists(path_grn_sub)

      red_path <- if (file.exists(path_red)) path_red else path_red_sub
      grn_path <- if (file.exists(path_grn)) path_grn else path_grn_sub

      if (red_found) {
        sz <- file.info(red_path)$size
        total_size <- total_size + sz
        message("[mETHYLotest Check]   ", bn, "_Red.idat : ",
                .fmt_bytes(sz))
      } else {
        message("[mETHYLotest Check]   ", bn, "_Red.idat : NOT FOUND")
        missing_idat <- c(missing_idat, paste0(bn, "_Red.idat"))
      }

      if (grn_found) {
        sz <- file.info(grn_path)$size
        total_size <- total_size + sz
        message("[mETHYLotest Check]   ", bn, "_Grn.idat : ",
                .fmt_bytes(sz))
      } else {
        message("[mETHYLotest Check]   ", bn, "_Grn.idat : NOT FOUND")
        missing_idat <- c(missing_idat, paste0(bn, "_Grn.idat"))
      }
    }

    if (length(missing_idat) > 0) {
      errors <- c(errors, paste0(
        "Missing ", length(missing_idat), " .idat file(s):\n",
        paste0("    - ", utils::head(missing_idat, 25L), collapse = "\n"),
        if (length(missing_idat) > 25L)
          paste0("\n    ... and ", length(missing_idat) - 25L, " more.")))
    } else {
      message("[mETHYLotest Check] OK  All ",
              length(basenames_uniq) * 2L,
              " .idat files found (total: ", .fmt_bytes(total_size), ").")
    }
  }

  # ========================================================================
  # 5. CONFIGURATION PARAMETER VALIDATION
  # ========================================================================

  message("[mETHYLotest Check] --- 9. Validating config parameters ---")

  # ── Import parameters ───────────────────────────────────────────────────
  valid_methods <- c("ChAMP", "minfi")
  if (!is.null(cfg$load_method) && !cfg$load_method %in% valid_methods)
    errors <- c(errors, paste0(
      "Invalid load_method: '", cfg$load_method,
      "'. Expected: ", paste(valid_methods, collapse = ", ")))

  # ── Normalisation ───────────────────────────────────────────────────────
  valid_norms <- c("BMIQ", "PBC", "SWAN", "illumina", "none")
  if (!is.null(cfg$norm_method) && !cfg$norm_method %in% valid_norms)
    errors <- c(errors, paste0(
      "Invalid norm_method: '", cfg$norm_method,
      "'. Expected: ", paste(valid_norms, collapse = ", ")))

  if (!is.null(cfg$norm_method) && cfg$norm_method == "SWAN" &&
      !is.null(cfg$load_method) && cfg$load_method == "ChAMP")
    errors <- c(errors,
                "SWAN normalisation requires load_method = 'minfi'.")

  # ── Imputation ──────────────────────────────────────────────────────────
  if (isTRUE(cfg$do_impute)) {
    valid_impute <- c("KNN", "Combine")
    if (!is.null(cfg$impute_method) &&
        !cfg$impute_method %in% valid_impute)
      errors <- c(errors, paste0(
        "Invalid impute_method: '", cfg$impute_method,
        "'. Expected: ", paste(valid_impute, collapse = ", ")))
    message("[mETHYLotest Check] OK  Imputation: ",
            cfg$impute_method, " (k=", cfg$impute_k, ")")
  }

  # ── Cell type correction ────────────────────────────────────────────────
  if (isTRUE(cfg$do_refbase)) {
    message("[mETHYLotest Check] OK  Cell type correction enabled",
            " (champ.refbase, blood only).")
  }

  # ── Batch correction ────────────────────────────────────────────────────
  if (isTRUE(cfg$perform_batch_correction)) {
    if (is.null(cfg$batch_cols) || length(cfg$batch_cols) == 0L) {
      warnings <- c(warnings,
                    "Batch correction enabled but no batch_cols specified.")
    } else {
      missing_batch <- setdiff(cfg$batch_cols, colnames(Pheno))
      if (length(missing_batch) > 0)
        errors <- c(errors, paste0(
          "Batch column(s) not in phenotype file: ",
          paste(missing_batch, collapse = ", ")))
      else
        message("[mETHYLotest Check] OK  Batch correction: ",
                paste(cfg$batch_cols, collapse = " + "),
                " | protect: ", cfg$biological_variable)
    }
  }

  # ── Analysis variable ──────────────────────────────────────────────────
  pheno_col <- if (!is.null(cfg$analysis_pheno)) cfg$analysis_pheno
  else "Sample_Group"

  if (!pheno_col %in% colnames(Pheno))
    errors <- c(errors, paste0(
      "Analysis variable '", pheno_col,
      "' not found in phenotype columns: ",
      paste(colnames(Pheno), collapse = ", ")))
  else
    message("[mETHYLotest Check] OK  Analysis variable: '", pheno_col, "'")

  # compare.group validation
  if (!is.null(cfg$compare_group) && length(cfg$compare_group) == 2L) {
    if (pheno_col %in% colnames(Pheno)) {
      avail_groups <- unique(as.character(Pheno[[pheno_col]]))
      bad_grp <- setdiff(cfg$compare_group, avail_groups)
      if (length(bad_grp) > 0)
        errors <- c(errors, paste0(
          "compare_group value(s) not found in '", pheno_col, "': ",
          paste(bad_grp, collapse = ", "),
          ". Available: ", paste(avail_groups, collapse = ", ")))
      else
        message("[mETHYLotest Check] OK  Comparison: ",
                cfg$compare_group[1], " vs ", cfg$compare_group[2])
    }
  }

  # ── DMP parameters ─────────────────────────────────────────────────────
  if (isTRUE(cfg$do_dmp)) {
    valid_adjust <- c("BH", "bonferroni", "holm", "hochberg",
                      "hommel", "BY", "fdr", "none")
    if (!is.null(cfg$dmp_adjust_method) &&
        !cfg$dmp_adjust_method %in% valid_adjust)
      errors <- c(errors, paste0(
        "Invalid dmp_adjust_method: '", cfg$dmp_adjust_method, "'."))
    message("[mETHYLotest Check] OK  DMP: adjPVal=",
            cfg$dmp_adj_p_val, ", method=", cfg$dmp_adjust_method)
  }

  # ── DMR parameters ─────────────────────────────────────────────────────
  if (isTRUE(cfg$do_dmr)) {
    valid_dmr <- c("Bumphunter", "DMRcate", "ProbeLasso")
    if (!is.null(cfg$dmr_method) && !cfg$dmr_method %in% valid_dmr)
      errors <- c(errors, paste0(
        "Invalid dmr_method: '", cfg$dmr_method,
        "'. Expected: ", paste(valid_dmr, collapse = ", ")))
    message("[mETHYLotest Check] OK  DMR: method=", cfg$dmr_method,
            ", minProbes=", cfg$dmr_min_probes)
  }

  # ── GSEA dependency ────────────────────────────────────────────────────
  if (isTRUE(cfg$do_gsea) && !isTRUE(cfg$do_dmp))
    warnings <- c(warnings,
                  "GSEA enabled but DMP is disabled. GSEA requires DMP results.")

  # ── CNA parameters ─────────────────────────────────────────────────────
  # Only check CNA control group if CNA is enabled
  if (isTRUE(cfg$do_cna)) {
    ctrl_grp <- cfg$cna_control_group
    if (!is.null(ctrl_grp) && nchar(ctrl_grp) > 0) {
      groups <- unique(as.character(pheno[[cfg$col_sample_group]]))
      if (!ctrl_grp %in% groups) {
        errors <- c(errors, sprintf(
          "CNA control group '%s' not found in '%s'. Available: %s",
          ctrl_grp, cfg$col_sample_group,
          paste(groups, collapse = ", ")))
      }
    }
  }

  # ── System resources ────────────────────────────────────────────────────
  if (!is.null(cfg$num_cores)) {
    total <- parallel::detectCores()
    if (cfg$num_cores > total)
      warnings <- c(warnings, paste0(
        "num_cores (", cfg$num_cores,
        ") exceeds detected cores (", total, ")."))
    message("[mETHYLotest Check] OK  Cores: ", cfg$num_cores,
            " / ", total, " available")
  }

  # ========================================================================
  # 6. FINAL SUMMARY
  # ========================================================================

  message("[mETHYLotest Check] ════════════════════════════════════════")

  # Print warnings
  if (length(warnings) > 0) {
    message("[mETHYLotest Check] ", length(warnings), " WARNING(s):")
    for (w in warnings)
      message("[mETHYLotest Check]   ! ", w)
  }

  # Print errors and stop
  if (length(errors) > 0) {
    stop("[mETHYLotest Check] ", length(errors), " ERROR(s):\n",
         paste0("  - ", errors, collapse = "\n"))
  }

  # Pipeline summary
  steps <- c(
    "champ.load",
    if (isTRUE(cfg$do_impute))                   "champ.impute",
    if (isTRUE(cfg$normalize_data))               paste0("champ.norm (",
                                                         cfg$norm_method, ")"),
    if (isTRUE(cfg$do_refbase))                   "champ.refbase",
    if (isTRUE(cfg$perform_batch_correction))     paste0("ComBat (",
                                                         paste(cfg$batch_cols,
                                                               collapse="+"),
                                                         ")"),
    if (isTRUE(cfg$do_dmp))                       "champ.DMP",
    if (isTRUE(cfg$do_dmr))                       paste0("champ.DMR (",
                                                         cfg$dmr_method, ")"),
    if (isTRUE(cfg$do_block))                     "champ.Block",
    if (isTRUE(cfg$do_gsea))                      "champ.GSEA",
    if (isTRUE(cfg$do_cna))                       "champ.CNA"
  )

  message("[mETHYLotest Check] All checks passed.")
  message("[mETHYLotest Check] Samples:  ", nrow(Pheno))
  message("[mETHYLotest Check] Arrays:   ",
          paste(names(plate_counts), collapse = " + "))
  message("[mETHYLotest Check] Pipeline: ",
          paste(steps, collapse = " -> "))
  message("[mETHYLotest Check] ════════════════════════════════════════")

  invisible(cfg)
}
