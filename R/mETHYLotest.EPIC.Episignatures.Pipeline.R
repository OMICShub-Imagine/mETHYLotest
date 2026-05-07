#' Episignature Scoring Pipeline (Standalone)
#'
#' @description
#' Computes Z-score-based episignature scores for test samples against
#' an internal control cohort. This is a standalone analysis that can be
#' run independently from the main EPIC pipeline.
#'
#' The function:
#' \enumerate{
#'   \item Loads IDAT data for test samples (from an existing project)
#'   \item Loads internal control IDATs (shipped with the package)
#'   \item Merges and harmonizes arrays
#'   \item For each known episignature, computes M-value Z-scores
#'   \item Generates per-sample scoring with coverage metrics
#'   \item Exports results (CSV, RDS) and an interactive HTML report
#' }
#'
#' @param project_directory Path to an existing mETHYLotest EPIC project.
#'   Must contain \code{Results/project_config.R}.
#'
#' @return Invisibly returns the results data.frame.
#' @export
#' @importFrom ChAMP champ.load
#' @importFrom readxl read_excel
#' @importFrom stats sd pnorm
#' @importFrom utils write.csv read.csv
#' @importFrom rmarkdown render
mETHYLotest.EPIC.Episignatures <- function(project_directory) {

  # ========================================================================
  # 1. PROJECT SETUP
  # ========================================================================

  config_file <- normalizePath(
    file.path(project_directory, "Results", "project_config.R"))
  if (!file.exists(config_file))
    stop("[Episignatures] Config not found: ", config_file)

  source(config_file)
  cfg <- project_config

  message("[Episignatures] ========================================")
  message("[Episignatures] Project: ", cfg$project_name)
  message("[Episignatures] ========================================")

  # Output directory
  episig_dir <- file.path(cfg$res_dir, "Episignatures_Scoring")
  if (!dir.exists(episig_dir))
    dir.create(episig_dir, recursive = TRUE)

  # ========================================================================
  # 2. PHENOTYPE
  # ========================================================================

  if (!file.exists(cfg$pheno_file))
    stop("[Episignatures] Phenotype not found: ", cfg$pheno_file)

  Pheno <- as.data.frame(readxl::read_excel(cfg$pheno_file))

  col_name  <- cfg$col_sample_name
  col_plate <- cfg$col_sample_plate
  col_group <- cfg$col_sample_group
  col_sid   <- cfg$col_sentrix_id
  col_pos   <- cfg$col_sentrix_pos

  message("[Episignatures] Samples: ", nrow(Pheno))
  message("[Episignatures] Groups: ",
          paste(unique(Pheno[[col_group]]), collapse = ", "))

  # ========================================================================
  # 3. ARRAY TYPE DETECTION
  # ========================================================================

  plate_to_champ <- function(plate) {
    switch(toupper(trimws(plate)),
           "450K"   = "450K",
           "EPICV1" = "EPICv1",
           "EPICV2" = "EPICv2",
           stop("[Episignatures] Unknown array: ", plate))
  }

  plate_values <- unique(Pheno[[col_plate]])
  idat_dir     <- cfg$idat_dir

  # ChAMP CSV writer
  write_champ_pheno <- function(df) {
    out <- df
    renames <- c(setNames("Sample_Name",  col_name),
                 setNames("Slide",        col_sid),
                 setNames("Array",        col_pos),
                 setNames("Sample_Group", col_group))
    for (from in names(renames)) {
      to <- renames[[from]]
      if (from %in% colnames(out) && from != to)
        colnames(out)[colnames(out) == from] <- to
    }
    csv_tmp <- file.path(idat_dir, "Pheno.csv")
    write.table(out, csv_tmp, row.names = FALSE,
                quote = FALSE, sep = ",")
    csv_tmp
  }

  # ========================================================================
  # 4. LOAD TEST DATA
  # ========================================================================

  message("[Episignatures] Loading test samples...")

  # Handle single-sample array types by duplicating temporarily
  duplicated_samples <- character(0)

  for (pv in plate_values) {
    idx <- which(Pheno[[col_plate]] == pv)
    if (length(idx) == 1L) {
      message("[Episignatures] Single sample for ", pv,
              ". Duplicating temporarily.")
      dup_row <- Pheno[idx, ]
      dup_name <- paste0(dup_row[[col_name]], "_dup")
      dup_row[[col_name]] <- dup_name
      Pheno <- rbind(Pheno, dup_row)
      duplicated_samples <- c(duplicated_samples, dup_name)
    }
  }

  loads_list <- list()
  for (pv in plate_values) {
    champ_at  <- plate_to_champ(pv)
    pheno_sub <- Pheno[Pheno[[col_plate]] == pv, , drop = FALSE]
    write_champ_pheno(pheno_sub)

    loads_list[[pv]] <- ChAMP::champ.load(
      directory = idat_dir, arraytype = champ_at, method = "ChAMP")

    message("[Episignatures] ", pv, ": ",
            nrow(loads_list[[pv]]$beta), " CpGs, ",
            ncol(loads_list[[pv]]$beta), " samples")
  }

  # Harmonize
  if (length(loads_list) == 1L) {
    myLoad <- loads_list[[1L]]
  } else {
    myLoad <- mETHYLotest.utils.HarmonizeArrays(loads_list)
  }

  rownames(myLoad$pd) <- myLoad$pd[["Sample_Name"]]

  # Clean temp CSV
  csv_tmp <- file.path(idat_dir, "Pheno.csv")
  if (file.exists(csv_tmp)) file.remove(csv_tmp)

  # ========================================================================
  # 5. LOAD INTERNAL CONTROLS
  # ========================================================================

  message("[Episignatures] Loading internal controls...")

  ctl_dir <- system.file("extdata", "idats_ctl", package = "mETHYLotest")
  if (!nzchar(ctl_dir))
    stop("[Episignatures] Control IDATs not found in package.")

  myLoad_ctl <- ChAMP::champ.load(
    directory = ctl_dir, arraytype = "EPICv1", method = "ChAMP")

  message("[Episignatures] Controls: ",
          ncol(myLoad_ctl$beta), " samples, ",
          nrow(myLoad_ctl$beta), " CpGs")

  # ========================================================================
  # 6. MERGE TEST + CONTROLS
  # ========================================================================

  message("[Episignatures] Merging datasets...")

  common_probes <- intersect(rownames(myLoad$beta),
                             rownames(myLoad_ctl$beta))
  message("[Episignatures] Common probes: ",
          format(length(common_probes), big.mark = ","))

  beta_combined <- cbind(
    myLoad$beta[common_probes, , drop = FALSE],
    myLoad_ctl$beta[common_probes, , drop = FALSE])

  control_samples <- colnames(myLoad_ctl$beta)
  test_samples    <- setdiff(colnames(myLoad$beta), duplicated_samples)

  message("[Episignatures] Test: ", length(test_samples),
          " | Controls: ", length(control_samples))

  # ========================================================================
  # 7. LOAD EPISIGNATURE DEFINITIONS
  # ========================================================================

  episig_tsv <- system.file("extdata", "episignatures.tsv",
                            package = "mETHYLotest")
  if (!nzchar(episig_tsv))
    stop("[Episignatures] episignatures.tsv not found in package.")

  episignatures <- read.delim(episig_tsv, header = TRUE,
                              stringsAsFactors = FALSE)
  message("[Episignatures] ", nrow(episignatures),
          " episignature(s) to score.")

  # ========================================================================
  # 8. SCORING LOOP
  # ========================================================================

  results_list <- list()

  for (i in seq_len(nrow(episignatures))) {
    sig_name <- episignatures$Abbreviation[i]

    sig_file <- system.file("extdata", "signatures",
                            paste0(sig_name, ".txt"),
                            package = "mETHYLotest")
    if (!file.exists(sig_file)) {
      warning("[Episignatures] File not found: ", sig_name)
      next
    }

    sig_cpgs    <- readLines(sig_file)
    sig_cpgs    <- trimws(sig_cpgs[nzchar(sig_cpgs)])
    total_count <- length(sig_cpgs)
    common_cpgs <- intersect(rownames(beta_combined), sig_cpgs)
    found_count <- length(common_cpgs)

    if (found_count < 3L) {
      warning("[Episignatures] < 3 probes for ", sig_name, ". Skipped.")
      next
    }

    # Subset
    dat_beta <- beta_combined[common_cpgs, , drop = FALSE]

    # Beta -> M-values
    eps <- 1e-5
    dat_clipped <- pmin(pmax(dat_beta, eps), 1 - eps)
    dat_m <- log2(dat_clipped / (1 - dat_clipped))

    # Control stats (M-values)
    ctrl_m  <- dat_m[, control_samples, drop = FALSE]
    mu_m    <- rowMeans(ctrl_m, na.rm = TRUE)
    sd_m    <- apply(ctrl_m, 1, sd, na.rm = TRUE)
    sd_m[sd_m == 0] <- 1e-6

    # Control mean beta (for delta-beta)
    ctrl_beta <- dat_beta[, control_samples, drop = FALSE]
    mu_beta   <- rowMeans(ctrl_beta, na.rm = TRUE)

    # Score each test sample
    for (sid in test_samples) {
      val_m    <- dat_m[, sid]
      val_beta <- dat_beta[, sid]

      z <- (val_m - mu_m) / sd_m
      p <- 2 * pnorm(-abs(z))
      db <- val_beta - mu_beta

      sig_idx    <- which(p < 0.05)
      n_sig      <- length(sig_idx)
      pct_sig    <- round(100 * n_sig / found_count, 2)
      coverage   <- round(100 * found_count / total_count, 2)
      global_sc  <- round(100 * n_sig / total_count, 2)

      # Details
      if (n_sig > 0L) {
        pnames <- names(p)[sig_idx]
        details <- paste(
          paste0(pnames, ": p=",
                 formatC(p[sig_idx], format = "e", digits = 1),
                 ", dB=", round(db[sig_idx], 3)),
          collapse = "; ")
      } else {
        details <- "None"
      }

      results_list[[length(results_list) + 1L]] <- data.frame(
        Sample         = sid,
        Signature      = sig_name,
        Total_Probes   = total_count,
        Found          = found_count,
        Missing        = total_count - found_count,
        Coverage_Pct   = coverage,
        Significant    = n_sig,
        Pct_Significant = pct_sig,
        Global_Score   = global_sc,
        Details        = details,
        stringsAsFactors = FALSE
      )
    }

    message("[Episignatures]   ", sig_name,
            " | probes: ", found_count, "/", total_count,
            " (", round(100 * found_count / total_count), "%)")
  }

  # ========================================================================
  # 9. AGGREGATE AND EXPORT
  # ========================================================================

  all_results <- do.call(rbind, results_list)

  csv_path <- file.path(episig_dir, "episignature_scores.csv")
  rds_path <- file.path(episig_dir, "episignature_scores.rds")

  write.csv(all_results, csv_path, row.names = FALSE)
  saveRDS(all_results, rds_path)

  message("[Episignatures] Results saved: ", episig_dir)
  message("[Episignatures] Samples scored: ",
          length(test_samples))
  message("[Episignatures] Signatures scored: ",
          length(unique(all_results$Signature)))

  # ========================================================================
  # 10. HTML REPORT
  # ========================================================================

  message("[Episignatures] Generating report...")

  rmd_template <- system.file("rmarkdown",
                              "mETHYLotest.EPIC.Episignature_Report.Rmd",
                              package = "mETHYLotest")
  if (!nzchar(rmd_template)) {
    warning("[Episignatures] Report template not found. Skipping.")
  } else {
    html_path <- file.path(episig_dir,
                           "mETHYLotest-Episignature_Report.html")
    tryCatch({
      rmarkdown::render(
        input       = rmd_template,
        output_file = basename(html_path),
        output_dir  = episig_dir,
        params      = list(rds_file = rds_path),
        envir       = new.env(),
        quiet       = TRUE)
      message("[Episignatures] Report: ", html_path)
    }, error = function(e)
      warning("[Episignatures] Report failed: ", e$message))
  }

  # ========================================================================
  # 11. SUMMARY
  # ========================================================================

  message("[Episignatures] ========================================")
  message("[Episignatures] Complete!")

  # Top hits
  top <- all_results[order(-all_results$Global_Score), ]
  top <- top[!duplicated(paste(top$Sample, top$Signature)), ]
  top5 <- head(top[top$Global_Score > 0, ], 10L)

  if (nrow(top5) > 0L) {
    message("[Episignatures] Top hits:")
    for (r in seq_len(nrow(top5)))
      message("[Episignatures]   ", top5$Sample[r],
              " x ", top5$Signature[r],
              " | Score: ", top5$Global_Score[r],
              "% | Sig: ", top5$Significant[r],
              "/", top5$Found[r])
  }

  message("[Episignatures] ========================================")

  invisible(all_results)
}
