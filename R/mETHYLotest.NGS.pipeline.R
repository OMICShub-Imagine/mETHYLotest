#' WGBS / Long Read Methylation Pipeline
#'
#' @description
#' Complete methylation analysis pipeline for WGBS and Long-Read data.
#' Handles data import, iterative QC with interactive Shiny UI,
#' control genome monitoring, comparative analysis, multi-scenario
#' differential methylation (CpG-level and region-level), segmentation,
#' signature generation/validation, genomic annotation, and HTML reporting.
#'
#' Passing \code{"demo"} as \code{project_directory} launches a
#' self-contained demonstration using synthetic CX reports bundled
#' with the package.
#'
#' @param project_directory Character. Path to the mETHYLotest project directory.
#'   If empty or invalid, launches the Project Setup UI.
#'   If \code{"demo"}, launches demo mode with bundled synthetic data.
#'
#' @return A named list (invisibly) containing:
#'   \describe{
#'     \item{meth}{methylBase object (united, after QC)}
#'     \item{diff_results}{Named list of methylDiff objects (one per scenario)}
#'     \item{tiles_results}{Named list of tiling window results (or NULL)}
#'     \item{seg_results}{Segmentation results (or NULL)}
#'     \item{validation_results}{Named list of signature validation results (or NULL)}
#'     \item{annotated}{Annotated DMC data (or NULL)}
#'     \item{config}{Pipeline configuration list}
#'     \item{profile}{Performance profiling data.frame}
#'   }
#'
#' @details
#' The pipeline proceeds through the following steps:
#' \enumerate{
#'   \item Configuration validation and phenotype import
#'   \item Data import with smart caching (skips if \code{myobj_raw.rds} exists)
#'   \item Control genome pre-computation (pUC19, Lambda, E. coli, PhiX)
#'   \item Iterative QC with interactive Shiny dashboard
#'   \item Comparative analysis (correlation, hierarchical clustering)
#'   \item Multi-scenario differential methylation (Basic/Batch/Covariate/Full)
#'   \item DMC export with per-group means (Excel + BED)
#'   \item Tiling windows / DMR analysis (optional)
#'   \item Methylation segmentation via HMM (optional)
#'   \item Signature generation (top 500 DMCs by combined score)
#'   \item DiffMeth Explorer UI (interactive, user signatures)
#'   \item Signature validation (SVM/PCA/Silhouette)
#'   \item Validation UI + HTML report
#'   \item Genomic annotation (annotatr)
#'   \item Final HTML report generation
#'   \item Performance profiling export
#' }
#'
#' Performance profiling (wall-clock time, peak RAM, disk footprint)
#' is recorded per step and exported as \code{Pipeline_Performance.xlsx}.
#'
#' Intermediate objects are cached as \code{.rds} files in
#' \code{data/interim/}. On re-runs, import is skipped entirely.
#'
#' @export
#'
#' @importFrom methylKit methRead getMethylationStats getCoverageStats
#'   clusterSamples filterByCoverage unite calculateDiffMeth
#'   getMethylDiff getSampleID getData reorganize percMethylation
#'   tileMethylCounts methSeg getCorrelation
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx
#' @importFrom parallel detectCores
#' @importFrom grDevices png pdf dev.off colorRampPalette
#' @importFrom stats aggregate
#' @importFrom utils object.size
#' @importFrom rmarkdown render
mETHYLotest.NGS.pipeline <- function(project_directory = "") {

  # ── Initialize profiling ──
  .timer <- new.env(parent = emptyenv())
  .timer$steps <- list()
  .timer$start_total <- proc.time()

  .start_step <- function(name) {
    .timer$current_step <<- name
    .timer$current_start <<- proc.time()
    tryCatch(gc(reset = TRUE), error = function(e) NULL)
  }

  .end_step <- function() {
    elapsed <- (proc.time() - .timer$current_start)[["elapsed"]]

    # Safe RAM measurement
    peak_ram <- tryCatch({
      gc_info <- gc()
      # gc() returns a matrix with "Ncells" and "Vcells" rows
      # Column names vary by R version
      if ("max used (Mb)" %in% colnames(gc_info)) {
        sum(gc_info[, "max used (Mb)"])
      } else if (ncol(gc_info) >= 6) {
        sum(gc_info[, 6])  # 6th column is typically max used Mb
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)

    .timer$steps[[.timer$current_step]] <<- list(
      elapsed_sec = round(elapsed, 2),
      peak_ram_mb = round(peak_ram, 1)
    )
    message(sprintf("[mETHYLotest] [PERF] %s: %.1fs | Peak RAM: %s",
                    .timer$current_step, elapsed,
                    if (is.na(peak_ram)) "N/A" else sprintf("%.0f MB", peak_ram)))
  }

  # ── Demo mode ──
  if (identical(tolower(trimws(project_directory)), "demo")) {
    .start_step("[DEMO] Setup")
    message("[mETHYLotest] ========================================")
    message("[mETHYLotest]           DEMO MODE")
    message("[mETHYLotest] ========================================")

    # Locate demo data
    demo_base <- system.file("extdata", "NGS_CXreport",
                             package = "mETHYLotest")
    if (demo_base == "")
      stop("[mETHYLotest] Demo data not found in package.")

    demo_excel <- file.path(demo_base, "demo_samples.xlsx")
    demo_cx    <- file.path(demo_base, "synthetic_CXreports")

    if (!file.exists(demo_excel))
      stop("[mETHYLotest] demo_samples.xlsx not found.")
    if (!dir.exists(demo_cx))
      stop("[mETHYLotest] synthetic_CXreports/ not found.")

    # Read and resolve DEMO:// paths
    pheno <- as.data.frame(readxl::read_excel(demo_excel))
    pheno$Path <- vapply(pheno$Path, function(p) {
      if (grepl("^DEMO://", p)) {
        fname <- sub("^DEMO://", "", p)
        file.path(demo_cx, fname)
      } else {
        p
      }
    }, character(1))

    # Verify all files exist
    ok <- file.exists(pheno$Path)
    message("[mETHYLotest] Demo samples: ", sum(ok), "/", nrow(pheno), " files found")
    if (any(!ok)) {
      warning("[mETHYLotest] Missing files:\n  ",
              paste(pheno$Path[!ok], collapse = "\n  "))
    }

    # Save resolved Excel to temp
    resolved_excel <- file.path(tempdir(), "mETHYLotest_demo_samples.xlsx")
    writexl::write_xlsx(pheno, resolved_excel)

    message("[mETHYLotest] Resolved metadata: ", resolved_excel)
    message("[mETHYLotest] Opening Project UI with demo data...")
    .end_step()

    project_directory <- mETHYLotest.NGS.ProjectUI(
      prefill_pheno = resolved_excel)
  }

  # ========================================================================
  # STARTUP
  # ========================================================================

  if (!nzchar(project_directory) || !dir.exists(project_directory))
    project_directory <- mETHYLotest.NGS.ProjectUI()

  config_path <- normalizePath(
    file.path(project_directory, "Results", "project_config.R"))
  if (!file.exists(config_path))
    stop("[mETHYLotest] Config not found: ", config_path)

  # ========================================================================
  # 1. PHENOTYPE IMPORT
  # ========================================================================
  .start_step("Import and filtering")

  # Validates config + phenotype + data files in one call
  message("[mETHYLotest] Checking config file")
  tryCatch({
    cfg <- mETHYLotest.NGS.Check_ProjectSettings(config_path)
  }, error = function(e){
    stop("[mETHYLotest] [CRITICAL] Unable to launch analysis : ", e)
  })

  # Load phenotype (already validated)
  Pheno <- as.data.frame(readxl::read_excel(cfg$pheno_file))
  
  # Save the phenotype table to a CSV for the Web App to read
  utils::write.csv(Pheno, file.path(res_dir, "Samples_Phenotype.csv"), row.names = FALSE)

  SampleIds <- as.list(as.character(Pheno[[cfg$col_sampleID]]))
  SamplePaths     <- as.list(Pheno[[cfg$col_file_path]])
  SampleTreatment <- as.vector(
    as.numeric(as.character(Pheno[[cfg$col_treatment]])))

  # ========================================================================
  # 2. DATA IMPORT (Smart Load)
  # ========================================================================

  rds_dir    <- cfg$rds_dir
  res_dir    <- cfg$res_dir
  interim_dir <- rds_dir
  output_rds <- file.path(rds_dir, "myobj_raw.rds")

  if (file.exists(output_rds)) {
    message("[mETHYLotest] Loading existing raw object: ", output_rds)
    myobj <- readRDS(output_rds)
  } else {
    message("[mETHYLotest] Importing from raw files...")

    # UI writes string "NA", methylKit expects R NA
    db_val <- cfg$dbtype
    if (is.character(db_val) && toupper(trimws(db_val)) %in% c("NA", "NONE", ""))
      db_val <- NA

    myobj <- methylKit::methRead(
      location   = SamplePaths,
      sample.id  = SampleIds,
      assembly   = cfg$assembly,
      treatment  = SampleTreatment,
      pipeline   = cfg$pipeline,
      context    = cfg$context,
      dbtype     = db_val,
      header     = cfg$header,
      skip       = cfg$skip,
      sep        = cfg$sep,
      resolution = cfg$resolution,
      mincov     = cfg$min_coverage
    )

    if (is.list(cfg$pipeline) && cfg$coord_offset != 0L) {
      message("[mETHYLotest] Applying offset: ", cfg$coord_offset)
      myobj <- mETHYLotest.NGS.ApplyOffsetToObj(
        myobj, offset = as.integer(cfg$coord_offset))
    }

    if (isTRUE(cfg$save_raw_obj)) {
      if (!dir.exists(rds_dir))
        dir.create(rds_dir, recursive = TRUE)
      saveRDS(myobj, output_rds)
    }
  }

  .end_step()

  # ========================================================================
  # 3. ITERATIVE QC LOOP
  # ========================================================================

  loop_iter     <- 1L
  current_cov   <- if (!is.null(cfg$qc_lo_count)) cfg$qc_lo_count else cfg$min_coverage
  current_hi    <- if (!is.null(cfg$qc_hi_perc)) cfg$qc_hi_perc else 99.9
  current_lo_p  <- cfg$qc_lo_perc
  excluded_ids  <- c()

  all_chrs      <- unique(methylKit::getData(myobj[[1]])$chr)
  kept_chrs     <- all_chrs

  # ── Directory for interim RDS files ──
  #interim_dir <- file.path(cfg$project_dir, "data", "interim")
  if (!dir.exists(interim_dir)) dir.create(interim_dir, recursive = TRUE)

  message("[mETHYLotest] === Iterative QC ===")

  .start_step("QC")

  message("[mETHYLotest] Interim data directory: ", interim_dir)

  # ============================================══════════════════════
  # PRE-COMPUTE: Control genome data from ALL original samples
  # ============================================══════════════════════
  control_chr_patterns <- c(
    # pUC19 (METHYLATED control — expected ~100%)
    "pUC19", "chrPUC", "NC_001773",
    # Lambda phage (unmethylated control — expected ~0%)
    "lambda", "Lambda", "chrL", "J02459", "NC_001416",
    # E. coli K-12 (UNMETHYLATED control — expected ~0%)
    "NC_000913",
    # PhiX (sequencing control)
    "phiX", "PhiX", "NC_001422",
    # T7 phage
    "NC_001604"
  )

  # Define expected methylation per control type
  # "high" = methylated control (pUC19), "low" = unmethylated control   # CpG Context
  control_expected <- list(
    "pUC19"      = "high",
    "chrPUC"     = "high",
    "NC_001773"  = "high",
    "lambda"     = "low",
    "Lambda"     = "low",
    "chrL"       = "low",
    "J02459"     = "low",
    "NC_001416"  = "low",
    "NC_000913"  = "low",
    "phiX"       = "low",
    "PhiX"       = "low",
    "NC_001422"  = "low",
    "NC_001604"  = "low"
  )

  all_controls_precomputed <- tryCatch({
    ctrl_list <- lapply(seq_along(myobj), function(i) {
      d   <- methylKit::getData(myobj[[i]])
      sid <- methylKit::getSampleID(myobj)[[i]]

      is_ctrl <- Reduce(`|`, lapply(control_chr_patterns, function(pat) {
        grepl(pat, d$chr, ignore.case = TRUE)
      }))

      ctrl_rows <- d[is_ctrl, ]
      if (nrow(ctrl_rows) == 0L) return(NULL)

      agg <- aggregate(
        cbind(numCs, coverage) ~ chr,
        data = ctrl_rows,
        FUN  = sum
      )
      agg$MethylationPercentage <- 100 * agg$numCs / agg$coverage

      cpg_counts <- as.data.frame(table(ctrl_rows$chr),
                                  stringsAsFactors = FALSE)
      colnames(cpg_counts) <- c("chr", "CpG_Count")
      merged <- merge(agg, cpg_counts, by = "chr")

      # Determine expected type for each chromosome
      expected_type <- vapply(merged$chr, function(ch) {
        for (pat in names(control_expected)) {
          if (grepl(pat, ch, ignore.case = TRUE)) {
            return(control_expected[[pat]])
          }
        }
        return("unknown")
      }, character(1))

      data.frame(
        Sample                = sid,
        Chromosome            = merged$chr,
        CpG_Count             = merged$CpG_Count,
        Coverage              = merged$coverage,
        MethylationPercentage = merged$MethylationPercentage,
        ExpectedMeth          = expected_type,
        stringsAsFactors      = FALSE
      )
    })

    result <- do.call(rbind, ctrl_list)
    if (!is.null(result) && nrow(result) > 0) {
      message("[mETHYLotest] Pre-computed control genome data for ",
              length(unique(result$Sample)), " samples, ",
              length(unique(result$Chromosome)), " control chr(s): ",
              paste(unique(result$Chromosome), collapse = ", "))
    } else {
      message("[mETHYLotest] No control genome chromosomes detected")
    }
    result
  }, error = function(e) {
    message("[mETHYLotest] Warning: could not pre-compute control data: ", e$message)
    NULL
  })

  # ============================================══════════════════════

  # ── Save original sample IDs (before any exclusion) ──
  original_sample_ids <- methylKit::getSampleID(myobj)
  message("[mETHYLotest] Original samples: ", length(original_sample_ids),
          " (", paste(original_sample_ids, collapse = ", "), ")")

  repeat {
    # ── Build unique QC directory (don't overwrite previous iteration) ──
    qc_dir_base <- file.path(res_dir, "QC", paste0("cov_", current_cov))
    if (dir.exists(qc_dir_base)) {
      suffix <- 2L
      repeat {
        qc_dir <- paste0(qc_dir_base, "_iter", suffix)
        if (!dir.exists(qc_dir)) break
        suffix <- suffix + 1L
      }
      message("[mETHYLotest] QC dir already exists, using: ", basename(qc_dir))
    } else {
      qc_dir <- qc_dir_base
    }

    message("[mETHYLotest] QC iteration ", loop_iter,
            " | cov=", current_cov, "x | hi=", current_hi, "%")

    temp_obj <- myobj

    # ── Filter chromosomes ──
    if (length(kept_chrs) < length(all_chrs)) {
      filtered <- lapply(temp_obj, function(s) {
        d <- methylKit::getData(s)
        s[d$chr %in% kept_chrs, ]
      })
      temp_obj <- new("methylRawList", filtered,
                      treatment = temp_obj@treatment)
    }

    # ── Exclude samples (manual + auto) ──
    if (length(excluded_ids) > 0L) {
      all_ids <- methylKit::getSampleID(temp_obj)
      keep    <- setdiff(all_ids, excluded_ids)
      if (length(keep) == 0L) stop("All samples excluded!")
      keep_tx <- SampleTreatment[match(keep,
                                       methylKit::getSampleID(myobj))]
      temp_obj <- methylKit::reorganize(temp_obj,
                                        sample.ids = keep,
                                        treatment = keep_tx)
    }

    # ── Filter coverage ──
    temp_filt <- methylKit::filterByCoverage(
      temp_obj,
      lo.count = current_cov,
      lo.perc  = current_lo_p,
      hi.count = NULL,
      hi.perc  = current_hi)

    rm(temp_obj); gc()

    # ============================================══════════════════
    # Remove samples that became empty after coverage filtering
    # ============================================══════════════════
    filt_nrows <- vapply(temp_filt, nrow, integer(1))
    empty_mask <- filt_nrows == 0L

    if (any(empty_mask)) {
      empty_now <- methylKit::getSampleID(temp_filt)[empty_mask]
      message("[mETHYLotest] Removing ", length(empty_now),
              " empty sample(s) after cov>=", current_cov,
              "x filter: ", paste(empty_now, collapse = ", "))

      excluded_ids <- unique(c(excluded_ids, empty_now))

      keep_mask <- !empty_mask
      if (sum(keep_mask) == 0L) {
        stop("[mETHYLotest] All samples are empty after coverage filtering at ",
             current_cov, "x! Lower the coverage threshold.")
      }

      remaining_tx <- temp_filt@treatment[keep_mask]
      if (length(unique(remaining_tx)) < 2L) {
        stop("[mETHYLotest] Only one treatment group remains after removing ",
             "empty samples. Cannot proceed with differential analysis. ",
             "Lower the coverage threshold or check your data.")
      }

      keep_ids <- methylKit::getSampleID(temp_filt)[keep_mask]
      temp_filt <- methylKit::reorganize(
        temp_filt,
        sample.ids = keep_ids,
        treatment  = remaining_tx)

      message("[mETHYLotest] Continuing with ", length(keep_ids),
              " samples: ", paste(keep_ids, collapse = ", "))
    }
    # ============================================══════════════════

    if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

    # ── Run QC analysis ──
    qc_data <- mETHYLotest.NGS.QC(
      methyl_obj      = temp_filt,
      output_base_dir = qc_dir,
      chromosomes     = kept_chrs,
      current_min_cov = current_cov,
      unite_destrand  = cfg$unite_destrand,
      save_summary    = TRUE)

    # ── Treatment lookup from Pheno (active samples only) ──
    active_treatments <- setNames(
      as.integer(as.character(Pheno[[cfg$col_treatment]])),
      as.character(Pheno[[cfg$col_sampleID]])
    )[methylKit::getSampleID(temp_filt)]

    # ── Enrich global meth with treatment group ──
    df_gm <- qc_data$df_global_meth
    if (!is.null(df_gm)) {
      df_gm$Group <- vapply(as.character(df_gm$Sample), function(sid) {
        tx <- active_treatments[sid]
        if (is.na(tx)) "Unknown"
        else if (tx == 0) "Control (0)"
        else "Case (1)"
      }, character(1))
    }

    # ── Launch QC UI ──
    ui_res <- mETHYLotest.NGS.QC.UI(
      df_meth           = qc_data$df_meth,
      df_controls       = qc_data$df_controls,
      df_stats          = qc_data$df_stats,
      df_global_meth    = df_gm,
      current_cov       = current_cov,
      current_hi_perc   = current_hi,
      current_lo_perc   = current_lo_p,
      current_excluded  = excluded_ids,
      all_chrs          = all_chrs,
      current_chrs_kept = kept_chrs,
      sample_treatments = active_treatments)

    if (ui_res$action == "proceed") {

      # ============================================════════════════
      # Save all interim data to data/interim/
      # ============================================════════════════
      analyzed_sample_ids <- methylKit::getSampleID(temp_filt)

      saveRDS(excluded_ids,
              file.path(interim_dir, "excluded_samples.rds"))

      qc_params <- list(
        final_min_cov       = current_cov,
        final_hi_perc       = current_hi,
        final_lo_perc       = current_lo_p,
        kept_chrs           = kept_chrs,
        excluded_samples    = excluded_ids,
        original_sample_ids = original_sample_ids,
        analyzed_sample_ids = analyzed_sample_ids,
        n_iterations        = loop_iter,
        timestamp           = Sys.time()
      )
      saveRDS(qc_params,
              file.path(interim_dir, "qc_final_params.rds"))

      if (!is.null(all_controls_precomputed)) {
        saveRDS(all_controls_precomputed,
                file.path(interim_dir, "control_genomes_all_samples.rds"))
      }

      if (!is.null(df_gm)) {
        saveRDS(df_gm,
                file.path(interim_dir, "global_methylation_per_sample.rds"))
      }

      message("[mETHYLotest] Original samples: ", length(original_sample_ids))
      message("[mETHYLotest] Analysed samples: ", length(analyzed_sample_ids),
              " (", paste(analyzed_sample_ids, collapse = ", "), ")")
      message("[mETHYLotest] Excluded samples: ", length(excluded_ids),
              " (", if (length(excluded_ids) > 0)
                paste(excluded_ids, collapse = ", ") else "none", ")")
      # ============================================════════════════

      # ── Generate QC Report ──
      mETHYLotest.NGS.GenerateQCReport(
        df_meth          = qc_data$df_meth,
        df_controls      = qc_data$df_controls,
        df_controls_all  = all_controls_precomputed,
        df_stats         = qc_data$df_stats,
        df_global_meth   = df_gm,
        min_cov          = current_cov,
        excluded_samples = excluded_ids,
        project_dir      = cfg$project_dir,
        output_file      = file.path(res_dir, "Final_QC_Report.html"))

      filtered.myobj <- temp_filt
      rm(myobj, qc_data); gc()
      break
    }

    if (ui_res$action == "update") {
      current_cov  <- ui_res$min_coverage
      current_hi   <- ui_res$hi_perc
      current_lo_p <- ui_res$lo_perc
      kept_chrs    <- ui_res$chrs_to_keep
      if (!is.null(ui_res$samples_to_exclude))
        excluded_ids <- unique(c(excluded_ids,
                                 ui_res$samples_to_exclude))
      rm(temp_filt, qc_data); gc()
      loop_iter <- loop_iter + 1L
    }
  }
  .end_step()

  # ========================================================================
  # 3b. COVERAGE NORMALIZATION (optional)
  # ========================================================================

  if (isTRUE(cfg$do_normalize_coverage)) {
    message("[mETHYLotest] === Coverage Normalization ===")
    .start_step("Coverage_Normalization")

    norm_method <- if (!is.null(cfg$normalize_cov_method))
      cfg$normalize_cov_method else "median"

    message("[mETHYLotest] Method: ", norm_method)
    message("[mETHYLotest] Samples before: ", length(filtered.myobj))

    # Log coverage stats before
    cov_before <- vapply(filtered.myobj, function(s) {
      median(methylKit::getData(s)$coverage)
    }, numeric(1))
    message("[mETHYLotest] Median coverage per sample (before): ",
            paste(round(cov_before, 1), collapse = ", "))

    filtered.myobj <- methylKit::normalizeCoverage(
      filtered.myobj,
      method = norm_method)

    # Log coverage stats after
    cov_after <- vapply(filtered.myobj, function(s) {
      median(methylKit::getData(s)$coverage)
    }, numeric(1))
    message("[mETHYLotest] Median coverage per sample (after):  ",
            paste(round(cov_after, 1), collapse = ", "))

    .end_step()
  } else {
    message("[mETHYLotest] Coverage normalization skipped.")
  }



  # ========================================================================
  # 4. COMPARATIVE ANALYSIS
  # ========================================================================

  message("[mETHYLotest] === Comparative Analysis ===")
  .start_step("Comparative_Analysis")

  comp_dir <- file.path(res_dir, "Comparative_Analysis")
  fig_dir  <- file.path(comp_dir, "figures")
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  if (length(filtered.myobj) < 2L)
    stop("[mETHYLotest] < 2 samples remaining.")

  destrand_val <- isTRUE(cfg$unite_destrand)
  message("[mETHYLotest] Unite (destrand=", destrand_val, ")...")

  meth <- methylKit::unite(filtered.myobj,
                           destrand = destrand_val,
                           mc.cores = cfg$num_cores)
  message("[mETHYLotest] Common CpGs: ", nrow(meth))
  saveRDS(meth, file.path(interim_dir, "methylBase_united.rds"))

  # Correlation
  dist_method    <- if (!is.null(cfg$cluster_dist))
    cfg$cluster_dist else "correlation"
  cluster_method <- if (!is.null(cfg$cluster_method))
    cfg$cluster_method else "ward"

  png(file.path(fig_dir, "Sample_Correlation.png"),
      width = 1000, height = 1000, res = 150)
  methylKit::getCorrelation(meth, plot = TRUE)
  dev.off()

  pdf(file.path(fig_dir, "Sample_Correlation.pdf"),
      width = 10, height = 10)
  methylKit::getCorrelation(meth, plot = TRUE)
  dev.off()

  # Clustering
  png(file.path(fig_dir, "Sample_Clustering.png"),
      width = 1000, height = 800, res = 150)
  methylKit::clusterSamples(meth, dist = dist_method,
                            method = cluster_method, plot = TRUE)
  dev.off()

  pdf(file.path(fig_dir, "Sample_Clustering.pdf"),
      width = 10, height = 8)
  methylKit::clusterSamples(meth, dist = dist_method,
                            method = cluster_method, plot = TRUE)
  dev.off()

  hc <- methylKit::clusterSamples(meth, dist = dist_method,
                                  method = cluster_method,
                                  plot = FALSE)
  saveRDS(hc, file.path(interim_dir, "clustering_hc_object.rds"))
  
  # PCA
  png(file.path(fig_dir, "QC_PCA.png"), width=800, height=800, res=150)
  methylKit::PCASamples(meth)
  dev.off()
  
  pca_res <- methylKit::PCASamples(meth, obj.return=TRUE)
  if (!is.null(pca_res) && !is.null(pca_res$x)) {
    pca_coords <- data.frame(Sample_Name = rownames(pca_res$x), pca_res$x)
    write.csv(pca_coords, file.path(fig_dir, "PCA_coords.csv"), row.names = FALSE)
  }
  
  # Export to QC dir for web interface
  tryCatch({
    qc_dir <- file.path(res_dir, "QC")
    if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive=TRUE)
    file.copy(file.path(fig_dir, "Sample_Correlation.png"), file.path(qc_dir, "Sample_Correlation.png"), overwrite=TRUE)
    file.copy(file.path(fig_dir, "Sample_Clustering.png"), file.path(qc_dir, "Sample_Clustering.png"), overwrite=TRUE)
    file.copy(file.path(fig_dir, "QC_PCA.png"), file.path(qc_dir, "QC_PCA.png"), overwrite=TRUE)
    file.copy(file.path(fig_dir, "PCA_coords.csv"), file.path(qc_dir, "PCA_coords.csv"), overwrite=TRUE)
    file.copy(file.path(interim_dir, "clustering_hc_object.rds"), file.path(qc_dir, "clustering_hc_object.rds"), overwrite=TRUE)
  }, error = function(e) warning("Failed to copy QC files: ", e$message))

  .end_step()

  # ========================================================================
  # 5a. CpG DIFFERENTIAL METHYLATION
  # ========================================================================

  message("[mETHYLotest] === Differential Methylation ===")

  .start_step("Differential_Methylation")

  pheno_full    <- readxl::read_excel(cfg$pheno_file)
  meth_ids      <- methylKit::getSampleID(meth)
  pheno_matched <- pheno_full[match(meth_ids,
                                    pheno_full[[cfg$col_sampleID]]), ]

  get_cov_df <- function(cols) {
    if (is.null(cols) || length(cols) == 0L) return(NULL)
    df <- as.data.frame(pheno_matched[, cols, drop = FALSE])
    row.names(df) <- meth_ids
    for (c in colnames(df))
      if (is.character(df[[c]])) df[[c]] <- as.factor(df[[c]])
    df
  }

  df_batch <- if (isTRUE(cfg$perform_batch_correction))
    get_cov_df(cfg$batch_cols) else NULL
  df_covs  <- if (isTRUE(cfg$perform_covariate_adj))
    get_cov_df(cfg$covariates) else NULL

  # Build scenarios
  scenarios <- list(
    Basic = list(covariates = NULL,
                 desc = "No correction"))

  if (!is.null(df_batch))
    scenarios$Batch_Correction <- list(
      covariates = df_batch, desc = "Batch only")
  if (!is.null(df_covs))
    scenarios$Covariate_Adj <- list(
      covariates = df_covs, desc = "Covariates only")
  if (!is.null(df_batch) && !is.null(df_covs))
    scenarios$Full_Model <- list(
      covariates = cbind(df_batch, df_covs),
      desc = "Batch + Covariates")

  message("[mETHYLotest] ", length(scenarios), " scenario(s).")

  overdispersion <- if (!is.null(cfg$diff_overdispersion))
    cfg$diff_overdispersion else "MN"
  diff_test      <- if (!is.null(cfg$diff_test))
    cfg$diff_test else "Chisq"
  diff_cutoff    <- if (!is.null(cfg$diff_cutoff))
    cfg$diff_cutoff else 10
  diff_qvalue    <- if (!is.null(cfg$diff_qvalue))
    cfg$diff_qvalue else 0.05

  diff_results <- list()

  # ── Safe cores ──
  safe_cores <- 1L
  if (.Platform$OS.type != "windows") {
    avail <- parallel::detectCores(logical = FALSE)
    if (is.na(avail)) avail <- 1L
    meth_gb <- as.numeric(object.size(meth)) / 1e9
    safe_cores <- min(
      if (!is.null(cfg$num_cores)) cfg$num_cores else 1L,
      avail,
      max(1L, floor(4 / max(meth_gb, 0.1))),
      6L
    )
  }
  message("[mETHYLotest] Using ", safe_cores, " core(s) for DiffMeth")

  # ── Safe calculateDiffMeth wrapper ──
  safe_calcDiffMeth <- function(meth, covs, od, tst, nc) {
    result <- NULL; ok <- TRUE
    tryCatch(
      withCallingHandlers({
        result <- methylKit::calculateDiffMeth(
          meth, covariates = covs, overdispersion = od,
          test = tst, mc.cores = nc)
      }, warning = function(w) {
        if (grepl("scheduled cores encountered errors", w$message)) {
          ok <<- FALSE; invokeRestart("muffleWarning")
        }
      }),
      error = function(e) { ok <<- FALSE }
    )
    if (!ok) return(NULL)
    result
  }

  diff_results <- list()

  for (model in names(scenarios)) {
    message("[mETHYLotest] Running: ", model, " [", scenarios[[model]]$desc, "]")
    covs <- scenarios[[model]]$covariates

    if (!is.null(covs) && min(table(meth@treatment)) <= (ncol(covs) + 1)) {
      warning("[mETHYLotest] Skipping ", model, ": overfitting risk.")
      next
    }

    # Try 1: safe_cores
    res <- safe_calcDiffMeth(meth, covs, overdispersion, diff_test, safe_cores)

    # Try 2: 1 core
    if (is.null(res) && safe_cores > 1L) {
      message("[mETHYLotest]   Retrying with mc.cores=1...")
      res <- safe_calcDiffMeth(meth, covs, overdispersion, diff_test, 1L)
    }

    # Try 3: simplified
    if (is.null(res)) {
      message("[mETHYLotest]   Retrying with overdispersion='none'...")
      res <- safe_calcDiffMeth(meth, covs, "none", "Chisq", 1L)
    }

    # Try 4: without covariates
    if (is.null(res) && !is.null(covs)) {
      message("[mETHYLotest]   Fallback without covariates...")
      res <- safe_calcDiffMeth(meth, NULL, "none", "Chisq", 1L)
    }

    if (!is.null(res)) {
      diff_results[[model]] <- res
      message("[mETHYLotest] ", model, ": success.")
    } else {
      message("[mETHYLotest] ", model, ": ALL ATTEMPTS FAILED.")
    }
  }

  if (length(diff_results) == 0L)
    stop("[mETHYLotest] All scenarios failed.")

  # Export
  results_dir <- file.path(res_dir, "Differential_Analysis")
  if (!dir.exists(results_dir))
    dir.create(results_dir, recursive = TRUE)

  if (length(diff_results) > 0L) {
    saveRDS(diff_results, file.path(interim_dir, "diff_results_list.rds"))

    # ============================================══════════════════
    # Pre-compute per-group methylation means from united object
    # ============================================══════════════════
    perc_meth  <- NULL
    meth_pos   <- NULL
    ctrl_idx   <- NULL
    case_idx   <- NULL

    tryCatch({
      perc_meth <- methylKit::percMethylation(meth)
      meth_data <- methylKit::getData(meth)
      meth_pos  <- paste0(meth_data$chr, ":", meth_data$start)
      ctrl_idx  <- which(meth@treatment == 0)
      case_idx  <- which(meth@treatment == 1)
      message("[mETHYLotest] Pre-computed methylation matrix: ",
              nrow(perc_meth), " CpGs x ", ncol(perc_meth), " samples",
              " (", length(ctrl_idx), " CTL, ", length(case_idx), " Test)")
    }, error = function(e) {
      message("[mETHYLotest] Warning: could not compute per-group means: ", e$message)
    })

    # Helper: enrich a DMC dataframe with group means and status
    enrich_dmc_df <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(df)

      # Status column
      if ("meth.diff" %in% colnames(df)) {
        df$Status <- ifelse(df$meth.diff > 0, "Hyper", "Hypo")
      }

      # Per-group means (requires pre-computed perc_meth)
      if (!is.null(perc_meth) && !is.null(meth_pos) &&
          all(c("chr", "start") %in% colnames(df))) {

        df_pos    <- paste0(df$chr, ":", df$start)
        match_idx <- match(df_pos, meth_pos)

        df$Mean_CTL  <- NA_real_
        df$Mean_Test <- NA_real_
        df$Delta     <- NA_real_

        valid <- !is.na(match_idx)
        if (any(valid) && length(ctrl_idx) > 0) {
          df$Mean_CTL[valid] <- round(
            rowMeans(perc_meth[match_idx[valid], ctrl_idx, drop = FALSE],
                     na.rm = TRUE), 2)
        }
        if (any(valid) && length(case_idx) > 0) {
          df$Mean_Test[valid] <- round(
            rowMeans(perc_meth[match_idx[valid], case_idx, drop = FALSE],
                     na.rm = TRUE), 2)
        }
        if (any(valid)) {
          df$Delta[valid] <- round(df$Mean_Test[valid] - df$Mean_CTL[valid], 2)
        }

        # Reorder: put new columns after meth.diff
        if ("meth.diff" %in% colnames(df)) {
          diff_pos <- which(colnames(df) == "meth.diff")
          new_cols <- c("Status", "Mean_CTL", "Mean_Test", "Delta")
          new_cols <- intersect(new_cols, colnames(df))
          other_cols <- setdiff(colnames(df), new_cols)
          # Insert after meth.diff
          before <- other_cols[seq_len(diff_pos)]
          after  <- if (diff_pos < length(other_cols))
            other_cols[(diff_pos + 1):length(other_cols)]
          else character(0)
          df <- df[, c(before, new_cols, after)]
        }
      }

      df
    }
    # ============================================══════════════════

    any_exported <- FALSE

    for (model in names(diff_results)) {
      safe   <- gsub("[^A-Za-z0-9_.-]", "_", model)
      dm_obj <- diff_results[[model]]

      # Diagnostic: check raw data first
      raw_df <- methylKit::getData(dm_obj)
      raw_df <- raw_df[!is.na(raw_df$qvalue) &
                         !is.na(raw_df$meth.diff), ]

      message("[mETHYLotest] ", model, " raw stats:")
      message("[mETHYLotest]   Total positions: ", nrow(raw_df))

      if (nrow(raw_df) > 0L) {
        message("[mETHYLotest]   Q-value range: [",
                signif(min(raw_df$qvalue, na.rm = TRUE), 3), " - ",
                signif(max(raw_df$qvalue, na.rm = TRUE), 3), "]")
        message("[mETHYLotest]   Meth.diff range: [",
                round(min(raw_df$meth.diff, na.rm = TRUE), 1), "% - ",
                round(max(raw_df$meth.diff, na.rm = TRUE), 1), "%]")
        message("[mETHYLotest]   Positions with q < ", diff_qvalue, ": ",
                sum(raw_df$qvalue < diff_qvalue, na.rm = TRUE))
        message("[mETHYLotest]   Positions with |diff| >= ", diff_cutoff, "%: ",
                sum(abs(raw_df$meth.diff) >= diff_cutoff, na.rm = TRUE))
        message("[mETHYLotest]   Positions passing BOTH: ",
                sum(raw_df$qvalue < diff_qvalue &
                      abs(raw_df$meth.diff) >= diff_cutoff, na.rm = TRUE))
      }

      # Try extraction with configured thresholds
      res_obj <- tryCatch(
        methylKit::getMethylDiff(
          dm_obj,
          difference = diff_cutoff,
          qvalue     = diff_qvalue),
        error = function(e) {
          message("[mETHYLotest]   getMethylDiff failed: ", e$message)
          NULL
        })

      if (is.null(res_obj)) next

      df_res <- methylKit::getData(res_obj)
      message("[mETHYLotest]   DMCs found: ", nrow(df_res))

      if (nrow(df_res) > 0L) {
        any_exported <- TRUE

        # ── Enrich with group means and status ──
        df_res <- enrich_dmc_df(df_res)

        # Excel
        writexl::write_xlsx(
          df_res,
          file.path(results_dir, paste0("DMC_", safe, ".xlsx")))

        # BED
        bed <- data.frame(
          chrom      = df_res$chr,
          chromStart = df_res$start - 1L,
          chromEnd   = df_res$end,
          name       = paste0("DMC_", seq_len(nrow(df_res)),
                              "_", df_res$Status),
          score      = pmin(round(abs(df_res$meth.diff) * 10), 1000L),
          strand     = df_res$strand)
        write.table(
          bed,
          file.path(results_dir, paste0("DMC_", safe, ".bed")),
          quote = FALSE, sep = "\t",
          row.names = FALSE, col.names = FALSE)

        message("[mETHYLotest]   Exported: ", safe,
                " (", nrow(df_res), " DMCs | ",
                sum(df_res$Status == "Hyper"), " Hyper / ",
                sum(df_res$Status == "Hypo"), " Hypo)")

      } else {
        message("[mETHYLotest]   No DMCs at diff>=", diff_cutoff,
                "% & q<", diff_qvalue,
                ". Trying relaxed thresholds...")

        # Relaxed export: just q-value, no diff cutoff
        res_relaxed <- tryCatch(
          methylKit::getMethylDiff(dm_obj, difference = 0,
                                   qvalue = diff_qvalue),
          error = function(e) NULL)

        n_relaxed <- if (!is.null(res_relaxed))
          nrow(methylKit::getData(res_relaxed)) else 0L

        message("[mETHYLotest]   With diff>=0% & q<", diff_qvalue,
                ": ", n_relaxed, " positions")

        if (n_relaxed > 0L) {
          df_relaxed <- methylKit::getData(res_relaxed)
          df_relaxed <- enrich_dmc_df(df_relaxed)
          writexl::write_xlsx(
            df_relaxed,
            file.path(results_dir,
                      paste0("DMC_relaxed_", safe, ".xlsx")))
          message("[mETHYLotest]   Relaxed export saved: DMC_relaxed_",
                  safe, ".xlsx")
          any_exported <- TRUE
        }
      }

      # Always export full results (all positions, no filter)
      raw_df <- enrich_dmc_df(raw_df)
      full_path <- file.path(results_dir,
                             paste0("Full_", safe, ".xlsx"))
      tryCatch(
        writexl::write_xlsx(raw_df, full_path),
        error = function(e) NULL)
        
      # Generate Volcano Plot
      tryCatch({
        if (all(c("meth.diff", "qvalue") %in% colnames(raw_df))) {
          plot_df <- raw_df
          plot_df$status <- "Unchanged"
          plot_df$status[plot_df$meth.diff > 0 & plot_df$qvalue < diff_qvalue] <- "Hyper"
          plot_df$status[plot_df$meth.diff < 0 & plot_df$qvalue < diff_qvalue] <- "Hypo"
          # handle qvalue == 0
          plot_df$qvalue[plot_df$qvalue == 0] <- 1e-300
          plot_df$logQ <- -log10(plot_df$qvalue)
          
          p_volc <- ggplot2::ggplot(plot_df, ggplot2::aes(x = meth.diff, y = logQ, color = status)) +
            ggplot2::geom_point(alpha = 0.6) +
            ggplot2::scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "Unchanged" = "gray")) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste("Volcano Plot:", safe), x = "Diff Meth (%)", y = "-log10(Q-value)")
          ggplot2::ggsave(file.path(results_dir, paste0("Volcano_", safe, ".png")), plot = p_volc, width = 8, height = 6)
        }
      }, error = function(e) warning("[mETHYLotest] Failed to generate Volcano plot: ", e$message))
      
      # Generate Distribution Plot
      tryCatch({
        if ("qvalue" %in% colnames(raw_df)) {
          p_dist <- ggplot2::ggplot(raw_df, ggplot2::aes(x = qvalue)) +
            ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "black") +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste("Q-Value Distribution:", safe), x = "Q-Value", y = "Count")
          ggplot2::ggsave(file.path(results_dir, paste0("Distribution_", safe, ".png")), plot = p_dist, width = 8, height = 6)
        }
      }, error = function(e) warning("[mETHYLotest] Failed to generate Distribution plot: ", e$message))

      message("[mETHYLotest]   Full results: Full_", safe, ".xlsx",
              " (", nrow(raw_df), " positions)")
    }

    if (!any_exported)
      warning("[mETHYLotest] No DMCs passed filters for any scenario.",
              " Check your thresholds (diff_cutoff=", diff_cutoff,
              "%, diff_qvalue=", diff_qvalue, ").",
              " Full unfiltered results have been exported.")

  } else {
    stop("[mETHYLotest] All scenarios failed.")
  }

  .end_step()

  # ========================================================================
  # 5b. TILING WINDOWS (Region-level DMR)
  # ========================================================================

  tiles_results <- NULL

  if (isTRUE(cfg$do_tiling)) {
    message("[mETHYLotest] === Tiling Windows ===")
    .start_step("Tiling_Windows")

    win_size  <- if (!is.null(cfg$tiling_win_size)) cfg$tiling_win_size else 1000L
    step_size <- if (!is.null(cfg$tiling_step_size)) cfg$tiling_step_size else 1000L
    min_cov   <- if (!is.null(cfg$tiling_min_cov)) cfg$tiling_min_cov else 3L

    message("[mETHYLotest] Window: ", win_size, "bp | Step: ",
            step_size, "bp | Min CpGs: ", min_cov)

    tiles_dir <- file.path(res_dir, "Tiling_Windows")
    if (!dir.exists(tiles_dir))
      dir.create(tiles_dir, recursive = TRUE)

    tryCatch({
      # Tile the methylBase object
      tiles <- methylKit::tileMethylCounts(
        meth,
        win.size  = as.integer(win_size),
        step.size = as.integer(step_size),
        cov.bases = as.integer(min_cov))

      message("[mETHYLotest] Tiled regions: ", nrow(tiles))

      # ── Safe cores for tiling (same logic as 5a) ──
      tiles_safe_cores <- 1L
      if (.Platform$OS.type != "windows") {
        avail <- parallel::detectCores(logical = FALSE)
        if (is.na(avail)) avail <- 1L
        tiles_gb <- as.numeric(object.size(tiles)) / 1e9
        tiles_safe_cores <- min(
          if (!is.null(cfg$num_cores)) cfg$num_cores else 1L,
          avail,
          max(1L, floor(4 / max(tiles_gb, 0.1))),
          6L
        )
      }
      message("[mETHYLotest] Using ", tiles_safe_cores,
              " core(s) for Tiling DiffMeth")

      # ── Differential methylation on tiles ──
      tiles_results <- list()

      for (model in names(scenarios)) {
        covs <- scenarios[[model]]$covariates

        if (!is.null(covs) && nrow(covs) <= (2 + ncol(covs))) {
          message("[mETHYLotest] Tiling skip ", model,
                  ": overfitting risk.")
          next
        }

        message("[mETHYLotest] Tiling running: ", model,
                " [", scenarios[[model]]$desc, "]")

        # Try 1: safe_cores
        res <- safe_calcDiffMeth(tiles, covs, overdispersion,
                                 diff_test, tiles_safe_cores)

        # Try 2: 1 core
        if (is.null(res) && tiles_safe_cores > 1L) {
          message("[mETHYLotest]   Tiling retrying with mc.cores=1...")
          res <- safe_calcDiffMeth(tiles, covs, overdispersion,
                                   diff_test, 1L)
        }

        # Try 3: simplified
        if (is.null(res)) {
          message("[mETHYLotest]   Tiling retrying with overdispersion='none'...")
          res <- safe_calcDiffMeth(tiles, covs, "none", "Chisq", 1L)
        }

        # Try 4: without covariates
        if (is.null(res) && !is.null(covs)) {
          message("[mETHYLotest]   Tiling fallback without covariates...")
          res <- safe_calcDiffMeth(tiles, NULL, "none", "Chisq", 1L)
        }

        if (is.null(res)) {
          message("[mETHYLotest] Tiling ", model, ": ALL ATTEMPTS FAILED.")
          next
        }

        dm_tiles <- res
        tiles_results[[model]] <- dm_tiles
        message("[mETHYLotest] Tiling ", model, ": success.")

        # Extract significant
        sig_tiles <- tryCatch(
          methylKit::getMethylDiff(dm_tiles,
                                   difference = diff_cutoff,
                                   qvalue = diff_qvalue),
          error = function(e) NULL)

        n_sig <- if (!is.null(sig_tiles))
          nrow(methylKit::getData(sig_tiles)) else 0L

        message("[mETHYLotest] Tiling ", model,
                ": ", nrow(dm_tiles), " regions, ",
                n_sig, " DMRs")

        safe <- gsub("[^A-Za-z0-9_.-]", "_", model)

        # ════════════════════════════════════════════════════════
        # Pre-build DMP GRanges for cross-referencing
        # ════════════════════════════════════════════════════════
        dmp_gr <- NULL
        if (model %in% names(diff_results)) {
          tryCatch({
            dmp_obj <- methylKit::getMethylDiff(
              diff_results[[model]],
              difference = diff_cutoff,
              qvalue     = diff_qvalue)
            dmp_data <- methylKit::getData(dmp_obj)
            if (nrow(dmp_data) > 0L) {
              dmp_gr <- GenomicRanges::GRanges(
                seqnames = dmp_data$chr,
                ranges   = IRanges::IRanges(
                  start = dmp_data$start,
                  end   = dmp_data$end))
              message("[mETHYLotest]   DMP reference: ",
                      length(dmp_gr), " DMPs for cross-check")
            }
          }, error = function(e)
            message("[mETHYLotest]   Could not build DMP reference: ",
                    e$message))
        }

        # ════════════════════════════════════════════════════════
        # Helper: annotate DMR df with DMP overlap + confidence
        # ════════════════════════════════════════════════════════
        annotate_dmr_confidence <- function(df) {
          if (is.null(df) || nrow(df) == 0L) return(df)

          df$n_DMP      <- 0L
          df$Confidence <- "Unsupported"

          if (!is.null(dmp_gr)) {
            tryCatch({
              dmr_gr <- GenomicRanges::GRanges(
                seqnames = df$chr,
                ranges   = IRanges::IRanges(
                  start = df$start,
                  end   = df$end))
              df$n_DMP <- GenomicRanges::countOverlaps(dmr_gr, dmp_gr)
            }, error = function(e)
              message("[mETHYLotest]   Overlap computation failed: ",
                      e$message))
          }

          df$Confidence <- ifelse(
            df$n_DMP >= 3L & abs(df$meth.diff) >= 25, "High",
            ifelse(
              df$n_DMP >= 1L & abs(df$meth.diff) >= 15, "Medium",
              ifelse(
                df$n_DMP >= 1L, "Low",
                "Unsupported")))

          df
        }

        # ════════════════════════════════════════════════════════
        # Export significant DMRs
        # ════════════════════════════════════════════════════════
        if (n_sig > 0L) {
          df_tiles <- methylKit::getData(sig_tiles)

          # Status Hyper/Hypo
          df_tiles$Status <- ifelse(df_tiles$meth.diff > 0,
                                    "Hyper", "Hypo")

          # DMP cross-reference + confidence
          df_tiles <- annotate_dmr_confidence(df_tiles)

          # Log confidence breakdown
          conf_table <- table(df_tiles$Confidence)
          message("[mETHYLotest]   DMR confidence: ",
                  paste(names(conf_table), conf_table,
                        sep = "=", collapse = " | "))
          message("[mETHYLotest]   DMR with >=1 DMP: ",
                  sum(df_tiles$n_DMP > 0L), "/", nrow(df_tiles),
                  " (", round(100 * sum(df_tiles$n_DMP > 0L) /
                                nrow(df_tiles), 1), "%)")

          # Export ALL DMRs (with confidence)
          writexl::write_xlsx(
            df_tiles,
            file.path(tiles_dir,
                      paste0("DMR_tiles_", safe, "_all.xlsx")))

          # Export supported only (>= 1 DMP)
          df_supported <- df_tiles[df_tiles$Confidence != "Unsupported", ]

          if (nrow(df_supported) > 0L) {
            writexl::write_xlsx(
              df_supported,
              file.path(tiles_dir,
                        paste0("DMR_tiles_", safe, ".xlsx")))

            # BED (supported only)
            bed <- data.frame(
              chrom      = df_supported$chr,
              chromStart = df_supported$start - 1L,
              chromEnd   = df_supported$end,
              name       = paste0("DMR_", seq_len(nrow(df_supported)),
                                  "_", df_supported$Status,
                                  "_", df_supported$Confidence),
              score      = pmin(round(abs(df_supported$meth.diff) * 10),
                                1000L),
              strand     = ".")
            write.table(
              bed,
              file.path(tiles_dir,
                        paste0("DMR_tiles_", safe, ".bed")),
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)

            message("[mETHYLotest]   Exported: ", safe,
                    " (", nrow(df_supported), " supported DMRs from ",
                    nrow(df_tiles), " total | ",
                    sum(df_supported$Status == "Hyper"), " Hyper / ",
                    sum(df_supported$Status == "Hypo"), " Hypo)")
                    
            # Heatmap of top 50 DMRs
            message("[mETHYLotest]   Generating Heatmap (Top 50 DMRs)...")
            tryCatch({
              sig_sorted <- df_tiles[order(df_tiles$qvalue), ]
              top50 <- head(sig_sorted, 50)
              
              tiles_perc <- methylKit::percMethylation(tiles)
              tiles_data <- methylKit::getData(tiles)
              tiles_pos <- paste0(tiles_data$chr, ":", tiles_data$start, "-", tiles_data$end)
              
              top50_pos <- paste0(top50$chr, ":", top50$start, "-", top50$end)
              idx_top <- match(top50_pos, tiles_pos)
              idx_top <- idx_top[!is.na(idx_top)]
              
              if (length(idx_top) > 0) {
                pm_top <- tiles_perc[idx_top, , drop=FALSE]
                rownames(pm_top) <- top50_pos[1:length(idx_top)]
                colnames(pm_top) <- keep_ids
                
                if (requireNamespace("pheatmap", quietly = TRUE)) {
                  annot_col <- data.frame(
                    Group = ifelse(seq_along(keep_ids) %in% case_idx, "Test", "Control"),
                    row.names = keep_ids
                  )
                  
                  pheatmap::pheatmap(
                    pm_top,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_colnames = TRUE,
                    annotation_col = annot_col,
                    show_rownames = FALSE,
                    main = sprintf("Top 50 DMRs Heatmap - %s", safe),
                    filename = file.path(tiles_dir, sprintf("Heatmap_tiles_%s.png", safe)),
                    width = 8, height = 6
                  )
                  write.csv(pm_top, file.path(tiles_dir, sprintf("Heatmap_tiles_data_%s.csv", safe)))
                } else {
                  message("[mETHYLotest]   pheatmap not installed.")
                }
              }
            }, error=function(e) message("[mETHYLotest]   Heatmap generation failed: ", e$message))

          } else {
            message("[mETHYLotest]   No supported DMRs ",
                    "(all ", nrow(df_tiles), " lack DMP evidence)")
          }

        } else {
          message("[mETHYLotest]   No DMRs at diff>=", diff_cutoff,
                  "% & q<", diff_qvalue, " for ", model)
        }

        # ════════════════════════════════════════════════════════
        # Per-group means (computed LAZILY, only for export)
        # ════════════════════════════════════════════════════════
        raw_tiles <- methylKit::getData(dm_tiles)
        raw_tiles <- raw_tiles[!is.na(raw_tiles$qvalue) &
                                 !is.na(raw_tiles$meth.diff), ]

        if (nrow(raw_tiles) > 0L) {
          raw_tiles$Status <- ifelse(raw_tiles$meth.diff > 0,
                                     "Hyper", "Hypo")

          # DMP cross-reference on full results too
          raw_tiles <- annotate_dmr_confidence(raw_tiles)

          # Compute group means only if manageable size
          if (nrow(raw_tiles) <= 500000L) {
            tryCatch({
              tiles_perc <- methylKit::percMethylation(tiles)
              tiles_data <- methylKit::getData(tiles)
              tiles_pos  <- paste0(tiles_data$chr, ":",
                                   tiles_data$start)

              df_pos    <- paste0(raw_tiles$chr, ":",
                                  raw_tiles$start)
              match_idx <- match(df_pos, tiles_pos)
              valid     <- !is.na(match_idx)

              raw_tiles$Mean_CTL  <- NA_real_
              raw_tiles$Mean_Test <- NA_real_
              raw_tiles$Delta     <- NA_real_

              if (any(valid) && length(ctrl_idx) > 0)
                raw_tiles$Mean_CTL[valid] <- round(
                  rowMeans(tiles_perc[match_idx[valid], ctrl_idx,
                                      drop = FALSE],
                           na.rm = TRUE), 2)
              if (any(valid) && length(case_idx) > 0)
                raw_tiles$Mean_Test[valid] <- round(
                  rowMeans(tiles_perc[match_idx[valid], case_idx,
                                      drop = FALSE],
                           na.rm = TRUE), 2)
              if (any(valid))
                raw_tiles$Delta[valid] <- round(
                  raw_tiles$Mean_Test[valid] -
                    raw_tiles$Mean_CTL[valid], 2)

              # Also enrich sig supported DMRs if they exist
              if (n_sig > 0L &&
                  exists("df_supported") &&
                  nrow(df_supported) > 0L) {

                sig_pos   <- paste0(df_supported$chr, ":",
                                    df_supported$start)
                sig_match <- match(sig_pos, tiles_pos)
                sig_valid <- !is.na(sig_match)

                df_supported$Mean_CTL  <- NA_real_
                df_supported$Mean_Test <- NA_real_
                df_supported$Delta     <- NA_real_

                if (any(sig_valid) && length(ctrl_idx) > 0)
                  df_supported$Mean_CTL[sig_valid] <- round(
                    rowMeans(tiles_perc[sig_match[sig_valid],
                                        ctrl_idx, drop = FALSE],
                             na.rm = TRUE), 2)
                if (any(sig_valid) && length(case_idx) > 0)
                  df_supported$Mean_Test[sig_valid] <- round(
                    rowMeans(tiles_perc[sig_match[sig_valid],
                                        case_idx, drop = FALSE],
                             na.rm = TRUE), 2)
                if (any(sig_valid))
                  df_supported$Delta[sig_valid] <- round(
                    df_supported$Mean_Test[sig_valid] -
                      df_supported$Mean_CTL[sig_valid], 2)

                # Re-export enriched supported DMRs
                writexl::write_xlsx(
                  df_supported,
                  file.path(tiles_dir,
                            paste0("DMR_tiles_", safe, ".xlsx")))
              }

              rm(tiles_perc, tiles_data, tiles_pos)
              gc()

              message("[mETHYLotest]   Group means computed for tiles.")

            }, error = function(e) {
              message("[mETHYLotest]   Warning: could not compute tile ",
                      "group means: ", e$message)
            })
          } else {
            message("[mETHYLotest]   Skipping group means: too many ",
                    "regions (", nrow(raw_tiles), " > 500k)")
          }
        }

        # Full export
        tryCatch(
          writexl::write_xlsx(
            raw_tiles,
            file.path(tiles_dir,
                      paste0("Full_tiles_", safe, ".xlsx"))),
          error = function(e) NULL)
        message("[mETHYLotest]   Full results: Full_tiles_", safe,
                ".xlsx (", nrow(raw_tiles), " regions)")
      }

      saveRDS(tiles_results,
              file.path(interim_dir, "tiles_diff_results.rds"))

    }, error = function(e)
      warning("[mETHYLotest] Tiling failed: ", e$message))
    .end_step()

  } else {
    message("[mETHYLotest] Tiling windows skipped.")
  }

  # ========================================================================
  # 5c. SEGMENTATION (methSeg)
  # ========================================================================

  seg_results <- NULL

  if (isTRUE(cfg$do_segmentation)) {
    message("[mETHYLotest] === Methylation Segmentation ===")
    .start_step("Segmentation")

    seg_dir <- file.path(res_dir, "Segmentation")
    if (!dir.exists(seg_dir))
      dir.create(seg_dir, recursive = TRUE)

    seg_G      <- if (!is.null(cfg$seg_G)) cfg$seg_G else 4L
    seg_min    <- if (!is.null(cfg$seg_min_seg)) cfg$seg_min_seg else 500L

    tryCatch({
      # methSeg works on a single methylDiff or methylRaw object
      # We use the first diff result available
      first_model <- names(diff_results)[1L]
      if (!is.null(first_model)) {

        message("[mETHYLotest] Segmentation on: ", first_model,
                " | G=", seg_G)

        seg_results <- methylKit::methSeg(
          diff_results[[first_model]],
          diagnostic.plot = FALSE,
          G = as.integer(seg_G))

        message("[mETHYLotest] Segments found: ", length(seg_results))

        # Convert to data.frame for export
        seg_df <- as.data.frame(seg_results)

        # Filter by minimum size
        if ("width" %in% colnames(seg_df)) {
          seg_df <- seg_df[seg_df$width >= seg_min, ]
          message("[mETHYLotest] After min size filter (",
                  seg_min, "bp): ", nrow(seg_df), " segments")
        }

        # Export
        if (nrow(seg_df) > 0L) {
          writexl::write_xlsx(
            seg_df, file.path(seg_dir, "segments.xlsx"))

          # BED
          bed <- data.frame(
            chrom      = seg_df$seqnames,
            chromStart = seg_df$start - 1L,
            chromEnd   = seg_df$end,
            name       = paste0("Seg_", seq_len(nrow(seg_df)),
                                "_state", seg_df$seg.group),
            score      = round(abs(seg_df$seg.mean) * 10),
            strand     = ".")
          write.table(
            bed, file.path(seg_dir, "segments.bed"),
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
        }

        saveRDS(seg_results, file.path(interim_dir, "methSeg_results.rds"))

        # Summary plot
        tryCatch({
          grDevices::pdf(file.path(seg_dir, "segment_summary.pdf"),
                         width = 12, height = 6)
          methylKit::methSeg(diff_results[[first_model]],
                             diagnostic.plot = TRUE,
                             G = as.integer(seg_G))
          grDevices::dev.off()
          message("[mETHYLotest] Segmentation diagnostic plot saved.")
        }, error = function(e) NULL)
        
        # Genome Track plot with ggplot2
        tryCatch({
          if (requireNamespace("ggplot2", quietly = TRUE)) {
            message("[mETHYLotest]   Generating Segmentation Track Plot (Genome Track)...")
            p <- ggplot2::ggplot(seg_df, ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = seg.mean, fill = seg.mean)) +
              ggplot2::geom_rect() +
              ggplot2::facet_wrap(~ seqnames, scales = "free_x") +
              ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50, limits=c(0, 100)) +
              ggplot2::theme_minimal() +
              ggplot2::labs(title = paste("Segmentation Track -", first_model), x = "Position", y = "Methylation Mean (%)", fill="Meth (%)") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
            
            ggplot2::ggsave(file.path(seg_dir, paste0("Segmentation_Track_", first_model, ".png")), plot = p, width = 12, height = 8, dpi = 300)
            message("[mETHYLotest]   Segmentation Track saved.")
          }
        }, error = function(e) message("[mETHYLotest]   Segmentation Track plot failed: ", e$message))
      }

    }, error = function(e)
      warning("[mETHYLotest] Segmentation failed: ", e$message))

    .end_step()

  } else {
    message("[mETHYLotest] Segmentation skipped.")
  }

  # ========================================================================
  # 6. SIGNATURES (default + user)
  # ========================================================================

  sig_dir <- file.path(res_dir, "Signatures")
  if (!dir.exists(sig_dir)) dir.create(sig_dir, recursive = TRUE)

  # ── A. Generate default signatures from DiffMeth results ───────────────
  max_sig_sites <- 500L

  for (model in names(diff_results)) {
    dm_obj <- diff_results[[model]]
    raw_df <- methylKit::getData(dm_obj)
    raw_df <- raw_df[!is.na(raw_df$qvalue) &
                       !is.na(raw_df$meth.diff), ]

    # Filter significant
    sig_df <- raw_df[raw_df$qvalue < diff_qvalue &
                       abs(raw_df$meth.diff) >= diff_cutoff, ]

    if (nrow(sig_df) == 0L) {
      message("[mETHYLotest] No significant DMCs for '", model,
              "'. Default signature skipped.")
      next
    }

    # Score = -log10(qvalue) * |meth.diff|
    sig_df$score <- -log10(pmax(sig_df$qvalue, 1e-300)) *
      abs(sig_df$meth.diff)
    sig_df <- sig_df[order(sig_df$score, decreasing = TRUE), ]

    if (nrow(sig_df) > max_sig_sites) {
      message("[mETHYLotest] '", model, "': ", nrow(sig_df),
              " DMCs. Selecting top ", max_sig_sites, ".")
      sig_df <- sig_df[seq_len(max_sig_sites), ]
    }

    # Write signature as chr:start-end (unique genomic IDs)
    sig_ids <- paste0(sig_df$chr, ":", sig_df$start, "-", sig_df$end)

    safe <- gsub("[^A-Za-z0-9_.-]", "_", model)
    sig_path <- file.path(sig_dir,
                          paste0("DMC_default_", safe, ".txt"))
    writeLines(sig_ids, sig_path)

    message("[mETHYLotest] Default signature: ", basename(sig_path),
            " | ", length(sig_ids), " sites",
            " | q: [",
            signif(min(sig_df$qvalue), 2), " - ",
            signif(max(sig_df$qvalue), 2),
            "] | |diff|: [",
            round(min(abs(sig_df$meth.diff)), 1), "% - ",
            round(max(abs(sig_df$meth.diff)), 1), "%]")
  }

  # ── B. Results UI (user can explore + create signatures) ───────────────
  tryCatch({
    mETHYLotest.NGS.DiffMeth.UI(
      diff_results,
      output_dir = file.path(res_dir, "DiffMeth_Figures"))
  }, error = function(e) {
    warning("[mETHYLotest] DiffMeth UI failed: ", e$message)
  })

  # ── C. Inventory all signatures ────────────────────────────────────────
  all_sigs  <- list.files(sig_dir, pattern = "\\.txt$",
                          full.names = TRUE)
  n_default <- sum(grepl("^DMC_default_", basename(all_sigs)))
  n_user    <- length(all_sigs) - n_default

  message("[mETHYLotest] Signatures: ",
          n_default, " default + ", n_user, " user = ",
          length(all_sigs), " total")

  # ── D. Signature Validation ────────────────────────────────────────────
  validation_results <- NULL

  if (length(all_sigs) > 0L) {
    message("[mETHYLotest] === Signature Validation ===")

    tryCatch({
      validation_results <- mETHYLotest.NGS.validate(
        meth              = meth,
        diff_results      = diff_results,
        signatures_folder = sig_dir,
        diff_cutoff       = diff_cutoff,
        diff_qvalue       = diff_qvalue
      )

      if (isTRUE(cfg$save_raw_obj)) {
        saveRDS(validation_results,
                file.path(interim_dir, "validation_results.rds"))
      }

      if (length(validation_results) > 0L) {
        # Shiny interactive UI
        message("[mETHYLotest] Launching Signature Validator UI...")
        tryCatch({
          mETHYLotest.NGS.Validate.UI(
            validation_results = validation_results)
        }, error = function(e) {
          warning("[mETHYLotest] Validator UI failed: ", e$message)
        })

        # Static HTML report
        message("[mETHYLotest] Generating validation report...")
        tryCatch({
          mETHYLotest.NGS.Validate.Report(
            validation_results = validation_results,
            output_dir         = file.path(res_dir, "Validation"),
            filename           = "NGS_Validation_Report.html")
        }, error = function(e) {
          warning("[mETHYLotest] Validation report failed: ", e$message)
        })
      } else {
        message("[mETHYLotest] No signatures passed validation.")
      }
    }, error = function(e) {
      warning("[mETHYLotest] Validation failed: ", e$message)
    })
  } else {
    message("[mETHYLotest] No signature files found. Skipping validation.")
  }

  # ========================================================================
  # 7. ANNOTATION
  # ========================================================================

  message("[mETHYLotest] === Annotation ===")

  .start_step("Annotation")

  annot_diff <- if (!is.null(cfg$annot_diff_cutoff))
    cfg$annot_diff_cutoff else 25
  annot_qval <- if (!is.null(cfg$annot_qval_cutoff))
    cfg$annot_qval_cutoff else 0.05

  annotated_data <- tryCatch({
    mETHYLotest.NGS.AnnotateDMCs(
      diff_obj    = diff_results,
      assembly    = cfg$assembly,
      output_dir  = file.path(res_dir, "Annotation"),
      diff_cutoff = annot_diff,
      qval_cutoff = annot_qval)
  }, error = function(e) {
    warning("[mETHYLotest] Annotation failed: ", e$message,
            "\n  Tip: install TxDb package with BiocManager::install('TxDb.Hsapiens.UCSC.",
            cfg$assembly, ".knownGene')")
    NULL
  })

  .end_step()

  # ========================================================================
  # 8. REPORT
  # ========================================================================

  message("[mETHYLotest] === Final Report ===")

  .start_step("Final_report")

  # ── Determine final sample lists ──
  final_sample_ids <- methylKit::getSampleID(meth)
  all_original_ids <- as.character(
    as.data.frame(readxl::read_excel(cfg$pheno_file))[[cfg$col_sampleID]])
  final_excluded <- setdiff(all_original_ids, final_sample_ids)

  tryCatch({
    mETHYLotest.NGS.GenerateReport(
      project_dir      = project_directory,
      output_file      = "Final_Report.html",
      analyzed_samples = final_sample_ids,
      excluded_samples = final_excluded
    )
  }, error = function(e)
    warning("[mETHYLotest] Report failed: ", e$message))

  .end_step()

  # ========================================================================
  # PROFILING SUMMARY
  # ========================================================================
  total_elapsed <- (proc.time() - .timer$start_total)[["elapsed"]]

  profile_df <- do.call(rbind, lapply(names(.timer$steps), function(nm) {
    s <- .timer$steps[[nm]]
    data.frame(
      Step        = nm,
      Time_sec    = s$elapsed_sec,
      Time_human  = sprintf("%dm %02ds", s$elapsed_sec %/% 60,
                            round(s$elapsed_sec) %% 60),
      Peak_RAM_MB = s$peak_ram_mb,
      stringsAsFactors = FALSE
    )
  }))
  profile_df$Pct_Total <- round(100 * profile_df$Time_sec / total_elapsed, 1)

  # Disk footprint
  disk_mb <- tryCatch({
    files <- list.files(cfg$project_dir, recursive = TRUE, full.names = TRUE)
    round(sum(file.info(files)$size, na.rm = TRUE) / 1e6, 1)
  }, error = function(e) NA)

  message("[mETHYLotest] ============================================")
  message("[mETHYLotest] PERFORMANCE PROFILE")
  message("[mETHYLotest] ============================================")
  for (i in seq_len(nrow(profile_df))) {
    message(sprintf("[mETHYLotest]   %-30s %8s  (%4.1f%%)  RAM: %6.0f MB",
                    profile_df$Step[i],
                    profile_df$Time_human[i],
                    profile_df$Pct_Total[i],
                    profile_df$Peak_RAM_MB[i]))
  }
  message("[mETHYLotest] ============================================")
  message(sprintf("[mETHYLotest]   TOTAL: %dm %02ds",
                  total_elapsed %/% 60, round(total_elapsed) %% 60))
  if (!is.na(disk_mb))
    message(sprintf("[mETHYLotest]   Disk footprint: %.1f MB", disk_mb))
  message("[mETHYLotest] ============================================")

  # Save profile
  profile_df$Total_sec <- total_elapsed
  profile_df$Disk_MB   <- disk_mb
  saveRDS(profile_df, file.path(rds_dir, "pipeline_profile.rds"))
  writexl::write_xlsx(profile_df,
                      file.path(res_dir, "Pipeline_Performance.xlsx"))

  # Final messages and return

  message("[mETHYLotest] ========================================")
  message("[mETHYLotest] Pipeline complete!")
  message("[mETHYLotest] Project: ", cfg$project_name)
  message("[mETHYLotest] Results: ", res_dir)
  message("[mETHYLotest] ========================================")

  invisible(list(
    meth              = meth,
    diff_results      = diff_results,
    tiles_results     = tiles_results,
    seg_results       = seg_results,
    validation_results = validation_results,
    annotated         = annotated_data,
    config            = cfg))
}

