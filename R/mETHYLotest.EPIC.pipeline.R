#' mETHYLotest Pipeline for EPIC / 450K Methylation Array (ChAMP wrapper)
#'
#' @description
#' Runs the full mETHYLotest EPIC analysis pipeline in order:
#' \enumerate{
#'   \item \code{champ.load} — Import IDAT files (with smart caching)
#'   \item Pre-normalisation QC
#'   \item \code{champ.impute} — Imputation
#'   \item \code{champ.norm} — Normalisation
#'   \item Post-normalisation QC
#'   \item \code{champ.refbase} — Cell type correction (blood only)
#'   \item \code{champ.runCombat} — Batch correction
#'   \item Post-ComBat QC
#'   \item \code{champ.DMP} — Differentially methylated positions
#'   \item Signature generation + validation (SVM/PCA/Silhouette)
#'   \item DMP Explorer UI + Validation UI
#'   \item DMP HTML reports
#'   \item \code{champ.DMR} — Differentially methylated regions
#'   \item \code{champ.Block} — Block methylation analysis
#'   \item \code{champ.GSEA} — Gene set enrichment analysis
#'   \item \code{champ.CNA} — Copy number aberration analysis
#'   \item Final HTML report + performance profiling
#' }
#'
#' Each step operates on the output of the previous step via a single
#' \code{beta_current} matrix that flows through the entire pipeline.
#'
#' Passing \code{"demo"} as \code{project_directory} launches a
#' self-contained demonstration using 450K IDATs from ChAMPdata.
#'
#' @param project_directory Character. Path to the project directory.
#'   If empty or invalid, launches the Project Setup UI.
#'   If \code{"demo"}, launches demo mode with ChAMPdata IDATs.
#'
#' @return A named list (invisibly) containing:
#'   \describe{
#'     \item{myLoad}{ChAMP load object (beta, pd, M, intensity...)}
#'     \item{beta_final}{Final beta matrix after all corrections}
#'     \item{beta_history}{Character string describing all transformations applied}
#'     \item{pheno}{Phenotype vector used for analysis}
#'     \item{pheno_col}{Name of the phenotype column used}
#'     \item{arraytype}{Detected array type (450K, EPICv1, EPICv2, mouse)}
#'     \item{myDMP}{DMP results (list of data.frames, or NULL)}
#'     \item{myDMR}{DMR results (or NULL)}
#'     \item{myBlock}{Block results (or NULL)}
#'     \item{myGSEA}{GSEA results (or NULL)}
#'     \item{myCNA}{CNA results (or NULL)}
#'     \item{myRefBase}{Cell type correction results (or NULL)}
#'     \item{config}{Pipeline configuration list}
#'   }
#'
#' @details
#' The pipeline includes performance profiling. Wall-clock time, peak RAM,
#' and disk footprint are recorded per step and exported as
#' \code{Pipeline_Performance.xlsx} in the Results directory.
#'
#' Intermediate objects are cached as \code{.rds} files. On re-runs,
#' \code{myLoad.rds} is detected and import is skipped entirely.
#'
#' @export
#' @importFrom ChAMP champ.load champ.impute champ.norm champ.refbase
#'   champ.runCombat champ.DMP champ.DMR champ.Block champ.GSEA champ.CNA
#' @import ChAMPdata
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx
#' @importFrom parallel detectCores
#' @importFrom grDevices pdf dev.off hcl.colors colorRampPalette
#' @importFrom graphics barplot legend par
#' @importFrom tools file_ext
#' @importFrom utils object.size write.csv
mETHYLotest.EPIC.pipeline <- function(project_directory = "") {

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


  # - DEMO MODE - #
  if (identical(tolower(trimws(project_directory)), "demo")) {
    message("[mETHYLotest] ========================================")
    message("[mETHYLotest]        EPIC DEMO MODE")
    message("[mETHYLotest] ========================================")

    if (!requireNamespace("ChAMPdata", quietly = TRUE))
      stop("[mETHYLotest] Package 'ChAMPdata' required for demo.\n",
           "  Install with: BiocManager::install('ChAMPdata')")

    champ_dir <- system.file("extdata", package = "ChAMPdata")
    if (champ_dir == "")
      stop("[mETHYLotest] ChAMPdata extdata not found.")

    demo_xlsx <- system.file("extdata", "EPIC_demo",
                             "demo_sample_sheet.xlsx",
                             package = "mETHYLotest")
    if (demo_xlsx == "")
      stop("[mETHYLotest] demo_sample_sheet.xlsx not found.")

    # ── Copy IDATs + sample sheet to a writable temp dir ──
    demo_dir <- file.path(tempdir(), "mETHYLotest_EPIC_demo")
    if (dir.exists(demo_dir)) unlink(demo_dir, recursive = TRUE)
    dir.create(demo_dir, recursive = TRUE)

    # Copy IDATs
    idats <- list.files(champ_dir, pattern = "\\.idat$", full.names = TRUE)
    file.copy(idats, demo_dir)
    message("[mETHYLotest] Copied ", length(idats), " IDAT files to temp dir")

    # Copy sample sheet as xlsx (for our UI validation)
    file.copy(demo_xlsx, file.path(demo_dir, "demo_sample_sheet.xlsx"))

    demo_pheno <- file.path(demo_dir, "demo_sample_sheet.xlsx")

    message("[mETHYLotest] Demo directory: ", demo_dir)
    message("[mETHYLotest] Sample sheet:   ", demo_pheno)

    project_directory <- mETHYLotest.EPIC.ProjectUI(
      prefill_pheno    = demo_pheno,
      prefill_idat_dir = demo_dir)
  }

  # ========================================================================
  # 1. PROJECT SETUP
  # ========================================================================

  if (!nzchar(project_directory) || !dir.exists(project_directory))
    project_directory <- mETHYLotest.EPIC.ProjectUI()

  config_path <- normalizePath(
    file.path(project_directory, "Results", "project_config.R"))
  if (!file.exists(config_path))
    stop("[mETHYLotest] Config not found: ", config_path)

  source(file = config_path)
  mETHYLotest.EPIC.Check_ProjectSettings(config_path)

  cfg <- project_config

  # ── Directories ──────────────────────────────────────────────────────────
  rds_dir       <- cfg$rds_dir
  res_dir       <- cfg$res_dir
  qc_dir        <- file.path(res_dir, "QC_Raw")
  qc_norm_dir   <- file.path(res_dir, "QC_Normalised")
  qc_combat_dir <- file.path(res_dir, "QC_Combat")
  celltype_dir  <- file.path(res_dir, "CellType")
  dmp_dir       <- file.path(res_dir, "DMP")
  dmr_dir       <- file.path(res_dir, "DMR")
  block_dir     <- file.path(res_dir, "Block")
  gsea_dir      <- file.path(res_dir, "GSEA")
  cna_dir       <- file.path(res_dir, "CNA")

  for (d in c(rds_dir, qc_dir))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # ========================================================================
  # 2. PHENOTYPE IMPORT
  # ========================================================================
  .start_step("Config_and_Import")

  if (!file.exists(cfg$pheno_file))
    stop("[mETHYLotest] Phenotype file not found: ", cfg$pheno_file)

  Pheno <- if (tools::file_ext(cfg$pheno_file) %in% c("xlsx", "xls")) {
    as.data.frame(readxl::read_excel(cfg$pheno_file))
  } else {
    read.csv(cfg$pheno_file, header = TRUE, stringsAsFactors = FALSE)
  }
  message("[mETHYLotest] Phenotype loaded: ", nrow(Pheno), " samples.")

  col_name  <- cfg$col_sample_name
  col_group <- cfg$col_sample_group
  col_plate <- cfg$col_sample_plate
  col_sid   <- cfg$col_sentrix_id
  col_pos   <- cfg$col_sentrix_pos

  message("[mETHYLotest] Groups: ",
          paste(unique(Pheno[[col_group]]), collapse = ", "))

  # ========================================================================
  # 3. ARRAY TYPE DETECTION
  # ========================================================================

  plate_to_champ <- function(plate) {
    switch(toupper(trimws(plate)),
           "450K"   = "450K",
           "EPICV1" = "EPICv1",
           "EPICV2" = "EPICv2",
           "MOUSE" = "mouse",
           stop("[mETHYLotest] Unknown array: ", plate))
  }

  plate_values <- unique(Pheno[[col_plate]])
  champ_types  <- vapply(plate_values, plate_to_champ, character(1))
  message("[mETHYLotest] Arrays: ",
          paste(plate_values, "->", champ_types, collapse = " | "))

  idat_dir      <- cfg$idat_dir
  pheno_csv_tmp <- file.path(idat_dir, "Pheno.csv")

  # ── ChAMP-compatible CSV writer ──────────────────────────────────────────
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
    write.table(out, file = pheno_csv_tmp,
                row.names = FALSE, quote = FALSE, sep = ",")
  }

  # ── champ.load base arguments ────────────────────────────────────────────
  champ_load_base <- list(
    directory      = idat_dir,
    method         = if (!is.null(cfg$load_method)) cfg$load_method
    else "ChAMP",
    autoimpute     = isTRUE(cfg$load_autoimpute),
    filterDetP     = isTRUE(cfg$filter_det_p),
    detPcut        = cfg$det_p_cut,
    SampleCutoff   = if (!is.null(cfg$sample_det_p_cut))
      cfg$sample_det_p_cut else 0.1,
    filterBeads    = isTRUE(cfg$filter_beads),
    beadCutoff     = cfg$bead_cutoff,
    filterNoCG     = isTRUE(cfg$filter_no_cg),
    filterSNPs     = isTRUE(cfg$filter_snps),
    population     = cfg$load_population,
    filterMultiHit = isTRUE(cfg$filter_multi_hit),
    filterXY       = isTRUE(cfg$filter_xy),
    force          = isTRUE(cfg$load_force)
  )

  # ========================================================================
  # 4. LOAD DATA BY ARRAY TYPE (Smart Load)
  # ========================================================================
  myLoad_rds <- file.path(rds_dir, "myLoad.rds")

  if (file.exists(myLoad_rds)) {
    message("[mETHYLotest] Loading existing myLoad object: ", myLoad_rds)
    myLoad <- readRDS(myLoad_rds)

    # Recover arraytype from cached object or config
    if (!is.null(cfg$technology) && !is.na(cfg$technology)) {
      arraytype <- cfg$technology
    } else {
      # Infer from plate values
      primary <- intersect(c("EPICv1", "EPICv2", "450K", "Mouse"),
                           plate_values)[1L]
      arraytype <- plate_to_champ(primary)
    }

    message("[mETHYLotest] Loaded from cache | CpGs: ",
            format(nrow(myLoad$beta), big.mark = ","),
            " | Samples: ", ncol(myLoad$beta))

  } else {
    message("[mETHYLotest] Importing from IDAT files...")

    loads_list <- list()

    for (plate_val in plate_values) {
      champ_at <- plate_to_champ(plate_val)
      message("[mETHYLotest] Loading ", plate_val, " (", champ_at, ") ...")

      pheno_sub <- Pheno[Pheno[[col_plate]] == plate_val, , drop = FALSE]
      write_champ_pheno(pheno_sub)

      loads_list[[plate_val]] <- do.call(
        ChAMP::champ.load,
        c(champ_load_base, list(arraytype = champ_at)))

      message("[mETHYLotest] ", plate_val,
              " | CpGs: ",    nrow(loads_list[[plate_val]]$beta),
              " | Samples: ", ncol(loads_list[[plate_val]]$beta))
    }

    # ====================================================================
    # 5. HARMONIZE
    # ====================================================================

    if (length(loads_list) == 1L) {
      myLoad    <- loads_list[[1L]]
      arraytype <- plate_to_champ(names(loads_list))
    } else {
      myLoad  <- mETHYLotest.utils.HarmonizeArrays(loads_list)
      primary <- intersect(c("EPICv1", "EPICv2", "450K"),
                           names(loads_list))[1L]
      arraytype <- plate_to_champ(primary)
    }

    message("[mETHYLotest] Loaded | CpGs: ",
            format(nrow(myLoad$beta), big.mark = ","),
            " | Samples: ", ncol(myLoad$beta))

    # Save for next run
    if (isTRUE(cfg$save_raw_obj)) {
      if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
      saveRDS(myLoad, myLoad_rds)
    message("[mETHYLotest] Saved: ", myLoad_rds)
    }
  }

  rownames(myLoad$pd) <- myLoad$pd[["Sample_Name"]]
  
  # Save the phenotype table to a CSV for the Web App to read
  utils::write.csv(myLoad$pd, file.path(res_dir, "Samples_Phenotype.csv"), row.names = FALSE)

  # Check if biological variable exists in Pheno
  bio_var <- cfg$biological_variable
  
  # --- Unified Technology Detection ---
  cfg$technology <- arraytype
  message("[mETHYLotest] Detected technology: ", cfg$technology)

  # ── Helper: remove samples ───────────────────────────────────────────────
  .remove_samples <- function(load_obj, to_remove, beta_mat = NULL) {
    if (length(to_remove) == 0L)
      return(list(load = load_obj, beta = beta_mat))

    message("[mETHYLotest] Removing ", length(to_remove),
            " sample(s): ", paste(to_remove, collapse = ", "))

    keep <- !colnames(load_obj$beta) %in% to_remove
    load_obj$pd   <- load_obj$pd[!rownames(load_obj$pd) %in% to_remove,
                                 , drop = FALSE]
    load_obj$beta <- load_obj$beta[, keep, drop = FALSE]

    for (slot in c("M", "intensity", "detP", "beadcount"))
      if (!is.null(load_obj[[slot]]))
        load_obj[[slot]] <- load_obj[[slot]][, keep, drop = FALSE]

    tryCatch({
      ph <- Pheno[!Pheno[[col_name]] %in% to_remove, , drop = FALSE]
      writexl::write_xlsx(ph, cfg$pheno_file)
    }, error = function(e)
      warning("[mETHYLotest] Could not update pheno: ", e$message))

    out <- list(load = load_obj)
    if (!is.null(beta_mat))
      out$beta <- beta_mat[, !colnames(beta_mat) %in% to_remove,
                           drop = FALSE]
    out
  }

  #if (isTRUE(cfg$save_raw_obj))
  #  saveRDS(myLoad, file.path(rds_dir, "myLoad.rds"))

  if (file.exists(pheno_csv_tmp)) file.remove(pheno_csv_tmp)

  # ========================================================================
  #
  #   PIPELINE FLOW : beta_current tracks the active matrix
  #
  #   champ.load → beta_current (raw)
  #        ↓ QC
  #   champ.impute → beta_current (imputed)
  #        ↓
  #   champ.norm → beta_current (normalised)
  #        ↓ QC
  #   champ.refbase → beta_current (cell-type corrected)
  #        ↓
  #   champ.runCombat → beta_current (batch corrected)
  #        ↓ QC
  #   champ.DMP / champ.DMR / champ.Block / champ.GSEA / champ.CNA
  #
  # ========================================================================

  beta_current <- myLoad$beta
  beta_label   <- "raw"

  .log_beta <- function(label) {
    message("[mETHYLotest] beta_current is now: ", label,
            " | ", format(nrow(beta_current), big.mark = ","),
            " probes x ", ncol(beta_current), " samples")
  }

  .log_beta(beta_label)

  .end_step()

  # ========================================================================
  # 6. PRE-NORMALISATION QC
  # ========================================================================

  message("[mETHYLotest] ═══ Step 1: Pre-normalisation QC ═══")

  .start_step("QC")

  samples_rem <- mETHYLotest.EPIC.QC(
    myLoad = myLoad, outputDir = qc_dir)

  if (length(samples_rem) > 0L) {
    res          <- .remove_samples(myLoad, samples_rem, beta_current)
    myLoad       <- res$load
    beta_current <- res$beta
    .log_beta(beta_label)
  }

  .end_step()

  # ========================================================================
  # 7. IMPUTATION (champ.impute)
  # ========================================================================

  if (isTRUE(cfg$do_impute)) {
    message("[mETHYLotest] ═══ Step 2: Imputation (champ.impute) ═══")

    .start_step("Imputation")

    message("[mETHYLotest] Input: beta_current [", beta_label, "]")

    impute_result <- ChAMP::champ.impute(
      beta         = beta_current,
      pd           = myLoad$pd,
      method       = if (!is.null(cfg$impute_method)) cfg$impute_method
      else "Combine",
      k            = if (!is.null(cfg$impute_k)) cfg$impute_k else 5L,
      ProbeCutoff  = if (!is.null(cfg$impute_probe_cutoff))
        cfg$impute_probe_cutoff else 0.2,
      SampleCutoff = if (!is.null(cfg$impute_sample_cutoff))
        cfg$impute_sample_cutoff else 0.1
    )

    removed_probes  <- setdiff(rownames(beta_current),
                               rownames(impute_result$beta))
    removed_samples <- setdiff(colnames(beta_current),
                               colnames(impute_result$beta))

    if (length(removed_probes) > 0L)
      message("[mETHYLotest] Imputation removed ",
              length(removed_probes), " probes.")
    if (length(removed_samples) > 0L)
      message("[mETHYLotest] Imputation removed ",
              length(removed_samples), " samples.")

    beta_current <- impute_result$beta
    myLoad$pd    <- impute_result$pd
    beta_label   <- "imputed"

    # Sync myLoad slots
    kept_p <- rownames(beta_current)
    kept_s <- colnames(beta_current)
    myLoad$beta <- beta_current
    for (slot in c("M", "intensity", "detP", "beadcount")) {
      if (!is.null(myLoad[[slot]])) {
        cr <- intersect(rownames(myLoad[[slot]]), kept_p)
        cc <- intersect(colnames(myLoad[[slot]]), kept_s)
        myLoad[[slot]] <- myLoad[[slot]][cr, cc, drop = FALSE]
      }
    }

    .log_beta(beta_label)

    if (isTRUE(cfg$save_raw_obj))
      saveRDS(impute_result, file.path(rds_dir, "myImpute.rds"))

    .end_step()

  } else {
    message("[mETHYLotest] Imputation skipped.")
  }

  # ========================================================================
  # 8. NORMALISATION (champ.norm)
  # ========================================================================

  if (isTRUE(cfg$normalize_data)) {
    message("[mETHYLotest] ═══ Step 3: Normalisation (",
            cfg$norm_method, ") ═══")

    .start_step("Normalization")

    message("[mETHYLotest] Input: beta_current [", beta_label, "]")

    if (!dir.exists(qc_norm_dir))
      dir.create(qc_norm_dir, recursive = TRUE)

    norm_time <- system.time({
      myNorm <- ChAMP::champ.norm(
        beta       = beta_current,
        rgSet      = myLoad$rgSet,
        mset       = myLoad$mset,
        resultsDir = qc_norm_dir,
        method     = cfg$norm_method,
        plotBMIQ   = isTRUE(cfg$plot_norm),
        arraytype  = arraytype,
        cores      = if (!is.null(cfg$norm_cores)) cfg$norm_cores
        else cfg$num_cores
      )
    })

    beta_current <- myNorm
    beta_label   <- paste0("normalised (", cfg$norm_method, ")")

    message("[mETHYLotest] Normalisation: ",
            round(norm_time[["elapsed"]], 1), " s")
    .log_beta(beta_label)

    if (isTRUE(cfg$save_raw_obj))
      saveRDS(myNorm, file.path(rds_dir, "myNorm.rds"))

    # Post-normalisation QC
    message("[mETHYLotest] ═══ Step 4: Post-normalisation QC ═══")
    samples_rem_n <- mETHYLotest.EPIC.QC(
      myLoad = myLoad, beta_matrix = beta_current,
      outputDir = qc_norm_dir)

    if (length(samples_rem_n) > 0L) {
      res          <- .remove_samples(myLoad, samples_rem_n, beta_current)
      myLoad       <- res$load
      beta_current <- res$beta
      .log_beta(beta_label)
    }

    .end_step()

  } else {
    message("[mETHYLotest] Normalisation skipped.")
  }

  # ========================================================================
  # 9. CELL TYPE CORRECTION (champ.refbase) — blood only
  # ========================================================================

  myRefBase <- NULL

  if (isTRUE(cfg$do_refbase)) {
    message("[mETHYLotest] ═══ Step 5: Cell type correction ═══")

    .start_step("Celltype_correction")

    message("[mETHYLotest] Input: beta_current [", beta_label, "]")

    beta_for_ref <- beta_current

    # Clamp boundaries
    n_zero <- sum(beta_for_ref <= 0, na.rm = TRUE)
    n_one  <- sum(beta_for_ref >= 1, na.rm = TRUE)
    if (n_zero > 0L || n_one > 0L) {
      message("[mETHYLotest] Clamping: ", n_zero, " vals <= 0, ",
              n_one, " vals >= 1")
      beta_for_ref[beta_for_ref <= 0] <- 0.001
      beta_for_ref[beta_for_ref >= 1] <- 0.999
    }

    # Remove NA probes
    na_probes <- rowSums(is.na(beta_for_ref)) > 0L
    if (any(na_probes)) {
      message("[mETHYLotest] Removing ", sum(na_probes), " NA probes.")
      beta_for_ref <- beta_for_ref[!na_probes, , drop = FALSE]
    }

    # Remove near-zero variance probes
    pvar <- apply(beta_for_ref, 1, var, na.rm = TRUE)
    if (any(pvar < 1e-10)) {
      message("[mETHYLotest] Removing ", sum(pvar < 1e-10),
              " near-zero variance probes.")
      beta_for_ref <- beta_for_ref[pvar >= 1e-10, , drop = FALSE]
    }

    myRefBase <- tryCatch({
      ChAMP::champ.refbase(beta = beta_for_ref, arraytype = arraytype)
    }, error = function(e) {
      warning("[mETHYLotest] refbase failed: ", e$message,
              ". Retrying with [0.01, 0.99]...")
      b2 <- beta_for_ref
      b2[b2 < 0.01] <- 0.01
      b2[b2 > 0.99] <- 0.99
      tryCatch(
        ChAMP::champ.refbase(beta = b2, arraytype = arraytype),
        error = function(e2) {
          warning("[mETHYLotest] refbase fallback failed: ", e2$message)
          NULL
        })
    })

    if (!is.null(myRefBase)) {
      corrected     <- myRefBase$CorrectedBeta
      shared_probes <- intersect(rownames(beta_current),
                                 rownames(corrected))
      beta_current  <- corrected[shared_probes, , drop = FALSE]
      beta_label    <- paste0(beta_label, " + refbase")

      .log_beta(beta_label)

      # Save cell fractions
      cell_frac <- myRefBase$CellFraction
      for (ct in colnames(cell_frac))
        message("[mETHYLotest]   ", ct,
                ": mean=", round(mean(cell_frac[, ct]), 3),
                " sd=",   round(sd(cell_frac[, ct]), 3))

      if (isTRUE(cfg$refbase_save_plots)) {
        if (!dir.exists(celltype_dir))
          dir.create(celltype_dir, recursive = TRUE)
        write.csv(cell_frac,
                  file.path(celltype_dir, "cell_fractions.csv"),
                  row.names = TRUE)
        tryCatch({
          cols <- grDevices::hcl.colors(ncol(cell_frac), "Dynamic")
          grDevices::pdf(
            file.path(celltype_dir, "cell_proportions.pdf"),
            width = 12, height = 7)
          graphics::par(mar = c(8, 4, 3, 8), xpd = TRUE)
          graphics::barplot(
            t(cell_frac), beside = FALSE, col = cols, las = 2,
            cex.names = 0.6, border = NA,
            main = "Cell Type Proportions", ylab = "Proportion")
          graphics::legend("topright", inset = c(-0.15, 0),
                           legend = colnames(cell_frac),
                           fill = cols, cex = 0.7, bty = "n")
          grDevices::dev.off()
        }, error = function(e)
          warning("[mETHYLotest] Plot error: ", e$message))
      }

      if (isTRUE(cfg$save_raw_obj))
        saveRDS(myRefBase, file.path(rds_dir, "myRefBase.rds"))

    } else {
      message("[mETHYLotest] Cell type correction failed, continuing.")
    }

    .end_step()

  } else {
    message("[mETHYLotest] Cell type correction skipped.")
  }

  # ========================================================================
  # 10. BATCH CORRECTION (champ.runCombat)
  # ========================================================================

  if (isTRUE(cfg$perform_batch_correction)) {

    message("[mETHYLotest] ═══ Step : Batch Effect correction ═══")

    .start_step("Batch_Effect_correction")

    batch_cols <- cfg$batch_cols
    bio_var    <- if (!is.null(cfg$biological_variable))
      cfg$biological_variable else col_group

    if (is.null(batch_cols) || length(batch_cols) == 0L) {
      warning("[mETHYLotest] Batch correction enabled but no batch_cols.")

    } else {

      # -- Map UI names to pd names --
      col_name_map <- c(
        "Sentrix_ID"       = "Slide",
        "Sentrix_Position" = "Array",
        "Sample_Name"      = "Sample_Name",
        "Sample_Group"     = "Sample_Group"
      )

      map_col <- function(col) {
        mapped <- col_name_map[col]
        ifelse(is.na(mapped), col, mapped)
      }

      batch_cols_pd <- map_col(batch_cols)
      bio_var_pd    <- map_col(bio_var)

      message("[mETHYLotest] ComBat column mapping:")
      for (i in seq_along(batch_cols))
        message("[mETHYLotest]   ", batch_cols[i], " -> ", batch_cols_pd[i])
      message("[mETHYLotest]   protect: ", bio_var, " -> ", bio_var_pd)

      # -- Confounding check (returns clean + confounded) --
      conf_result <- mETHYLotest.EPIC.checkConfounding(
        pd           = myLoad$pd,
        variablename = bio_var_pd,
        batchname    = batch_cols_pd
      )

      clean_batch_pd <- conf_result$clean
      confounded_pd  <- conf_result$confounded

      # Reverse-map for UI display
      reverse_map <- setNames(names(col_name_map), col_name_map)
      reverse_col <- function(col) {
        rev <- reverse_map[col]
        ifelse(is.na(rev), col, rev)
      }

      if (length(confounded_pd) > 0L) {
        confounded_ui <- vapply(confounded_pd, reverse_col, character(1))
        warning(
          "\n",
          paste0(strrep("=", 60), "\n"),
          "[mETHYLotest] CONFOUNDING WARNING\n",
          paste0(strrep("-", 60), "\n"),
          "  Biological variable: ", bio_var, " (", bio_var_pd, ")\n",
          "  EXCLUDED (confounded): ",
          paste(confounded_ui, collapse = ", "), "\n",
          if (length(clean_batch_pd) > 0L)
            paste0("  KEPT (proceeding): ",
                   paste(vapply(clean_batch_pd, reverse_col,
                                character(1)),
                         collapse = ", "), "\n")
          else
            "  No unconfounded variables remain. ComBat SKIPPED.\n",
          paste0(strrep("=", 60), "\n")
        )
      }

      if (length(clean_batch_pd) > 0L) {
        message("[mETHYLotest] === ComBat (",
                paste(clean_batch_pd, collapse = " + "), ") ===")
        message("[mETHYLotest] Input: beta_current [", beta_label, "]")

        if (!dir.exists(qc_combat_dir))
          dir.create(qc_combat_dir, recursive = TRUE)

        for (i in seq_along(clean_batch_pd)) {
          bv    <- clean_batch_pd[i]
          bv_ui <- reverse_col(bv)

          message("[mETHYLotest] ComBat pass ", i, "/",
                  length(clean_batch_pd), ": ", bv_ui, " (", bv, ")")

          tryCatch({
            beta_current <- ChAMP::champ.runCombat(
              beta         = beta_current,
              pd           = myLoad$pd,
              variablename = bio_var_pd,
              batchname    = bv,
              logitTrans   = isTRUE(cfg$combat_logit_transform)
            )

            beta_label <- paste0(beta_label, " + ComBat(", bv_ui, ")")
            .log_beta(beta_label)

          }, error = function(e)
            warning("[mETHYLotest] ComBat failed for '", bv_ui,
                    "': ", e$message, ". Continuing."))
        }

        if (isTRUE(cfg$save_raw_obj))
          saveRDS(beta_current, file.path(rds_dir, "myCombat.rds"))

        # Post-ComBat QC
        message("[mETHYLotest] === Post-ComBat QC ===")
        samples_rem_cb <- mETHYLotest.EPIC.QC(
          myLoad = myLoad, beta_matrix = beta_current,
          outputDir = qc_combat_dir)

        if (length(samples_rem_cb) > 0L) {
          res          <- .remove_samples(myLoad, samples_rem_cb,
                                          beta_current)
          myLoad       <- res$load
          beta_current <- res$beta
          .log_beta(beta_label)
        }

      } else {
        warning("[mETHYLotest] ALL batch variables confounded.",
                " ComBat entirely skipped.")
      }
    }
    .end_step()
  } else {
    message("[mETHYLotest] Batch correction skipped.")
  }

  # ========================================================================
  # 11. FINAL BETA & ANALYSIS VARIABLE
  # ========================================================================

  beta_final <- beta_current

  message("[mETHYLotest] ═══ Final beta matrix ═══")
  message("[mETHYLotest] History: ", beta_label)
  message("[mETHYLotest] Dimensions: ",
          format(nrow(beta_final), big.mark = ","),
          " probes x ", ncol(beta_final), " samples")

  pheno_col <- if (!is.null(cfg$analysis_pheno)) cfg$analysis_pheno else col_group

  if (!pheno_col %in% colnames(myLoad$pd))
    stop("[mETHYLotest] Analysis variable '", pheno_col,
         "' not in phenotype data.")

  pheno_vec   <- myLoad$pd[[pheno_col]]
  compare_grp <- cfg$compare_group

  message("[mETHYLotest] Pheno: '", pheno_col, "' | Levels: ",
          paste(unique(as.character(pheno_vec)), collapse = ", "))

  # ========================================================================
  # 12. DMP + Signatures + Validation + Reports
  # ========================================================================

  myDMP              <- NULL
  validation_results <- NULL
  sig_dir            <- file.path(cfg$res_dir, "Episignatures")

  if (isTRUE(cfg$do_dmp)) {
    message("[mETHYLotest] ═══ DMP (champ.DMP) ═══")
    .start_step("DMP")
    message("[mETHYLotest] Input: beta_final [", beta_label, "]")

    if (!dir.exists(dmp_dir)) dir.create(dmp_dir, recursive = TRUE)

    myDMP <- ChAMP::champ.DMP(
      beta          = beta_final,
      pheno         = pheno_vec,
      compare.group = compare_grp,
      adjPVal       = if (!is.null(cfg$dmp_adj_p_val))
        cfg$dmp_adj_p_val else 0.05,
      adjust.method = if (!is.null(cfg$dmp_adjust_method))
        cfg$dmp_adjust_method else "BH",
      arraytype     = cfg$technology
    )

    if (isTRUE(cfg$save_raw_obj))
      saveRDS(myDMP, file.path(rds_dir, "myDMP.rds"))

    for (comp in names(myDMP)) {
      safe <- gsub("[^A-Za-z0-9_.-]", "_", comp)
      write.csv(myDMP[[comp]],
                file.path(dmp_dir, paste0("DMP_", safe, ".csv")),
                row.names = TRUE)
                
      # Generate Volcano Plot
      tryCatch({
        df_volc <- myDMP[[comp]]
        if (all(c("logFC", "adj.P.Val") %in% colnames(df_volc))) {
          df_volc$status <- "Unchanged"
          df_volc$status[df_volc$logFC > 0 & df_volc$adj.P.Val < 0.05] <- "Hyper"
          df_volc$status[df_volc$logFC < 0 & df_volc$adj.P.Val < 0.05] <- "Hypo"
          df_volc$logQ <- -log10(df_volc$adj.P.Val)
          
          p_volc <- ggplot2::ggplot(df_volc, ggplot2::aes(x = logFC, y = logQ, color = status)) +
            ggplot2::geom_point(alpha = 0.6) +
            ggplot2::scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "Unchanged" = "gray")) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste("Volcano Plot:", comp), x = "logFC", y = "-log10(adj.P.Val)")
          ggplot2::ggsave(file.path(dmp_dir, paste0("Volcano_", safe, ".png")), plot = p_volc, width = 8, height = 6)
        }
      }, error = function(e) warning("[mETHYLotest] Failed to generate Volcano plot: ", e$message))
      
      # Generate Distribution Plot
      tryCatch({
        df_dist <- myDMP[[comp]]
        if ("P.Value" %in% colnames(df_dist)) {
          p_dist <- ggplot2::ggplot(df_dist, ggplot2::aes(x = P.Value)) +
            ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "black") +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste("P-Value Distribution:", comp), x = "P-Value", y = "Count")
          ggplot2::ggsave(file.path(dmp_dir, paste0("Distribution_", safe, ".png")), plot = p_dist, width = 8, height = 6)
        }
      }, error = function(e) warning("[mETHYLotest] Failed to generate Distribution plot: ", e$message))
    }
    message("[mETHYLotest] DMP: ", length(myDMP), " comparison(s).")

    # ── A. Generate default signatures from DMP results ──────────────────
    if (!dir.exists(sig_dir)) dir.create(sig_dir, recursive = TRUE)

    adj_threshold  <- if (!is.null(cfg$dmp_adj_p_val))
      cfg$dmp_adj_p_val else 0.05
    max_sig_cpgs   <- 500L
    min_delta_beta <- 0.1

    for (comp in names(myDMP)) {
      df_comp <- myDMP[[comp]]

      if (!"deltaBeta" %in% colnames(df_comp)) {
        if ("logFC" %in% colnames(df_comp))
          df_comp$deltaBeta <- df_comp$logFC
        else next
      }

      df_sig <- df_comp[df_comp$adj.P.Val < adj_threshold, , drop = FALSE]

      if (nrow(df_sig) == 0L) {
        message("[mETHYLotest] No significant CpGs for '", comp,
                "'. Default signature skipped.")
        next
      }

      df_sig <- df_sig[abs(df_sig$deltaBeta) >= min_delta_beta, ,
                       drop = FALSE]

      if (nrow(df_sig) == 0L) {
        message("[mETHYLotest] No CpGs passing |dB| >= ", min_delta_beta,
                " for '", comp, "'. Relaxing to P-value only.")
        df_sig <- df_comp[df_comp$adj.P.Val < adj_threshold, ,
                          drop = FALSE]
      }

      df_sig$score <- -log10(pmax(df_sig$adj.P.Val, 1e-300)) *
        abs(df_sig$deltaBeta)
      df_sig <- df_sig[order(df_sig$score, decreasing = TRUE), ]

      if (nrow(df_sig) > max_sig_cpgs) {
        message("[mETHYLotest] '", comp, "': ", nrow(df_sig),
                " CpGs passed. Selecting top ", max_sig_cpgs, ".")
        df_sig <- df_sig[seq_len(max_sig_cpgs), ]
      }

      safe     <- gsub("[^A-Za-z0-9_.-]", "_", comp)
      sig_path <- file.path(sig_dir,
                            paste0("DMP_default_", safe, ".txt"))
      writeLines(rownames(df_sig), sig_path)

      message("[mETHYLotest] Default signature: ", basename(sig_path),
              " | ", nrow(df_sig), " CpGs",
              " | adj.P: [",
              formatC(min(df_sig$adj.P.Val), format = "e", digits = 1),
              " - ",
              formatC(max(df_sig$adj.P.Val), format = "e", digits = 1),
              "] | |dB|: [",
              round(min(abs(df_sig$deltaBeta)), 3), " - ",
              round(max(abs(df_sig$deltaBeta)), 3), "]")
    }

    # ── B. DMP Explorer UI (user creates additional signatures) ──────────
    if (!is.null(myDMP) && length(myDMP) > 0L) {
      message("[mETHYLotest] Launching DMP Explorer UI...")
      mETHYLotest.EPIC.DMP.UI(
        dmp_results          = myDMP,
        beta                 = beta_final,
        pheno                = pheno_vec,
        arraytype            = cfg$technology,
        path_to_episignature = sig_dir
      )
    }

    # ── C. Inventory all signatures (default + user) ─────────────────────
    all_sigs  <- list.files(sig_dir, pattern = "\\.txt$",
                            full.names = TRUE)
    n_default <- sum(grepl("^DMP_default_", basename(all_sigs)))
    n_user    <- length(all_sigs) - n_default

    message("[mETHYLotest] Signatures: ",
            n_default, " default + ", n_user, " user = ",
            length(all_sigs), " total")

    # ── D. Signature Validation + UI + Report ────────────────────────────
    if (length(all_sigs) > 0L) {

      message("[mETHYLotest] ═══ Signature Validation ═══")

      tryCatch({
        validation_results <- mETHYLotest.EPIC.validate(
          beta              = beta_final,
          pheno             = pheno_vec,
          signatures_folder = sig_dir,
          arraytype         = cfg$technology,
          k_folds           = 5
        )

        message("[mETHYLotest] Validation complete: ",
                length(validation_results), " signature(s).")

        if (isTRUE(cfg$save_raw_obj))
          saveRDS(validation_results,
                  file.path(rds_dir, "validation_results.rds"))

        if (length(validation_results) > 0L) {
          message("[mETHYLotest] Launching Signature Validator UI...")
          mETHYLotest.EPIC.Validate.UI(
            validation_results = validation_results
          )

          message("[mETHYLotest] Generating validation report...")
          tryCatch({
            mETHYLotest.EPIC.Validate.Report(
              validation_results = validation_results,
              output_dir         = file.path(cfg$res_dir, "Validation"),
              filename           = "Validation_Benchmark_Report.html"
            )
          }, error = function(e)
            warning("[mETHYLotest] Validation report failed: ", e$message))
        }

      }, error = function(e)
        warning("[mETHYLotest] Validation failed: ", e$message))

    } else {
      message("[mETHYLotest] No signatures found. Validation skipped.")
    }

    # ── E. DMP HTML Reports (standard + signatures) ──────────────────────
    message("[mETHYLotest] ═══ Generating DMP Reports ═══")
    tryCatch({
      mETHYLotest.EPIC.DMP.report(
        dmp_results      = myDMP,
        beta             = beta_final,
        pheno            = pheno_vec,
        output_directory = file.path(cfg$res_dir, "DMP_Reports"),
        signatures_folder = if (length(all_sigs) > 0L) sig_dir
        else NULL,
        signature_comparison_source = names(myDMP)[1L],
        adjPVal          = adj_threshold,
        arraytype        = cfg$technology,
        top_n_heatmap    = 50,
        top_n_windows    = 20
      )
    }, error = function(e)
      warning("[mETHYLotest] DMP report failed: ", e$message))
    .end_step()
  } else {
    message("[mETHYLotest] DMP skipped.")
  }

  mETHYLotest.EPIC.Report(
    project_dir        = cfg$project_dir,
    beta               = beta_final,
    pheno              = pheno_vec,
    cfg                = cfg,
    dmp_list           = myDMP,
    signatures_dir     = sig_dir,
    validation_results = validation_results,
    top_n              = 50
  )

  # ========================================================================
  # 13. DMR (champ.DMR)
  # ========================================================================

  myDMR <- NULL

  if (isTRUE(cfg$do_dmr)) {
    dmr_method <- if (!is.null(cfg$dmr_method)) cfg$dmr_method
    else "Bumphunter"
    message("[mETHYLotest] ═══ DMR (", dmr_method, ") ═══")
    message("[mETHYLotest] Input: beta_final [", beta_label, "]")

    if (!dir.exists(dmr_dir)) dir.create(dmr_dir, recursive = TRUE)

    dmr_args <- list(
      beta          = beta_final,
      pheno         = pheno_vec,
      compare.group = compare_grp,
      arraytype     = cfg$technology,
      method        = dmr_method,
      minProbes     = if (!is.null(cfg$dmr_min_probes))
        cfg$dmr_min_probes else 7L,
      adjPvalDmr    = if (!is.null(cfg$dmr_adj_p_val))
        cfg$dmr_adj_p_val else 0.05,
      cores         = if (!is.null(cfg$dmr_cores)) cfg$dmr_cores
      else cfg$num_cores
    )

    if (dmr_method == "Bumphunter") {
      dmr_args$maxGap     <- if (!is.null(cfg$dmr_bh_max_gap))
        cfg$dmr_bh_max_gap else 300L
      dmr_args$cutoff     <- if (!is.null(cfg$dmr_bh_cutoff))
        cfg$dmr_bh_cutoff else NULL
      dmr_args$pickCutoff <- isTRUE(if (!is.null(cfg$dmr_bh_pick_cutoff))
        cfg$dmr_bh_pick_cutoff else TRUE)
      dmr_args$smooth     <- isTRUE(if (!is.null(cfg$dmr_bh_smooth))
        cfg$dmr_bh_smooth else TRUE)
      dmr_args$B          <- if (!is.null(cfg$dmr_bh_B))
        cfg$dmr_bh_B else 250L
    } else if (dmr_method == "DMRcate") {
      dmr_args$lambda <- if (!is.null(cfg$dmr_cate_lambda))
        cfg$dmr_cate_lambda else 1000
      dmr_args$C      <- if (!is.null(cfg$dmr_cate_c))
        cfg$dmr_cate_c else 2L
    } else if (dmr_method == "ProbeLasso") {
      dmr_args$minDmrSep    <- if (!is.null(cfg$dmr_pl_min_sep))
        cfg$dmr_pl_min_sep else 1000L
      dmr_args$minDmrSize   <- if (!is.null(cfg$dmr_pl_min_size))
        cfg$dmr_pl_min_size else 50L
      dmr_args$adjPvalProbe <- if (!is.null(cfg$dmr_pl_adj_p_probe))
        cfg$dmr_pl_adj_p_probe else 0.05
    }

    myDMR <- do.call(ChAMP::champ.DMR, dmr_args)

    # ------------------------------------------------------------------------
    # Sauvegarde des résultats DMR
    # ------------------------------------------------------------------------
    if (!is.null(myDMR)) {
      message("[mETHYLotest] DMR result structure: ", class(myDMR))

      # --- 1. Sauvegarde brute (objet complet) ---
      if (isTRUE(cfg$save_raw_obj)) {
        rds_path <- file.path(rds_dir, "myDMR.rds")
        saveRDS(myDMR, rds_path)
        message("[mETHYLotest] Raw DMR object saved: ", rds_path)
      }

      # --- 2. Dossier DMR Results ---
      if (!dir.exists(dmr_dir)) dir.create(dmr_dir, recursive = TRUE)

      # Convertit le résultat en DataFrame si c’est une liste imbriquée
      if (is.list(myDMR) && !is.data.frame(myDMR)) {
        # ChAMP::champ.DMR renvoie une liste de data.frames : un par comparaison
        for (nm in names(myDMR)) {
          df <- myDMR[[nm]]
          if (!is.data.frame(df)) next

          safe <- gsub("[^A-Za-z0-9_.-]", "_", nm)
          csv_path <- file.path(dmr_dir, paste0("DMR_", safe, ".csv"))
          rds_path <- file.path(dmr_dir, paste0("DMR_", safe, ".rds"))

          # Écriture CSV
          tryCatch({
            utils::write.csv(df, csv_path, row.names = TRUE)
            saveRDS(df, rds_path)
            message("[mETHYLotest] Saved: ", basename(csv_path),
                    " (", nrow(df), " regions)")
          }, error = function(e) {
            warning("[mETHYLotest] Could not save DMR ", nm, ": ", e$message)
          })
        }

        # Sauvegarde globale (liste complète)
        saveRDS(myDMR, file.path(dmr_dir, "DMR_All_Results.rds"))

      } else if (is.data.frame(myDMR)) {
        # Single data.frame returned (rare cas)
        csv_path <- file.path(dmr_dir, "DMR_Results.csv")
        utils::write.csv(myDMR, csv_path, row.names = TRUE)
        saveRDS(myDMR, file.path(dmr_dir, "DMR_Results.rds"))
        message("[mETHYLotest] Single DMR table saved: ", csv_path)
      }

      message("[mETHYLotest] DMR saved successfully in: ", dmr_dir)
    }

    message("[mETHYLotest] DMR complete.")
  } else {
    message("[mETHYLotest] DMR skipped.")
  }

  # ========================================================================
  # 14. BLOCK (champ.Block)
  # ========================================================================

  myBlock <- NULL

  if (isTRUE(cfg$do_block)) {
    message("[mETHYLotest] ═══ Block (champ.Block) ═══")
    message("[mETHYLotest] Input: beta_final [", beta_label, "]")

    if (!dir.exists(block_dir)) dir.create(block_dir, recursive = TRUE)

    # Prep args
    block_args <- list(
      beta      = beta_final,
      pheno     = pheno_vec,
      arraytype = cfg$technology,
      maxClusterGap = if (!is.null(cfg$block_bp_span))
        cfg$block_bp_span else 250000L,
      B         = if (!is.null(cfg$block_B))
        cfg$block_B else 500L,
      bpSpan    = if (!is.null(cfg$block_bp_span))
        cfg$block_bp_span else 250000L,
      minNum    = if (!is.null(cfg$block_min_num))
        cfg$block_min_num else 5L,
      cores     = if (!is.null(cfg$block_cores))
        cfg$block_cores else cfg$num_cores
    )

    # --- tryCatch if errors ---
    myBlock <- tryCatch({
      do.call(ChAMP::champ.Block, block_args)
    }, error = function(e) {
      message("[mETHYLotest] [ERROR] champ.Block a échoué : ", e$message)
      warning("Analyse Block sautée en raison d'une erreur.")
      return(NULL)
    })

    # --- Saving results ---
    if (!is.null(myBlock)) {
      message("[mETHYLotest] Block result structure: ", class(myBlock))

      # 1. Save rds
      if (isTRUE(cfg$save_raw_obj)) {
        rds_interim <- file.path(rds_dir, "myBlock.rds")
        saveRDS(myBlock, rds_interim)
        message("[mETHYLotest] Raw Block object saved: ", rds_interim)
      }

      # 2. Export
      if (is.list(myBlock) && !is.data.frame(myBlock)) {

        for (nm in names(myBlock)) {
          df <- myBlock[[nm]]
          if (!is.data.frame(df)) next

          # Cleaning name file
          safe_name <- gsub("[^A-Za-z0-9_.-]", "_", nm)
          csv_path <- file.path(block_dir, paste0("Block_", safe_name, ".csv"))
          rds_path <- file.path(block_dir, paste0("Block_", safe_name, ".rds"))

          tryCatch({
            utils::write.csv(df, csv_path, row.names = TRUE)
            saveRDS(df, rds_path)
            message("[mETHYLotest] Saved: ", basename(csv_path), " (", nrow(df), " blocks)")
          }, error = function(e) {
            message("[mETHYLotest] Warning: Could not save Block CSV for ", nm, ": ", e$message)
          })
        }

        # Save the results
        saveRDS(myBlock, file.path(block_dir, "Block_All_Results.rds"))

      } else if (is.data.frame(myBlock)) {
        # If only one df
        csv_path <- file.path(block_dir, "Block_Results.csv")
        utils::write.csv(myBlock, csv_path, row.names = TRUE)
        saveRDS(myBlock, file.path(block_dir, "Block_Results.rds"))
        message("[mETHYLotest] Single Block table saved: ", csv_path)
      }

      message("[mETHYLotest] Block saved successfully in: ", block_dir)
    } else {
      message("[mETHYLotest] No Block results generated (Analysis failed or returned NULL).")
    }

    message("[mETHYLotest] Block process complete.")

  } else {
    message("[mETHYLotest] Block skipped.")
  }

  # ========================================================================
  # 15. GSEA (champ.GSEA)
  # ========================================================================

  myGSEA <- NULL

  if (isTRUE(cfg$do_gsea)) {
    # Safety check: GSEA requires DMP results to function
    if (is.null(myDMP) || length(myDMP) == 0) {
      message("[mETHYLotest] [WARNING] GSEA requires a valid DMP object. Step skipped.")
    } else {
      message("[mETHYLotest] ═══ GSEA (champ.GSEA) ═══")
      message("[mETHYLotest] Input: beta_final [", beta_label, "]")

      # Create results directory if it doesn't exist
      if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)

      # Prepare function arguments
      # Note: We use myDMP[[1L]] because champ.DMP usually returns a list of results
      gsea_args <- list(
        beta      = beta_final,
        DMP       = myDMP[[1L]],
        DMR       = myDMR,        # Can be NULL, champ.GSEA handles it
        arraytype = cfg$technology,
        adjPval   = if (!is.null(cfg$gsea_adj_p_val)) cfg$gsea_adj_p_val else 0.05,
        method    = if (!is.null(cfg$gsea_method)) cfg$gsea_method else "fisher"
      )

      # --- Execution with tryCatch safety ---
      myGSEA <- tryCatch({
        do.call(ChAMP::champ.GSEA, gsea_args)
      }, error = function(e) {
        message("[mETHYLotest] [ERROR] champ.GSEA failed: ", e$message)

        # Specific check for the common strsplit error
        if (grepl("strsplit", e$message)) {
          message("[mETHYLotest] Tip: This error often means the Gene columns in your DMP are empty or malformed.")
        }
        return(NULL) # Return NULL so the pipeline can continue to the next step
      })

      # --- Results Saving Logic ---
      if (!is.null(myGSEA)) {
        message("[mETHYLotest] GSEA result structure detected: ", class(myGSEA))

        # 1. Save raw object (Interim/Debug)
        if (isTRUE(cfg$save_raw_obj)) {
          rds_interim <- file.path(rds_dir, "myGSEA.rds")
          saveRDS(myGSEA, rds_interim)
          message("[mETHYLotest] Raw GSEA object saved: ", rds_interim)
        }

        # 2. Export individual result tables (DMP enrichment, DMR enrichment, etc.)
        if (is.list(myGSEA)) {
          for (nm in names(myGSEA)) {
            df_res <- myGSEA[[nm]]

            # Only save if the element is a data frame or a matrix
            if (is.data.frame(df_res) || is.matrix(df_res)) {
              # Clean filename
              safe_name <- gsub("[^A-Za-z0-9_.-]", "_", nm)
              csv_path  <- file.path(gsea_dir, paste0("GSEA_", safe_name, ".csv"))
              rds_path  <- file.path(gsea_dir, paste0("GSEA_", safe_name, ".rds"))

              tryCatch({
                utils::write.csv(df_res, csv_path, row.names = TRUE)
                saveRDS(df_res, rds_path)
                message("[mETHYLotest] Saved GSEA table: ", nm, " (", nrow(df_res), " pathways)")
              }, error = function(e) {
                message("[mETHYLotest] Warning: Could not save GSEA table ", nm, ": ", e$message)
              })
            }
          }

          # Save the full combined list object in the results directory
          saveRDS(myGSEA, file.path(gsea_dir, "GSEA_All_Results.rds"))
        }

        message("[mETHYLotest] GSEA results successfully saved in: ", gsea_dir)
      } else {
        message("[mETHYLotest] No GSEA results were generated due to previous errors.")
      }

      message("[mETHYLotest] GSEA process complete.")
    }
  } else {
    message("[mETHYLotest] GSEA skipped.")
  }

  # ========================================================================
  # 16. CNA (champ.CNA)
  # ========================================================================

  myCNA <- NULL

  if (isTRUE(cfg$do_cna)) {
    # 1. Basic Requirement Check
    if (is.null(myLoad$intensity)) {
      message("[mETHYLotest] [WARNING] CNA requires intensity data (minfi method). Step skipped.")
    } else {
      message("[mETHYLotest] ═══ CNA (champ.CNA) ═══")

      # 2. Data Validation
      ctrl_name <- if (!is.null(cfg$cna_control_group)) cfg$cna_control_group else "CTL"

      # Check if the control group exists in your phenotype vector
      if (!(ctrl_name %in% pheno_vec)) {
        message("[mETHYLotest] [ERROR] Control group '", ctrl_name, "' not found in pheno_vec.")
        message("[mETHYLotest] Available groups: ", paste(unique(pheno_vec), collapse = ", "))
        message("[mETHYLotest] CNA skipped to prevent crash.")
      } else {

        if (!dir.exists(cna_dir)) dir.create(cna_dir, recursive = TRUE)

        # 3. Prepare Arguments based on your provided list
        # Note: 'genome' and 'cores' are removed to avoid "unused argument" errors
        cna_args <- list(
          intensity      = myLoad$intensity,
          pheno          = pheno_vec,
          control        = TRUE,  # Explicitly set to calculate differences
          controlGroup   = ctrl_name,
          sampleCNA      = isTRUE(if (!is.null(cfg$cna_sample_cna)) cfg$cna_sample_cna else TRUE),
          groupFreqPlots = isTRUE(if (!is.null(cfg$cna_group_freq_plots)) cfg$cna_group_freq_plots else TRUE),
          freqThreshold  = if (!is.null(cfg$cna_freq_threshold)) cfg$cna_freq_threshold else 0.3,
          PDFplot        = TRUE,
          Rplot          = FALSE, # Set to FALSE if running on a server without X11
          resultsDir     = paste0(cna_dir, "/"),
          arraytype      = cfg$technology
        )

        # 4. Execution with tryCatch safety
        myCNA <- tryCatch({
          do.call(ChAMP::champ.CNA, cna_args)
        }, error = function(e) {
          message("[mETHYLotest] [ERROR] champ.CNA failed: ", e$message)
          return(NULL)
        })

        # 5. Results Saving Logic
        if (!is.null(myCNA)) {
          message("[mETHYLotest] CNA results detected. Saving data...")

          # Save raw object for R
          if (isTRUE(cfg$save_raw_obj)) {
            saveRDS(myCNA, file.path(rds_dir, "myCNA.rds"))
          }

          # Export list elements (DataFrames) to CSV
          if (is.list(myCNA)) {
            for (nm in names(myCNA)) {
              df_res <- myCNA[[nm]]
              if (is.data.frame(df_res) || is.matrix(df_res)) {
                safe_name <- gsub("[^A-Za-z0-9_.-]", "_", nm)
                csv_path  <- file.path(cna_dir, paste0("CNA_", safe_name, ".csv"))
                utils::write.csv(df_res, csv_path, row.names = TRUE)
              }
            }
            saveRDS(myCNA, file.path(cna_dir, "CNA_All_Results.rds"))
          }
          message("[mETHYLotest] CNA CSV tables and RDS saved in: ", cna_dir)
        } else {
          message("[mETHYLotest] No CNA results generated (check if probe names match the array type).")
        }
      }
      message("[mETHYLotest] CNA process complete.")
    }
  } else {
    message("[mETHYLotest] CNA skipped.")
  }

  # ========================================================================
  # 17. RESULTS
  # ========================================================================

  results <- list(
    myLoad        = myLoad,
    beta_final    = beta_final,
    beta_history  = beta_label,
    pheno         = pheno_vec,
    pheno_col     = pheno_col,
    arraytype     = arraytype,
    myDMP         = myDMP,
    myDMR         = myDMR,
    myBlock       = myBlock,
    myGSEA        = myGSEA,
    myCNA         = myCNA,
    myRefBase     = myRefBase,
    config        = cfg
  )

  if (isTRUE(cfg$save_raw_obj))
    saveRDS(results, file.path(rds_dir, "pipeline_results.rds"))

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

  # ========================================================================
  # JSON EXPORTS FOR FUTURE WORKS
  # ========================================================================
  tryCatch({
    # Capture warnings
    warns <- names(warnings())
    if (is.null(warns)) warns <- list()
    jsonlite::write_json(warns, file.path(res_dir, "warnings.json"))
  }, error = function(e) {
    message("[mETHYLotest] Warning: Could not save warnings.json: ", e$message)
  })

  message("[mETHYLotest] ============================================")
  message("[mETHYLotest] Pipeline complete!")
  message("[mETHYLotest] Project:  ", cfg$project_name)
  message("[mETHYLotest] Results:  ", res_dir)
  message("[mETHYLotest] Samples:  ", ncol(beta_final))
  message("[mETHYLotest] Probes:   ", format(nrow(beta_final), big.mark = ","))
  message("[mETHYLotest] History:  ", beta_label)
  message("[mETHYLotest] ============================================")

  invisible(results)
}
