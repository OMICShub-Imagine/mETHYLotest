#' Validate Epigenetic Signatures
#'
#' @description
#' Performs supervised (SVM) and unsupervised (PCA, Silhouette) validation
#' on user-defined or auto-generated CpG signatures.
#'
#' Signatures are \code{.txt} files containing one CpG ID per line,
#' typically created via the DMP Explorer UI.
#'
#' @param beta              Numeric matrix (probes x samples).
#' @param pheno             Phenotype vector, length == ncol(beta).
#' @param signatures_folder Path to folder containing .txt signature files.
#' @param arraytype         Array type (e.g. "EPICv1", "EPICv2", "450K").
#'                          Passed from \code{cfg$technology}.
#' @param k_folds           Number of CV folds (default 5).
#'
#' @return Named list of validation results (one entry per signature).
#' @export
#' @importFrom caret train trainControl defaultSummary twoClassSummary
#' @importFrom e1071 svm
#' @importFrom cluster silhouette
#' @importFrom stats prcomp dist
#' @importFrom utils data
mETHYLotest.EPIC.validate <- function(beta,
                                     pheno,
                                     signatures_folder,
                                     arraytype = "EPICv1",
                                     k_folds = 5) {

  for (pkg in c("caret", "e1071", "cluster", "ChAMPdata"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required for validation.", pkg))

  if (!dir.exists(signatures_folder))
    stop("[Validation] Signature folder not found: ", signatures_folder)

  sig_files <- list.files(signatures_folder, pattern = "\\.txt$",
                          full.names = TRUE)
  if (length(sig_files) == 0L)
    stop("[Validation] No .txt files found in: ", signatures_folder)

  if (ncol(beta) != length(pheno))
    stop("[Validation] Mismatch: ncol(beta)=", ncol(beta),
         " vs length(pheno)=", length(pheno))

  # ── Annotation ──────────────────────────────────────────────────────────
  message("[Validation] Loading annotation for: ", arraytype)
  probe_features <- NULL
  tryCatch({
    pkg_data <- switch(arraytype,
                       "EPICv1" = "probe.features.epicv1",
                       "EPIC"   = "probe.features.epicv1",
                       "EPICv2" = "probe.features.epicv2",
                       "450K"   = "probe.features",
                       NULL)
    if (!is.null(pkg_data)) {
      local_env <- new.env(parent = emptyenv())
      utils::data(list = pkg_data, package = "ChAMPdata",
                  envir = local_env)
      if (exists("probe.features", envir = local_env))
        probe_features <- get("probe.features", envir = local_env)
    }
  }, error = function(e)
    warning("[Validation] Annotation load failed: ", e$message))

  if (is.null(probe_features)) {
    probe_features <- data.frame(
      gene = rep("N/A", nrow(beta)),
      row.names = rownames(beta))
  } else {
    probe_features <- probe_features[, "gene", drop = FALSE]
  }

  # ── Phenotype normalisation ─────────────────────────────────────────────
  pheno_char <- as.character(pheno)
  is_ctl <- grepl("ctl|control|healthy|normal",
                  pheno_char, ignore.case = TRUE)
  if (any(is_ctl)) {
    message("[Validation] Normalising: ", sum(is_ctl),
            " samples renamed to 'CTL'.")
    pheno_char[is_ctl] <- "CTL"
  }

  pheno_fac <- as.factor(pheno_char)
  levels(pheno_fac) <- make.names(levels(pheno_fac))

  # ── CV strategy ─────────────────────────────────────────────────────────
  min_class <- min(table(pheno_fac))
  train_method  <- "cv"
  train_number  <- k_folds
  train_metric  <- "Accuracy"
  train_probs   <- FALSE
  train_summary <- caret::defaultSummary

  if (min_class < k_folds) {
    message("[Validation] Small class (n=", min_class,
            "). Switching to LOOCV.")
    train_method <- "LOOCV"
    train_number <- NULL

    if (min_class >= 5 && length(levels(pheno_fac)) == 2L) {
      train_metric  <- "ROC"
      train_probs   <- TRUE
      train_summary <- caret::twoClassSummary
    }
  } else {
    if (length(levels(pheno_fac)) == 2L) {
      train_metric  <- "ROC"
      train_probs   <- TRUE
      train_summary <- caret::twoClassSummary
    } else {
      train_probs <- TRUE
    }
  }

  # ── Process each signature ──────────────────────────────────────────────
  results <- list()
  message("[Validation] Processing ", length(sig_files), " signature(s)...")

  for (f in sig_files) {
    sig_name <- tools::file_path_sans_ext(basename(f))
    cpgs     <- trimws(readLines(f))
    cpgs     <- cpgs[nzchar(cpgs)]

    valid_cpgs <- intersect(cpgs, rownames(beta))

    if (length(valid_cpgs) < 2L) {
      warning("[Validation] Skip ", sig_name, ": < 2 valid CpGs.")
      next
    }

    # Gene mapping
    genes_vec <- setNames(rep("Unknown", length(valid_cpgs)), valid_cpgs)
    common    <- intersect(valid_cpgs, rownames(probe_features))
    if (length(common) > 0L)
      genes_vec[common] <- probe_features[common, "gene"]

    # Data prep (samples x CpGs)
    dat <- t(beta[valid_cpgs, , drop = FALSE])

    train_x <- as.data.frame(dat)
    colnames(train_x) <- make.names(colnames(train_x))

    ml_df_raw       <- data.frame(dat)
    ml_df_raw$Class <- pheno_fac

    beta_subset <- beta[valid_cpgs, , drop = FALSE]

    # ── 1. SVM ────────────────────────────────────────────────────────────
    ctrl <- caret::trainControl(
      method          = train_method,
      number          = train_number,
      classProbs      = train_probs,
      summaryFunction = train_summary,
      savePredictions = "final",
      allowParallel   = FALSE
    )

    set.seed(42L)
    svm_fit <- tryCatch({
      caret::train(
        x          = train_x,
        y          = pheno_fac,
        method     = "svmLinear",
        trControl  = ctrl,
        preProcess = c("center", "scale"),
        tuneLength = 5,
        metric     = train_metric
      )
    }, error = function(e) {
      warning("[Validation] SVM failed for ", sig_name, ": ", e$message)
      NULL
    })

    # ── 2. PCA ────────────────────────────────────────────────────────────
    pca_res <- tryCatch(
      stats::prcomp(dat, scale. = TRUE),
      error = function(e) NULL
    )

    # ── 3. Silhouette ─────────────────────────────────────────────────────
    avg_sil <- 0
    if (length(unique(pheno_fac)) > 1L) {
      tryCatch({
        dist_mat <- stats::dist(dat)
        sil      <- cluster::silhouette(as.numeric(pheno_fac), dist_mat)
        avg_sil  <- mean(sil[, 3L], na.rm = TRUE)
        if (is.na(avg_sil)) avg_sil <- 0
      }, error = function(e) NULL)
    }

    results[[sig_name]] <- list(
      n_cpgs      = length(valid_cpgs),
      cpgs        = valid_cpgs,
      genes_map   = genes_vec,
      svm_model   = svm_fit,
      pca         = pca_res,
      silhouette  = avg_sil,
      raw_data    = ml_df_raw,
      beta_matrix = beta_subset,
      metric_used = train_metric
    )

    perf_val <- "NA"
    if (!is.null(svm_fit)) {
      val      <- svm_fit$results[[train_metric]]
      perf_val <- round(max(val, na.rm = TRUE), 2)
    }

    message("  -> ", sig_name,
            " | CpGs: ", length(valid_cpgs),
            " | ", train_metric, ": ", perf_val,
            " | Silhouette: ", round(avg_sil, 2))
  }

  message("[Validation] Done. ", length(results),
          " signature(s) validated.")
  results
}
