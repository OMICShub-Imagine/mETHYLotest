#' Validate Epigenetic Signatures for WGBS Long-Read Data
#'
#' @description
#' Performs supervised (SVM) and unsupervised (PCA, Silhouette) validation
#' on CpG signatures from methylKit objects.
#'
#' Signatures are .txt files with one locus ID per line (format: chr:start-end).
#'
#' @param meth              methylBase object (united).
#' @param diff_results      Named list of methylDiff objects.
#' @param signatures_folder Path to folder containing .txt signature files.
#' @param diff_cutoff       Methylation difference cutoff (default 10).
#' @param diff_qvalue       Q-value cutoff (default 0.05).
#' @param k_folds           Number of CV folds (default 5).
#'
#' @return Named list of validation results (one entry per signature).
#' @export
mETHYLotest.NGS.validate <- function(meth,
                                        diff_results      = NULL,
                                        signatures_folder = NULL,
                                        diff_cutoff       = 10,
                                        diff_qvalue       = 0.05,
                                        k_folds           = 5) {

  for (pkg in c("caret", "e1071", "cluster")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("[Validation] Package '%s' is required. Install with: install.packages('%s')",
                   pkg, pkg))
    }
  }

  if (!inherits(meth, "methylBase")) {
    stop("[Validation] 'meth' must be a methylKit methylBase object.")
  }

  if (is.null(signatures_folder) || !dir.exists(signatures_folder)) {
    stop("[Validation] Signature folder not found: ", signatures_folder)
  }

  sig_files <- list.files(signatures_folder, pattern = "\\.txt$",
                          full.names = TRUE)
  if (length(sig_files) == 0L) {
    stop("[Validation] No .txt files found in: ", signatures_folder)
  }

  # ══════════════════════════════════════════════════════
  # 1. Extract beta matrix and phenotype from methylBase
  # ══════════════════════════════════════════════════════
  message("[Validation] Extracting data from methylBase...")

  sample_ids <- methylKit::getSampleID(meth)
  treatment  <- meth@treatment
  meth_data  <- methylKit::getData(meth)

  # Beta values (0-1)
  beta_pct <- methylKit::percMethylation(meth, rowids = TRUE)
  beta <- beta_pct / 100

  if (!is.matrix(beta)) beta <- as.matrix(beta)

  # Rownames as chr:start-end (matching signature format)
  rn <- paste0(meth_data$chr, ":", meth_data$start, "-", meth_data$end)
  rownames(beta) <- rn

  # Phenotype: treatment 0 -> "CTL", 1 -> "Test"
  pheno_char <- ifelse(treatment == 0, "CTL", "Test")
  pheno_fac  <- as.factor(pheno_char)
  levels(pheno_fac) <- make.names(levels(pheno_fac))

  message(sprintf("[Validation]   %s CpGs x %d samples (%s)",
                  format(nrow(beta), big.mark = ","),
                  ncol(beta),
                  paste(table(pheno_char), collapse = " CTL / ") ))

  # ══════════════════════════════════════════════════════
  # 2. Build gene annotation from diff_results if available
  # ══════════════════════════════════════════════════════
  # WGBS doesn't have array annotation, but we can store
  # genomic coordinates as "gene" equivalent
  probe_info <- data.frame(
    locus = rn,
    chr   = as.character(meth_data$chr),
    start = meth_data$start,
    end   = meth_data$end,
    row.names = rn,
    stringsAsFactors = FALSE
  )

  # If annotatr or another annotation is available, try to get genes
  has_annotatr <- requireNamespace("annotatr", quietly = TRUE)

  # ══════════════════════════════════════════════════════
  # 3. Determine CV strategy
  # ══════════════════════════════════════════════════════
  min_class <- min(table(pheno_fac))
  n_samples <- length(pheno_fac)

  train_method  <- "cv"
  train_number  <- k_folds
  train_metric  <- "Accuracy"
  train_probs   <- FALSE
  train_summary <- caret::defaultSummary

  if (min_class < k_folds) {
    message("[Validation] Small class (n=", min_class, "). Switching to LOOCV.")
    train_method <- "LOOCV"
    train_number <- NULL

    if (min_class >= 3 && length(levels(pheno_fac)) == 2L) {
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

  message(sprintf("[Validation] CV strategy: %s (metric: %s)",
                  train_method, train_metric))

  # ══════════════════════════════════════════════════════
  # 4. Process each signature
  # ══════════════════════════════════════════════════════
  results <- list()
  message(sprintf("[Validation] Processing %d signature(s)...", length(sig_files)))

  for (f in sig_files) {
    sig_name <- tools::file_path_sans_ext(basename(f))
    cpgs     <- trimws(readLines(f, warn = FALSE))
    cpgs     <- cpgs[nzchar(cpgs)]

    if (length(cpgs) == 0L) {
      warning("[Validation] Skip ", sig_name, ": empty file.")
      next
    }

    # ── Match CpGs to beta rownames ──
    # Signatures can be in format chr:start-end or chr:start
    # Beta rownames are chr:start-end
    valid_cpgs <- intersect(cpgs, rownames(beta))

    # If no direct match, try chr:start format
    if (length(valid_cpgs) == 0L) {
      # Build chr:start lookup from beta rownames
      beta_start <- sub("-[0-9]+$", "", rownames(beta))
      sig_start  <- sub("-[0-9]+$", "", cpgs)

      match_idx <- match(sig_start, beta_start)
      matched   <- !is.na(match_idx)

      if (any(matched)) {
        valid_cpgs <- rownames(beta)[match_idx[matched]]
        message(sprintf("  [%s] Matched %d/%d via chr:start format.",
                        sig_name, length(valid_cpgs), length(cpgs)))
      }
    }

    if (length(valid_cpgs) < 2L) {
      warning("[Validation] Skip ", sig_name, ": only ",
              length(valid_cpgs), " valid CpGs (need >= 2).")
      next
    }

    message(sprintf("  [%s] %d/%d CpGs matched...",
                    sig_name, length(valid_cpgs), length(cpgs)),
            appendLF = FALSE)

    # ── Gene mapping (genomic coordinates) ──
    genes_vec <- setNames(rep("N/A", length(valid_cpgs)), valid_cpgs)
    if (all(valid_cpgs %in% rownames(probe_info))) {
      # Use coordinates as "gene" for now
      genes_vec <- setNames(
        paste0(probe_info[valid_cpgs, "chr"], ":",
               probe_info[valid_cpgs, "start"]),
        valid_cpgs
      )
    }

    # ── Data prep (samples x CpGs) ──
    beta_subset <- beta[valid_cpgs, , drop = FALSE]

    # Remove CpGs with zero variance (cause SVM to fail)
    cpg_var <- apply(beta_subset, 1, var, na.rm = TRUE)
    cpg_var[is.na(cpg_var)] <- 0
    good_var <- cpg_var > 1e-10

    if (sum(good_var) < 2L) {
      warning("[Validation] Skip ", sig_name, ": < 2 variable CpGs.")
      next
    }

    if (sum(!good_var) > 0) {
      message(sprintf(" (%d zero-var removed)", sum(!good_var)),
              appendLF = FALSE)
      valid_cpgs  <- valid_cpgs[good_var]
      beta_subset <- beta_subset[good_var, , drop = FALSE]
    }

    # Impute NAs with row means
    if (any(is.na(beta_subset))) {
      rm <- rowMeans(beta_subset, na.rm = TRUE)
      rm[is.na(rm)] <- 0.5
      na_idx <- which(is.na(beta_subset), arr.ind = TRUE)
      if (nrow(na_idx) > 0) {
        beta_subset[na_idx] <- rm[na_idx[, 1]]
      }
    }

    dat <- t(beta_subset)  # samples x CpGs

    train_x <- as.data.frame(dat)
    colnames(train_x) <- make.names(colnames(train_x))

    ml_df_raw       <- data.frame(dat)
    ml_df_raw$Class <- pheno_fac

    # ── 4a. SVM ──
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

    # ── 4b. PCA ──
    pca_res <- tryCatch(
      stats::prcomp(dat, scale. = TRUE, center = TRUE),
      error = function(e) NULL
    )

    # ── 4c. Silhouette ──
    avg_sil <- 0
    sil_obj <- NULL
    if (length(unique(pheno_fac)) > 1L && nrow(dat) >= 3) {
      tryCatch({
        dist_mat <- stats::dist(dat)
        sil_obj  <- cluster::silhouette(as.numeric(pheno_fac), dist_mat)
        avg_sil  <- mean(sil_obj[, 3L], na.rm = TRUE)
        if (is.na(avg_sil)) avg_sil <- 0
      }, error = function(e) NULL)
    }

    # ── 4d. Extract performance metrics ──
    perf_val   <- NA
    accuracy   <- NA
    auc_val    <- NA
    sens_val   <- NA
    spec_val   <- NA
    conf_mat   <- NULL
    fold_accs  <- NULL

    if (!is.null(svm_fit)) {
      # Main metric
      perf_val <- round(max(svm_fit$results[[train_metric]], na.rm = TRUE), 4)

      # Confusion matrix from predictions
      if (!is.null(svm_fit$pred)) {
        pred_df <- svm_fit$pred
        tryCatch({
          cm <- caret::confusionMatrix(pred_df$pred, pred_df$obs)
          conf_mat <- cm$table
          accuracy <- round(cm$overall["Accuracy"], 4)

          if (length(levels(pheno_fac)) == 2L) {
            sens_val <- round(cm$byClass["Sensitivity"], 4)
            spec_val <- round(cm$byClass["Specificity"], 4)
          }
        }, error = function(e) NULL)
      }

      # Per-fold accuracy (if CV)
      if (train_method == "cv" && !is.null(svm_fit$resample)) {
        if ("Accuracy" %in% colnames(svm_fit$resample)) {
          fold_accs <- svm_fit$resample$Accuracy
        }
      }

      # AUC
      if (train_metric == "ROC") {
        auc_val <- perf_val
      }
    }

    # ── Store results ──
    results[[sig_name]] <- list(
      # Core data
      n_cpgs          = length(valid_cpgs),
      cpgs            = valid_cpgs,
      genes_map       = genes_vec,
      beta_matrix     = beta_subset,
      raw_data        = ml_df_raw,

      # Models
      svm_model       = svm_fit,
      pca             = pca_res,

      # Metrics
      silhouette      = avg_sil,
      metric_used     = train_metric,
      accuracy        = accuracy,
      auc             = if (!is.na(auc_val)) auc_val else NULL,
      sensitivity     = if (!is.na(sens_val)) sens_val else NULL,
      specificity     = if (!is.na(spec_val)) spec_val else NULL,
      confusion       = conf_mat,
      fold_accuracies = fold_accs
    )

    # ── Console summary ──
    sil_str  <- sprintf("%.3f", avg_sil)
    perf_str <- if (!is.na(perf_val)) sprintf("%.3f", perf_val) else "NA"
    acc_str  <- if (!is.na(accuracy)) sprintf("%.1f%%", accuracy * 100) else "NA"

    message(sprintf(" %s: %s | Sil: %s | Acc: %s",
                    train_metric, perf_str, sil_str, acc_str))
  }

  # ══════════════════════════════════════════════════════
  # 5. Summary
  # ══════════════════════════════════════════════════════
  message(sprintf("\n[Validation] === Complete: %d/%d signatures validated ===",
                  length(results), length(sig_files)))

  if (length(results) > 0) {
    summ <- do.call(rbind, lapply(names(results), function(nm) {
      r <- results[[nm]]
      data.frame(
        Signature  = nm,
        CpGs       = r$n_cpgs,
        Silhouette = sprintf("%.3f", r$silhouette),
        Metric     = r$metric_used,
        Value      = if (!is.null(r$svm_model)) {
          sprintf("%.3f", max(r$svm_model$results[[r$metric_used]], na.rm = TRUE))
        } else { "NA" },
        Accuracy   = if (!is.na(r$accuracy)) sprintf("%.1f%%", r$accuracy * 100) else "NA",
        stringsAsFactors = FALSE
      )
    }))
    message("\n", paste(utils::capture.output(print(summ, row.names = FALSE)),
                        collapse = "\n"))
  }

  return(results)
}
