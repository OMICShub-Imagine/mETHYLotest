#' Differential Methylation Position (DMP) Pipeline
#'
#' @description
#' Wrapper for \code{ChAMP::champ.DMP}. Handles phenotype normalisation (CTL)
#' and flexible exports (CSV, Excel, RDS).
#'
#' @param beta          Numeric matrix (probes x samples).
#' @param pheno         Phenotype vector, length == ncol(beta).
#' @param compare.group List of 2-element character vectors for pairwise
#'   comparisons. \code{NULL} = all unique pairs derived from \code{pheno}.
#' @param adjPVal       Adjusted p-value threshold (default \code{0.05}).
#' @param p.val         Raw p-value threshold (default \code{0.05}).
#' @param adjust.method Multiple testing correction method (default \code{"BH"}).
#' @param arraytype     \code{"EPICv1"}, \code{"EPICv2"}, or \code{"450K"}.
#' @param output_directory Default output directory.
#' @param ctl_label_standardize Logical. Normalise CTL/control/normal to
#'   \code{"CTL"} (default \code{TRUE}).
#' @param save_csv      \code{TRUE} / path / \code{FALSE}.
#' @param save_excel    \code{TRUE} / path / \code{FALSE}.
#' @param save_RDS      \code{TRUE} / path / \code{NULL}.
#'
#' @return Named list of DMP data frames (one per comparison), invisibly.
#' @export
#' @importFrom ChAMP champ.DMP
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#'   createStyle addStyle
#' @importFrom utils write.csv combn
mETHYLotest.EPIC.DMP <- function(beta             = NULL,
                                pheno            = NULL,
                                compare.group    = NULL,
                                adjPVal          = 0.05,
                                p.val            = 0.05,    # fixed: was missing
                                adjust.method    = "BH",
                                arraytype        = "EPICv1",
                                output_directory = file.path(getwd(), "DMP_Results"),
                                ctl_label_standardize = TRUE,
                                save_csv         = TRUE,
                                save_excel       = TRUE,
                                save_RDS         = NULL) {

  # ── Input checks ───────────────────────────────────────────────────────────
  if (is.null(beta) || is.null(pheno))
    stop("'beta' and 'pheno' are required.")
  if (ncol(beta) != length(pheno))
    stop("Dimension mismatch: ncol(beta) = ", ncol(beta),
         ", length(pheno) = ", length(pheno), ".")

  if (!dir.exists(output_directory))
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  # ── Path resolution ────────────────────────────────────────────────────────
  resolve_path <- function(arg, default_dir) {
    if (is.null(arg) || isFALSE(arg)) return(NULL)
    if (isTRUE(arg))                  return(default_dir)
    if (is.character(arg)) {
      if (!dir.exists(arg)) dir.create(arg, recursive = TRUE,
                                       showWarnings = FALSE)
      return(arg)
    }
    warning("[mETHYLotest DMP] Unrecognised save argument — ignored.")
    NULL
  }

  path_csv   <- resolve_path(save_csv,   output_directory)
  path_excel <- resolve_path(save_excel, output_directory)
  path_rds   <- resolve_path(save_RDS,   output_directory)

  message("[mETHYLotest DMP] Starting DMP analysis.")

  # ── Phenotype normalisation ────────────────────────────────────────────────
  if (ctl_label_standardize)
    pheno <- ifelse(grepl("ctl|control|normal", pheno, ignore.case = TRUE),
                    "CTL", pheno)

  levels_avail <- unique(pheno[!is.na(pheno)])

  if (length(levels_avail) < 2L)
    stop("[mETHYLotest DMP] At least 2 phenotype groups required. Found: ",
         paste(levels_avail, collapse = ", "))

  # ── compare.group: build all pairwise pairs if NULL ───────────────────────
  # ChAMP expects a list of length-2 character vectors.
  # fixed: was unique(pheno) which returns all levels, not pairwise lists
  if (is.null(compare.group)) {
    if (length(levels_avail) == 2L) {
      compare.group <- as.character(levels_avail)  # c("CTL", "ATRX") — vecteur simple
    } else {
      # > 2 groupes : laisser NULL, ChAMP fait toutes les paires automatiquement
      compare.group <- NULL
      message("[mETHYLotest DMP] ", length(levels_avail),
              " groups detected — ChAMP will perform all pairwise comparisons: ",
              paste(levels_avail, collapse = ", "))
    }
  }

  message("[mETHYLotest DMP] adjPVal = ", adjPVal,
          " | p.val = ", p.val,
          " | method = ", adjust.method,
          " | arraytype = ", arraytype)

  # ── Run ChAMP ─────────────────────────────────────────────────────────────
  dmp_results <- ChAMP::champ.DMP(
    beta          = beta,
    pheno         = pheno,
    compare.group = compare.group,
    adjPVal       = adjPVal,
    adjust.method = adjust.method,
    arraytype     = arraytype
  )

  if (length(dmp_results) == 0L) {
    warning("[mETHYLotest DMP] No DMPs identified.")
    return(invisible(dmp_results))
  }

  message("[mETHYLotest DMP] ", length(dmp_results),
          " comparison(s) returned results.")

  # ── Export ─────────────────────────────────────────────────────────────────
  for (group in names(dmp_results)) {

    safe_name <- gsub("[^[:alnum:]_]", "_", group)
    df        <- dmp_results[[group]]

    # Normalise column names
    if (!"deltaBeta" %in% colnames(df) && "logFC" %in% colnames(df))
      df$deltaBeta <- df$logFC
    if ("cgi" %in% colnames(df))
      colnames(df)[colnames(df) == "cgi"] <- "CGI"

    # CpG ID as first column
    df <- data.frame(CpG = rownames(df), df, check.names = FALSE)

    n_sig   <- sum(df$adj.P.Val < adjPVal, na.rm = TRUE)
    n_hyper <- sum(df$adj.P.Val < adjPVal & df$deltaBeta > 0, na.rm = TRUE)
    n_hypo  <- sum(df$adj.P.Val < adjPVal & df$deltaBeta < 0, na.rm = TRUE)

    message("[mETHYLotest DMP] ", group,
            " | total: ", nrow(df),
            " | significant: ", n_sig,
            " (Hyper: ", n_hyper, ", Hypo: ", n_hypo, ")")

    # ── CSV ──────────────────────────────────────────────────────────────────
    if (!is.null(path_csv)) {
      fpath <- file.path(path_csv, paste0("DMP_", safe_name, ".csv"))
      utils::write.csv(df, file = fpath, row.names = FALSE)
      message("[mETHYLotest DMP] CSV saved: ", fpath)
    }

    # ── RDS (per comparison) ─────────────────────────────────────────────────
    if (!is.null(path_rds)) {
      fpath <- file.path(path_rds, paste0("DMP_", safe_name, ".rds"))
      saveRDS(df, file = fpath)
    }

    # ── Excel ─────────────────────────────────────────────────────────────────
    if (!is.null(path_excel)) {
      xlsx_path <- file.path(path_excel, paste0("DMP_", safe_name, ".xlsx"))

      df_sig   <- df[df$adj.P.Val < adjPVal, ]
      df_hyper <- df_sig[df_sig$deltaBeta >  0, ]
      df_hypo  <- df_sig[df_sig$deltaBeta <  0, ]

      # Beta subset: significant probes only (not all probes — avoids huge files)
      # fixed: was all probes in the DMP result which can be > 800k rows
      sig_probes  <- df_sig$CpG
      valid_sig   <- intersect(sig_probes, rownames(beta))
      beta_export <- if (length(valid_sig) > 0L) {
        sub <- beta[valid_sig, , drop = FALSE]
        data.frame(CpG = rownames(sub), sub, check.names = FALSE)
      } else {
        data.frame(CpG = character(0))
      }

      df_summary <- data.frame(
        Category = c("Total Probes Tested", "Significant (adj.P < threshold)",
                     "Hyper-methylated", "Hypo-methylated"),
        Count    = c(nrow(df), n_sig, n_hyper, n_hypo),
        Criteria = c("—",
                     paste("adj.P.Val <", adjPVal),
                     "deltaBeta > 0",
                     "deltaBeta < 0"),
        stringsAsFactors = FALSE
      )

      wb <- openxlsx::createWorkbook()
      sheets <- list(
        "Summary"          = df_summary,
        "Significant_DMPs" = df_sig,
        "Hyper_Methylated" = df_hyper,
        "Hypo_Methylated"  = df_hypo,
        "Beta_Significant" = beta_export,
        "Full_Raw_Data"    = df
      )

      bold_style <- openxlsx::createStyle(textDecoration = "bold")

      for (sname in names(sheets)) {
        openxlsx::addWorksheet(wb, sname)
        openxlsx::writeData(wb, sname, sheets[[sname]])
        openxlsx::addStyle(wb, sname, style = bold_style,
                           rows = 1L,
                           cols = seq_len(ncol(sheets[[sname]])))
      }

      openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
      message("[mETHYLotest DMP] Excel saved: ", xlsx_path)
    }
  }

  # ── Global RDS ────────────────────────────────────────────────────────────
  if (!is.null(path_rds)) {
    fpath_all <- file.path(path_rds, "DMP_All_Results.rds")
    saveRDS(dmp_results, file = fpath_all)
    message("[mETHYLotest DMP] Full RDS saved: ", fpath_all)
  }

  invisible(dmp_results)
}
