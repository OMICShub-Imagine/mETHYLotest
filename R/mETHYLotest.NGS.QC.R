#' Quality Control and Descriptive Analysis for Methylation Data
#'
#' @description
#' Performs a comprehensive QC on a methylRawList object. It generates plots and
#' aggregates statistics, optionally saving all metrics into a summary Excel file.
#'
#' @param methyl_obj A methylRawList object (filtered or not)
#' @param output_base_dir The base directory for figures
#' @param chromosomes Vector of chromosomes to analyze
#' @param current_min_cov Integer, the coverage threshold used (for logging)
#' @param save_summary Boolean. If TRUE, saves all stats tables to an Excel file.
#'
#' @return A list containing the dataframes used for the QC UI (df_meth, df_controls, df_stats)
#' @export
mETHYLotest.NGS.QC <- function(methyl_obj,
                                  output_base_dir,
                                  chromosomes = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"),
                                  current_min_cov = 1, unite_destrand = FALSE,
                                  save_summary = TRUE) {

  message(paste("--- Starting QC Analysis in:", output_base_dir, "---"))

  # --- 0. SÉCURITÉ : Vérification des données vides ---
  rows_count <- sapply(methyl_obj, nrow)
  empty_indices <- which(rows_count == 0)
  all_sample_ids <- methylKit::getSampleID(methyl_obj)

  if (length(empty_indices) > 0) {
    empty_ids <- all_sample_ids[empty_indices]
    if (length(empty_indices) == length(methyl_obj)) {
      stop(paste("CRITICAL ERROR: All samples are empty (0 CpGs found). Check min_cov."))
    }
    warning(paste("WARNING: Skipping empty samples:", paste(empty_ids, collapse = ", ")))
    methyl_obj <- methyl_obj[-empty_indices]
  }

  message(paste("Proceeding QC with", length(methyl_obj), "valid samples..."))

  # --- 1. Gestion des Dossiers ---
  dir_meth    <- file.path(output_base_dir, "Methylation_Stats")
  dir_cov     <- file.path(output_base_dir, "Coverage_Stats")
  dir_pos     <- file.path(output_base_dir, "PositionStats")
  dir_per_chr <- file.path(output_base_dir, "Per_Chromosome_Methylation")
  dir_qc      <- file.path(output_base_dir, "QC")

  directories <- c(dir_meth, dir_cov, dir_pos, dir_per_chr, dir_qc)
  for(d in directories) if(!dir.exists(d)) dir.create(d, recursive = TRUE)

  # --- 2. Plots natifs methylKit ---
  message("Generating methylKit native histograms...")
  current_sample_ids <- methylKit::getSampleID(methyl_obj)

  for(i in 1:length(methyl_obj)){
    sample_id <- current_sample_ids[i]
    if(nrow(methyl_obj[[i]]) < 5) next

    tryCatch({
      png(file.path(dir_meth, paste0(sample_id, "_methylation_stats.png")))
      methylKit::getMethylationStats(methyl_obj[[i]], plot=TRUE, both.strands=FALSE)
      dev.off()
    }, error = function(e) { if(dev.cur() > 1) dev.off() })

    tryCatch({
      png(file.path(dir_cov, paste0(sample_id, "_coverage_stats.png")))
      methylKit::getCoverageStats(methyl_obj[[i]], plot=TRUE, both.strands=FALSE)
      dev.off()
    }, error = function(e) { if(dev.cur() > 1) dev.off() })
  }

  # --- 3. Analyses Avancées ---

  # A. Pourcentage de méthylation par chromosome
  df_meth <- tryCatch({
    mETHYLotest.NGS.get_chromosome_methylation_stats(
      methyl_obj, output_dir = dir_per_chr, chromosomes = chromosomes, file_ext = ".png"
    )
  }, error = function(e) return(NULL))

  # B. Statistiques des contrôles
  df_controls <- tryCatch({
    mETHYLotest.NGS.get_control_stats(
      methyl_obj, save_plot_to = file.path(dir_qc, "control_methylation.png")
    )
  }, error = function(e) return(NULL))

  # C. Statistiques de positions
  stats_df <- tryCatch({
    mETHYLotest.NGS.get_position_counts_with_common(
      methyl_obj, output_dir = dir_pos, file_name = "positions_stats.png", unite_destrand = unite_destrand, min_cov = current_min_cov
    )
  }, error = function(e) return(NULL))

  # D. Statistiques globales C/T
  c_stats <- tryCatch({
    mETHYLotest.NGS.get_C_counts(
      methyl_obj = methyl_obj, output_dir = dir_pos, plot_name = "Coverage_Stats.png", excel_name = NULL # On gère l'excel global après
    )
  }, error = function(e) return(NULL))

  # --- 4. Sauvegarde Excel Récapitulatif ---
  if (save_summary) {
    message("Saving comprehensive QC Summary Excel...")

    # Création d'une liste nommée pour les onglets Excel
    sheets_list <- list()

    if (!is.null(df_meth))     sheets_list[["Chromosomes"]] <- df_meth
    if (!is.null(stats_df))    sheets_list[["Positions"]]   <- stats_df
    if (!is.null(c_stats))     sheets_list[["Coverage_CT"]] <- c_stats
    if (!is.null(df_controls)) sheets_list[["Controls"]]    <- df_controls

    # Ajout d'un onglet "Metadata" avec les paramètres utilisés
    meta_df <- data.frame(
      Parameter = c("Min Coverage Used", "Samples Count", "Date"),
      Value = c(current_min_cov, length(methyl_obj), as.character(Sys.time()))
    )
    sheets_list[["Run_Info"]] <- meta_df

    if (length(sheets_list) > 0) {
      summary_file <- file.path(output_base_dir, paste0("QC_Summary_Cov", current_min_cov, ".xlsx"))
      tryCatch({
        writexl::write_xlsx(sheets_list, path = summary_file)
        message(" -> Saved: ", summary_file)
      }, error = function(e) warning("Failed to save QC Summary Excel: ", e$message))
    }
  }

  # --- 5. Global methylation per sample ---
  df_global_meth <- tryCatch({
    gm_list <- lapply(seq_along(methyl_obj), function(i) {
      d   <- methylKit::getData(methyl_obj[[i]])
      sid <- methylKit::getSampleID(methyl_obj)[[i]]
      if (nrow(d) == 0) return(NULL)
      total_C <- sum(as.numeric(d$numCs))
      total_T <- sum(as.numeric(d$numTs))
      total   <- total_C + total_T
      data.frame(
        Sample          = sid,
        Global_Meth_Pct = if (total > 0) round(100 * total_C / total, 2) else NA_real_,
        Total_Cs        = total_C,
        Total_Ts        = total_T,
        Total_Reads     = total,
        N_Positions     = nrow(d),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, gm_list)
  }, error = function(e) NULL)

  # --- 6. Export JSON for future works ---
  message("Exporting JSON files for future works...")
  tryCatch({
    qc_data <- list()
    chromosomes_set <- character(0)
    current_sample_ids <- methylKit::getSampleID(methyl_obj)
    
    for (i in seq_along(methyl_obj)) {
      sid <- current_sample_ids[i]
      data <- methylKit::getData(methyl_obj[[i]])
      
      if (nrow(data) == 0) next
      
      # Statistiques globales
      total_c <- sum(as.numeric(data$numCs))
      total_t <- sum(as.numeric(data$numTs))
      total_cov <- total_c + total_t
      global_meth_pct <- if (total_cov > 0) 100 * total_c / total_cov else 0
      
      # Histogramme de méthylation
      meth_pct <- 100 * data$numCs / (data$numCs + data$numTs)
      h_meth <- hist(meth_pct, breaks=seq(0, 100, by=10), plot=FALSE)
      meth_hist <- data.frame(
        bin = h_meth$mids,
        count = h_meth$counts
      )
      
      # Histogramme de couverture (log10)
      cov_val <- data$numCs + data$numTs
      h_cov <- hist(log10(cov_val + 1), breaks=20, plot=FALSE)
      cov_hist <- data.frame(
        bin = h_cov$mids,
        count = h_cov$counts
      )
      
      # % Methylation par chromosome
      chr_meth <- aggregate(meth_pct, by=list(chr=data$chr), FUN=mean, na.rm=TRUE)
      names(chr_meth) <- c("chr", "mean_meth")
      chromosomes_set <- unique(c(chromosomes_set, as.character(chr_meth$chr)))
      
      qc_data[[sid]] <- list(
        global_meth_pct = global_meth_pct,
        total_positions = nrow(data),
        total_reads = total_cov,
        meth_hist = meth_hist,
        cov_hist = cov_hist,
        chr_meth = chr_meth
      )
    }
    
    json_out <- list(
      samples = qc_data,
      chromosomes = chromosomes_set
    )
    
    json_file <- file.path(dir_qc, "qc_results.json")
    jsonlite::write_json(json_out, json_file, pretty = TRUE, auto_unbox = TRUE)
  }, error = function(e) {
    warning("Failed to export JSON files: ", e$message)
  })

  message("QC Analysis completed.")
  return(list(
    df_meth        = df_meth,
    df_controls    = df_controls,
    df_stats       = stats_df,
    df_global_meth = df_global_meth
  ))
}
