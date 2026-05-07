#' Generate HTML Reports for DMP Signatures
#'
#' @description
#' Generates HTML reports. Can works in two modes:
#' 1. Standard Mode: Iterates over dmp_results, filters by adjPVal, creates reports.
#' 2. Signature Mode: If 'signatures_folder' is provided, generates reports for each .txt signature file found.
#'
#' @param dmp_results List of dataframes (ChAMP output).
#' @param beta Numeric matrix of Beta values.
#' @param pheno Vector of phenotype labels.
#' @param output_directory Path to save HTML files.
#' @param signatures_folder (Optional) Path to folder containing .txt signature files.
#' @param signature_comparison_source (Optional) Name of the comparison in 'dmp_results' to use for statistics when in Signature Mode.
#' @param adjPVal Threshold for standard mode.
#' @param arraytype "EPICv1", "EPICv2", or "450K".
#' @param rmd_template_path Path to the Rmd file (optional override).
#'
#' @export
mETHYLotest.EPIC.DMP.report <- function(dmp_results,
                                       beta,
                                       pheno,
                                       output_directory = file.path(getwd(), "EpicPipeline", "DMP_Reports"),
                                       signatures_folder = NULL,
                                       signature_comparison_source = NULL,
                                       adjPVal = 0.05,
                                       arraytype = "EPICv1",
                                       top_n_heatmap = 50,
                                       top_n_windows = 20,
                                       rmd_template_path = NULL) {

  # --- Checks ---
  if (is.null(dmp_results) || length(dmp_results) == 0) stop("'dmp_results' empty.")
  if (is.null(beta) || is.null(pheno)) stop("'beta' and 'pheno' are required.")

  # Packages check
  req_pkgs <- c("rmarkdown", "knitr", "ggplot2", "heatmaply", "DT", "dplyr", "tidyr", "htmltools", "ChAMPdata")
  missing <- req_pkgs[!req_pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) stop(paste("Missing packages:", paste(missing, collapse = ", ")))

  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

  # --- 1. CHARGEMENT DE L'ANNOTATION ---
  message("[mETHYLotest] Loading full annotation...")
  probe_features <- NULL
  tryCatch({
    if(arraytype %in% c("EPICv1", "EPIC")) {
      utils::data("probe.features.epicv1", package = "ChAMPdata", envir = environment())
      if(exists("probe.features")) probe_features <- probe.features
    } else if (arraytype == "EPICv2") {
      utils::data("probe.features.epicv2", package = "ChAMPdata", envir = environment())
      if(exists("probe.features")) probe_features <- probe.features
    } else if (arraytype == "450K") {
      utils::data("probe.features", package = "ChAMPdata", envir = environment())
      if(exists("probe.features")) probe_features <- probe.features
    }
  }, error = function(e) warning("Annotation loading failed."))

  if(is.null(probe_features)) stop("Could not load ChAMPdata annotation. Check 'arraytype'.")

  # Optimisation : On allège l'annotation
  probe_features <- probe_features[, c("CHR", "MAPINFO", "gene")]
  probe_features$CpG <- rownames(probe_features)

  # --- 2. LOCATE RMD TEMPLATE ---
  if (is.null(rmd_template_path)) {
    rmd_template <- system.file("rmarkdown",
                                "mETHYLotest.EPIC.signature_report.Rmd",
                                package = "mETHYLotest")
    if (rmd_template == "")
      stop("[mETHYLotest] Template not found in package. ",
           "Check inst/rmarkdown/mETHYLotest.EPIC.signature_report.Rmd")
  } else {
    rmd_template <- rmd_template_path
  }

  if (!file.exists(rmd_template))
    stop("[mETHYLotest] RMD template not found: ", rmd_template)

  # --- 3. WORK ---
  work_list <- list()

  # Mode SIGNATURE (Prioritaire si dossier fourni)
  if (!is.null(signatures_folder)) {
    if(!dir.exists(signatures_folder)) stop("Provided 'signatures_folder' does not exist.")

    sig_files <- list.files(signatures_folder, pattern = "\\.txt$", full.names = TRUE)
    if(length(sig_files) == 0) stop("No .txt files found in signatures_folder.")

    # Détermination de la source des stats (Comparison)
    stats_source <- signature_comparison_source
    if(is.null(stats_source)) {
      stats_source <- names(dmp_results)[1]
      warning(paste("No 'signature_comparison_source' provided. Using the first comparison available for statistics:", stats_source))
    }
    if(!stats_source %in% names(dmp_results)) stop(paste("Comparison", stats_source, "not found in dmp_results."))

    message(paste("[mETHYLotest] Mode: Signature Folder. Processing", length(sig_files), "signatures using stats from:", stats_source))

    for (f in sig_files) {
      sig_name <- tools::file_path_sans_ext(basename(f))
      cpgs <- readLines(f)
      cpgs <- trimws(cpgs[cpgs != ""]) # Clean empty lines

      work_list[[sig_name]] <- list(
        type = "signature",
        comparison_name = stats_source,
        target_cpgs = cpgs
      )
    }

  } else {
    # Mode STANDARD (Parcours de toutes les comparaisons)
    message("[mETHYLotest] Mode: Standard (All Comparisons).")
    for (comp_name in names(dmp_results)) {
      work_list[[comp_name]] <- list(
        type = "standard",
        comparison_name = comp_name,
        target_cpgs = NULL # Sera déterminé par adjPVal
      )
    }
  }

  # --- 4. BOUCLE DE GÉNÉRATION ---

  for (report_name in names(work_list)) {
    task <- work_list[[report_name]]
    comp_name <- task$comparison_name

    message(paste0("Processing: ", report_name, " (Source: ", comp_name, ")..."))

    # 1. Récupération des données brutes
    dmp_full <- dmp_results[[comp_name]]

    # --- Nettoyage Colonnes ---
    if (!"CpG" %in% colnames(dmp_full)) dmp_full$CpG <- rownames(dmp_full)
    if (!"deltaBeta" %in% colnames(dmp_full) && "logFC" %in% colnames(dmp_full)) dmp_full$deltaBeta <- dmp_full$logFC
    if ("cgi" %in% colnames(dmp_full)) colnames(dmp_full)[colnames(dmp_full) == "cgi"] <- "CGI"
    if (!"CGI" %in% colnames(dmp_full)) dmp_full$CGI <- "N/A"
    if (!"feature" %in% colnames(dmp_full)) dmp_full$feature <- "N/A"

    # 2. Filtrage (Signature vs Standard)
    if (task$type == "signature") {
      # On ne garde que les CpGs de la signature
      # Intersection pour éviter les erreurs si un CpG n'est pas dans le résultat
      valid_cpgs <- intersect(task$target_cpgs, rownames(dmp_full))

      if(length(valid_cpgs) == 0) {
        warning(paste("  -> No matching CpGs found in dmp_results for signature:", report_name))
        next
      }
      dmp_sig <- dmp_full[valid_cpgs, ]

    } else {
      # Standard : filtre par P-value
      dmp_sig <- dmp_full[dmp_full$adj.P.Val < adjPVal, ]
    }

    if (nrow(dmp_sig) == 0) {
      warning(sprintf("  -> No DMPs found for %s. Skipping.", report_name))
      next
    }

    # Ajout du Status (Hyper/Hypo)
    dmp_sig$Status <- ifelse(dmp_sig$deltaBeta > 0, "Hyper", "Hypo")

    # 3. Beta Matrix Subset
    # On prend tous les rownames de l'annotation qui sont dans beta (contexte global pour heatmap/windows)
    valid_beta_probes <- intersect(rownames(beta), rownames(probe_features))
    beta_subset <- beta[valid_beta_probes, , drop = FALSE]

    # Construction de l'objet de données
    report_data <- list(
      signature_name = report_name, # Nom du rapport = Nom de la signature OU de la comparaison
      dmp_full = dmp_full,
      dmp_sig = dmp_sig,
      full_annotation = probe_features,
      beta_matrix = beta_subset,
      pheno = pheno,
      adjPVal = adjPVal,
      top_n_heatmap = top_n_heatmap,
      top_n_windows = top_n_windows
    )

    # Sauvegarde temporaire pour Rmarkdown
    safe_name <- gsub("[^[:alnum:]_]", "_", report_name)
    temp_rds <- file.path(output_directory, paste0(".temp_", safe_name, ".rds"))
    saveRDS(report_data, temp_rds)

    html_output <- file.path(output_directory, paste0("Report_", safe_name, ".html"))

    tryCatch({
      rmarkdown::render(
        input = rmd_template,
        output_file = basename(html_output),
        output_dir = output_directory,
        params = list(data_path = temp_rds),
        envir = new.env(),
        quiet = TRUE
      )
      message(sprintf("  -> [Success] Created: %s", basename(html_output)))
    }, error = function(e) {
      message(paste("  -> [Error] Failed to render:", e$message))
    }, finally = {
      if (file.exists(temp_rds)) file.remove(temp_rds)
    })
  }

  message("[mETHYLotest] Reporting complete.")
}
