#' Annotate Differentially Methylated Cytosines (DMCs)
#'
#' @description
#' Annotates one or multiple `methylDiff` objects using the `annotatr` package.
#' Maps DMCs to gene features (promoters, exons, introns) and CpG features
#' (islands, shores, shelves). Generates summary plots and exports results to Excel.
#'
#' @param diff_obj A single `methylDiff` object or a named list of `methylDiff` objects.
#' @param assembly Character. Genome assembly (e.g., "hg38", "GRCh38", "mm10").
#' @param output_dir Character. Absolute path to save Excel and Plot files.
#' @param diff_cutoff Numeric. Minimum methylation difference percentage (default: 25).
#' @param qval_cutoff Numeric. Maximum q-value/FDR for significance (default: 0.05).
#'
#' @return A named list of dataframes containing the annotated results.
#' @export
#' @importFrom methylKit getMethylDiff
#' @importFrom writexl write_xlsx
#' @importFrom ggplot2 ggsave
mETHYLotest.NGS.AnnotateDMCs <- function(diff_obj,
                                            assembly,
                                            output_dir,
                                            diff_cutoff = 25,
                                            qval_cutoff = 0.05) {

  message("\n=== Starting Genomic Annotation ===")

  if (!requireNamespace("annotatr", quietly = TRUE)) {
    stop("Package 'annotatr' is required. Install via BiocManager::install('annotatr')")
  }
  if (!requireNamespace("methylKit", quietly = TRUE)) stop("Package 'methylKit' is required.")
  if (!requireNamespace("writexl", quietly = TRUE)) stop("Package 'writexl' is required.")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Input Standardization
  diff_list <- list()
  if (inherits(diff_obj, "methylDiff")) {
    diff_list[["Annotated_DMCs"]] <- diff_obj
  } else if (is.list(diff_obj)) {
    diff_list <- diff_obj
  } else {
    stop("Error: 'diff_obj' must be a 'methylDiff' object or a list of 'methylDiff' objects.")
  }

  # Assembly Mapping
  raw_assembly <- tolower(trimws(assembly))
  assembly_map <- c(
    "grch38" = "hg38", "grch37" = "hg19",
    "grcm38" = "mm10", "grcm39" = "mm39"
  )

  my_genome <- ifelse(raw_assembly %in% names(assembly_map), assembly_map[raw_assembly], raw_assembly)

  supported_genomes <- annotatr::builtin_genomes()
  if (!(my_genome %in% supported_genomes)) {
    stop(paste("Error: Assembly '", my_genome, "' is not supported by annotatr.\n",
               "Supported genomes include: ", paste(supported_genomes, collapse=", "), sep=""))
  }

  message("Genome assembly mapped: ", my_genome)

  # Database Construction
  annot_types <- c(paste0(my_genome, "_basicgenes"), paste0(my_genome, "_cpgs"))

  annotations <- tryCatch({
    annotatr::build_annotations(genome = my_genome, annotations = annot_types)
  }, error = function(e) {
    stop(paste("Failed to build annotations. Missing TxDb package?\nError:", e$message))
  })

  annotated_results <- list()

  # Annotation Loop
  for (model_name in names(diff_list)) {

    safe_name <- gsub("[ ()/]", "_", model_name)
    message(paste("\n-> Processing:", model_name))

    sig_diff_obj <- methylKit::getMethylDiff(diff_list[[model_name]],
                                             difference = diff_cutoff,
                                             qvalue = qval_cutoff)

    if (nrow(sig_diff_obj) == 0) {
      message("   No significant DMCs found. Skipping.")
      next
    }

    dmc_gr <- as(sig_diff_obj, "GRanges")

    message("   Mapping DMCs to genomic features...")
    dmc_annotated <- annotatr::annotate_regions(
      regions = dmc_gr,
      annotations = annotations,
      ignore.strand = TRUE,
      quiet = TRUE
    )

    # Export to Dataframe & Excel
    df_annotated <- data.frame(dmc_annotated)
    annotated_results[[model_name]] <- df_annotated

    xlsx_path <- file.path(output_dir, paste0("Annotated_DMCs_", safe_name, ".xlsx"))
    writexl::write_xlsx(df_annotated, xlsx_path)

    # Generate Plots
    tryCatch({
      p_genes <- annotatr::plot_annotation(
        annotated_regions = dmc_annotated,
        annotation_order = c(paste0(my_genome, "_promoters"), paste0(my_genome, "_5UTRs"),
                             paste0(my_genome, "_exons"), paste0(my_genome, "_introns"),
                             paste0(my_genome, "_3UTRs"), paste0(my_genome, "_intergenic")),
        plot_title = paste("Gene Features -", model_name)
      )

      p_cpgs <- annotatr::plot_annotation(
        annotated_regions = dmc_annotated,
        annotation_order = c(paste0(my_genome, "_cpg_islands"), paste0(my_genome, "_cpg_shores"),
                             paste0(my_genome, "_cpg_shelves"), paste0(my_genome, "_cpg_inter")),
        plot_title = paste("CpG Features -", model_name)
      )

      ggplot2::ggsave(filename = file.path(output_dir, paste0("Plot_Genes_", safe_name, ".png")), plot = p_genes, width = 8, height = 6)
      ggplot2::ggsave(filename = file.path(output_dir, paste0("Plot_CpGs_", safe_name, ".png")), plot = p_cpgs, width = 8, height = 6)

    }, error = function(e) {
      message("   Warning: Failed to generate plots. ", e$message)
    })
  }

  message("\n=== Genomic Annotation Completed ===")
  return(invisible(annotated_results))
}
