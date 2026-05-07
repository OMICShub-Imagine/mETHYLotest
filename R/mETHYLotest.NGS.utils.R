#' Calculate and optionally plot methylation stats per chromosome
#'
#' Computes the weighted average methylation percentage for each chromosome.
#' If specific chromosomes are requested, it filters the data. Otherwise, it uses all available chromosomes.
#' Returns a consolidated dataframe and optionally saves bar plots.
#'
#' @param methyl_obj A `methylRawList` object or list of `methylRaw` objects (from methylKit).
#' @param output_dir Directory path to save the plots. If NULL (default), plots are not generated.
#' @param chromosomes Character vector of chromosomes to analyze (e.g. c("chr1", "chr2")).
#'                    If NULL (default), all chromosomes present in the data are preserved.
#' @param file_ext Image file extension (e.g., ".png", ".pdf").
#'
#' @return A data.frame containing Sample, Chromosome, and Average Methylation %.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom purrr map_dfr
#' @importFrom methylKit getData
#' @export
mETHYLotest.NGS.get_chromosome_methylation_stats <- function(methyl_obj,
                                                             output_dir = NULL,
                                                             chromosomes = NULL,
                                                             file_ext = ".png") {

  # 1. Calculate statistics for all samples into one master dataframe
  all_stats_df <- purrr::map_dfr(methyl_obj, function(x) {

    df <- methylKit::getData(x)
    sample_id <- x@sample.id

    # Filter logic: Only filter if 'chromosomes' is not NULL
    if (!is.null(chromosomes)) {
      df <- df %>% dplyr::filter(.data$chr %in% chromosomes)
    }

    # If df is empty after filtering (e.g. wrong chromosome names), return empty safely
    if (nrow(df) == 0) return(data.frame())

    # Calculate weighted average
    df %>%
      dplyr::group_by(.data$chr) %>%
      dplyr::summarise(
        avg_methylation = 100 * sum(.data$numCs) / (sum(.data$numCs) + sum(.data$numTs)),
        .groups = 'drop'
      ) %>%
      dplyr::mutate(Sample = sample_id)
  })

  # Safety check: If result is empty, stop here to avoid errors downstream
  if (nrow(all_stats_df) == 0) {
    warning("No data found for the specified chromosomes. Returning empty dataframe.")
    return(all_stats_df)
  }

  # 2. Handle Factor Levels (Ordering)
  if (!is.null(chromosomes)) {
    plot_levels <- chromosomes
  } else {
    plot_levels <- sort(unique(all_stats_df$chr))
  }

  # Ensure only levels actually present in data are used to avoid errors if a chromosome is missing
  plot_levels <- intersect(plot_levels, unique(all_stats_df$chr))
  all_stats_df$chr <- factor(all_stats_df$chr, levels = plot_levels)


  # 3. Plotting logic
  if (!is.null(output_dir)) {

    if (!dir.exists(output_dir)) {
      message("Creating output directory: ", output_dir)
      dir.create(output_dir, recursive = TRUE)
    }

    message("Generating plots per chromosome...")

    global_max <- max(all_stats_df$avg_methylation, na.rm = TRUE)
    plot_y_limit <- if (is.finite(global_max)) global_max * 1.05 else 100

    unique_samples <- unique(all_stats_df$Sample)

    for (sid in unique_samples) {

      data_to_plot <- all_stats_df %>% dplyr::filter(.data$Sample == sid)

      # Skip if no data for this specific sample
      if(nrow(data_to_plot) == 0) next

      p <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = .data$chr, y = .data$avg_methylation, fill = .data$chr)) +
        ggplot2::geom_bar(stat = "identity", color = "black") +
        ggplot2::labs(
          title = "Methylation Percentage per Chromosome",
          subtitle = paste("Sample:", sid),
          x = "Chromosome",
          y = "Average Methylation (%)"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, plot_y_limit)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none"
        )

      filename <- file.path(output_dir, paste0(sid, "_per_chromosome_methylation", file_ext))

      tryCatch({
        ggplot2::ggsave(filename = filename, plot = p, width = 12, height = 7)
      }, error = function(e) {
        warning(paste("Could not save plot for", sid, ":", e$message))
      })
    }
    message("Plots saved to: ", output_dir)
  }

  return(all_stats_df)
}



#' Extract (and optionally plot) methylation stats for control genomes
#'
#' Scans samples for specific control chromosomes (e.g., pUC19, E. coli)
#' and returns a summary table of methylation levels, coverage, and position counts.
#'
#' @param methyl_obj A `methylRawList` object from methylKit.
#' @param target_chromosomes Character vector of chromosomes to look for. Default: c("pUC19", "NC_000913.3").
#' @param save_plot_to Path to save the plot (e.g. "figures/controls.png"). If NULL (default), no plot is saved.
#' @param sort_samples Logical. If TRUE, sorts sample names alphabetically in the output.
#'
#' @return A data.frame containing statistics for each sample/chromosome pair.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom purrr map_dfr
#' @importFrom methylKit getData
#' @export
mETHYLotest.NGS.get_control_stats <- function(methyl_obj,
                                              target_chromosomes = c("pUC19", "NC_000913.3"),
                                              save_plot_to = NULL,
                                              sort_samples = TRUE) {

  # Iterate over samples
  results_df <- purrr::map_dfr(methyl_obj, function(sample_obj) {

    sample_id <- sample_obj@sample.id
    raw_data  <- methylKit::getData(sample_obj)

    # Check each target chromosome for this sample
    purrr::map_dfr(target_chromosomes, function(chrom) {

      # Fast subsetting (Base R is faster here than dplyr filter)
      chr_data <- raw_data[raw_data$chr == chrom, ]

      if (nrow(chr_data) > 0) {
        total_Cs <- sum(chr_data$numCs)
        total_Ts <- sum(chr_data$numTs)
        total_bases <- total_Cs + total_Ts

        perc_meth <- if (total_bases > 0) (total_Cs / total_bases) * 100 else 0
        num_positions <- nrow(chr_data)
        num_reads     <- sum(chr_data$coverage)

      } else {
        perc_meth     <- NA_real_
        num_positions <- 0
        num_reads     <- 0
      }

      data.frame(
        Sample = sample_id,
        Chromosome = chrom,
        MethylationPercentage = perc_meth,
        NumberOfPositions = num_positions,
        NumberOfReads = num_reads,
        stringsAsFactors = FALSE
      )
    })
  })

  # Handle sorting
  if (sort_samples && nrow(results_df) > 0) {
    sorted_levels <- sort(unique(results_df$Sample))
    results_df$Sample <- factor(results_df$Sample, levels = sorted_levels)
  }

  # 2. Plotting logic
  if (!is.null(save_plot_to) && nrow(results_df) > 0) {

    message(paste("Saving control stats plot to:", save_plot_to))

    out_dir <- dirname(save_plot_to)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Dynamic width
    n_samples  <- length(unique(results_df$Sample))
    calc_width <- max(8, n_samples * 0.8)

    p <- ggplot2::ggplot(results_df, ggplot2::aes(x = .data$Sample, y = .data$MethylationPercentage, fill = .data$Chromosome)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9), color = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = round(.data$MethylationPercentage, 2)),
        position = ggplot2::position_dodge(width = 0.9),
        vjust = -0.5,
        size = 3
      ) +
      ggplot2::labs(
        title = "Methylation Percentage: Control Genomes",
        x = "Sample",
        y = "Methylation %"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top"
      )

    ggplot2::ggsave(filename = save_plot_to, plot = p, width = calc_width, height = 6)
  }

  return(results_df)
}



#' Calculate and plot number of analyzed positions (Individual vs Common)
#'
#' Counts the number of genomic positions (CpGs) covered/analyzed for each sample.
#' Also calculates the intersection (positions covered in ALL samples).
#'
#' @param methyl_obj A methylRawList object
#' @param output_dir Path to save the plot (optional)
#' @param file_name Name of the output plot file
#' @param unite_destrand Logical, passed to methylKit::unite (merges strands if TRUE)
#' @param min_cov Integer, for display in title only
#'
#' @return A dataframe with counts per sample and common intersection
#' @importFrom methylKit getData unite
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual labs theme_minimal theme element_text ggsave
#' @importFrom dplyr bind_rows
#' @export
mETHYLotest.NGS.get_position_counts_with_common <- function(methyl_obj,
                                                            output_dir = NULL,
                                                            file_name = "positions_count_comparison.png",
                                                            unite_destrand = FALSE,
                                                            min_cov = NULL) {

  message("Counting analyzed positions per sample...")

  # --- 1. S4-Safe Individual counts ---
  # Conversion explicite en liste pour itérer proprement sur l'objet S4
  methyl_list <- as.list(methyl_obj)

  if (length(methyl_list) == 0) {
    warning("The methyl_obj provided is empty. Cannot generate stats.")
    return(NULL)
  }

  list_of_dfs <- lapply(methyl_list, function(x) {

    # [FIX IMPORTANT] On récupère la dataframe interne pour être sûr du compte
    # nrow(x) direct peut parfois échouer sur des objets S4 vides ou mal formés
    df_internal <- methylKit::getData(x)
    n_rows <- nrow(df_internal)

    # Sécurisation de l'ID échantillon
    s_id <- as.character(x@sample.id)
    if (length(s_id) == 0) s_id <- "Unknown"

    # Debug info pour vérifier si les données sont bien lues
    if (n_rows == 0) {
      message(sprintf("  Warning: Sample '%s' has 0 positions.", s_id))
    } else {
      # message(sprintf("  Sample '%s': %d positions.", s_id, n_rows)) # Décommenter pour verbeux
    }

    data.frame(
      Sample = s_id,
      PositionCount = n_rows,
      Type = "Individual",
      stringsAsFactors = FALSE
    )
  })

  # Agrégation des résultats individuels
  stats_df <- dplyr::bind_rows(list_of_dfs)
  final_df <- stats_df
  intersection_exists <- FALSE

  # --- 2. Calculate Common Positions (Intersection) ---
  if (length(methyl_obj) > 1) {
    message("Calculating intersection (sites covered in all samples)...")

    # tryCatch obligatoire : unite() renvoie une erreur fatale s'il n'y a aucune intersection
    united_obj <- tryCatch({
      methylKit::unite(methyl_obj, destrand = unite_destrand)
    }, error = function(e) {
      warning("Intersection calculation failed (likely no common bases found).")
      return(NULL)
    })

    if (!is.null(united_obj)) {
      # Là aussi, on utilise getData() pour être sûr
      df_united <- methylKit::getData(united_obj)
      common_count <- nrow(df_united)

      if(common_count > 0) {
        intersection_exists <- TRUE
        common_row <- data.frame(
          Sample = "Common Intersection",
          PositionCount = common_count,
          Type = "Common",
          stringsAsFactors = FALSE
        )
        final_df <- dplyr::bind_rows(stats_df, common_row)
        message(sprintf("  Intersection found: %d common positions.", common_count))
      } else {
        message("  Intersection result is empty (0 common positions).")
      }
    }
  }

  # --- 3. Plotting Logic ---
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    message("Generating comparison plot...")

    # Gestion des facteurs pour l'ordre d'affichage (Intersection toujours à la fin)
    sample_levels <- sort(unique(stats_df$Sample))
    all_levels <- if(intersection_exists) c(sample_levels, "Common Intersection") else sample_levels
    final_df$Sample <- factor(final_df$Sample, levels = all_levels)

    # Titre dynamique
    plot_title <- "Number of Analyzed Positions"
    if (!is.null(min_cov)) plot_title <- paste0(plot_title, " (Coverage >= ", min_cov, ")")

    p <- ggplot2::ggplot(final_df, ggplot2::aes(x = .data$Sample, y = .data$PositionCount, fill = .data$Type)) +
      ggplot2::geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      ggplot2::geom_text(ggplot2::aes(label = format(.data$PositionCount, big.mark = ",")),
                         vjust = -0.5, size = 3) +
      ggplot2::scale_fill_manual(values = c("Individual" = "steelblue", "Common" = "firebrick")) +
      ggplot2::labs(
        title = plot_title,
        subtitle = paste("Total Samples:", length(methyl_obj)),
        x = "Sample / Set",
        y = "Count of CpG Sites"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none"
      )

    # Calcul de la largeur de l'image
    calc_width <- max(6, length(all_levels) * 0.7)

    full_path <- file.path(output_dir, file_name)
    ggplot2::ggsave(filename = full_path, plot = p, width = calc_width, height = 6)
    message("Plot saved to: ", full_path)
  }

  return(final_df)
}



#' Calculate, plot and export total C vs T counts per sample
#'
#' Aggregates the total number of Methylated Cytosines (numCs) and
#' Unmethylated Cytosines (numTs) for each sample.
#' Generates a stacked bar plot and an Excel report with percentages.
#'
#' @param methyl_obj A `methylRawList` object.
#' @param output_dir Directory path to save outputs. If NULL, nothing is saved.
#' @param plot_name Name of the output plot file (e.g., "total_C_counts.png").
#' @param excel_name Name of the output Excel file (e.g., "total_C_counts.xlsx").
#'                   If NULL, no Excel file is generated.
#' @param scale_y_log10 Logical. If TRUE, uses a log10 scale for the plot.
#'
#' @return A data.frame (Long format) used for plotting.
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import writexl
#' @importFrom purrr map_dfr
#' @importFrom methylKit getData
#' @export
mETHYLotest.NGS.get_C_counts <- function(methyl_obj,
                                         output_dir = NULL,
                                         plot_name = "total_C_counts.png",
                                         excel_name = "total_C_counts.xlsx",
                                         scale_y_log10 = FALSE) {

  message("Aggregating total C (Methylated) and T (Unmethylated) counts...")

  # 1. Calculate sums per sample (Long Format)
  stats_df <- purrr::map_dfr(methyl_obj, function(x) {

    df <- methylKit::getData(x)

    if(nrow(df) == 0) return(data.frame())

    total_meth   <- sum(as.numeric(df$numCs))
    total_unmeth <- sum(as.numeric(df$numTs))

    rbind(
      data.frame(Sample = x@sample.id, Type = "Methylated Cs", Count = total_meth, stringsAsFactors = FALSE),
      data.frame(Sample = x@sample.id, Type = "Unmethylated Cs", Count = total_unmeth, stringsAsFactors = FALSE)
    )
  })

  if (nrow(stats_df) == 0) {
    warning("No data found in methyl_obj.")
    return(stats_df)
  }

  # Ensure output directory exists if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 2. Excel Export Logic (Wide Format calculation)
  if (!is.null(output_dir) && !is.null(excel_name)) {
    message("Generating Excel report...")

    # Transformation en format large : une ligne par sample
    excel_df <- stats_df %>%
      tidyr::pivot_wider(names_from = Type, values_from = Count) %>%
      dplyr::mutate(
        Total_Reads = `Methylated Cs` + `Unmethylated Cs`,
        # Ajout d'une condition pour éviter division par zéro
        Global_Methylation_Percent = ifelse(Total_Reads > 0, (`Methylated Cs` / Total_Reads) * 100, 0)
      ) %>%
      dplyr::arrange(Sample)

    excel_path <- file.path(output_dir, excel_name)
    writexl::write_xlsx(excel_df, path = excel_path)
    message("Excel saved to: ", excel_path)
  }

  # 3. Plotting Logic
  if (!is.null(output_dir) && !is.null(plot_name)) {

    message("Generating plot...")

    # Ordering
    sample_levels <- sort(unique(stats_df$Sample))
    stats_df$Sample <- factor(stats_df$Sample, levels = sample_levels)

    my_colors <- c("Methylated Cs" = "#E31A1C", "Unmethylated Cs" = "#1F78B4")

    # CORRECTION ICI : Remplacement de size par linewidth
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = .data$Sample, y = .data$Count, fill = .data$Type)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = my_colors) +
      ggplot2::labs(
        title = "Total Cytosine Calls per Sample",
        subtitle = "Stacked: Unmethylated (Ts) + Methylated (Cs)",
        x = "Sample",
        y = "Total Count (Reads)",
        fill = "Call Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "top"
      )

    if (scale_y_log10) {
      p <- p + ggplot2::scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
    } else {
      p <- p + ggplot2::scale_y_continuous(labels = scales::label_comma())
    }

    n_samples <- length(unique(stats_df$Sample))
    calc_width <- max(6, n_samples * 0.7)

    plot_path <- file.path(output_dir, plot_name)
    ggplot2::ggsave(filename = plot_path, plot = p, width = calc_width, height = 7)
    message("Plot saved to: ", plot_path)
  }

  return(stats_df)
}

#' Apply coordinate offset to methylRawList
#'
#' @description
#' Adds a fixed offset to start and end coordinates of all samples
#' in a methylRawList object. Used to convert 0-based BED coordinates
#' to 1-based coordinates expected by methylKit.
#'
#' @param myobj A \code{methylRawList} object from \code{methylKit::methRead}.
#' @param offset Integer. Value to add to start and end (default: 1).
#'
#' @return A corrected \code{methylRawList} object.
#' @export
#' @importFrom methylKit getData
mETHYLotest.NGS.ApplyOffsetToObj <- function(myobj, offset = 1L) {

  if (offset == 0L) {
    message("[mETHYLotest] Offset is 0, skipping.")
    return(myobj)
  }

  correct_one <- function(obj) {
    obj[["start"]] <- obj[["start"]] + offset
    obj[["end"]]   <- obj[["end"]]   + offset
    obj
  }

  message(sprintf("[mETHYLotest] Applying +%d coordinate offset...",
                  offset))
  corrected <- lapply(myobj, correct_one)

  new("methylRawList",
      corrected,
      treatment = myobj@treatment)
}
