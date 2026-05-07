#' Quick Common Positions Analysis: Controls as Reference
#'
#' Fast approach: builds the common CpG set from Control samples,
#' then checks how well each Test sample covers it.
#' Produces a single summary figure.
#'
#' @param methyl_obj A methylRawList object (filtered).
#' @param output_file Path for the output plot (PNG).
#' @param treatment_labels Named vector: 0 = "Control", 1 = "Test" (optional).
#'
#' @return A data.frame with coverage stats per sample (invisibly).
#' @export
mETHYLotest.NGS.QuickCommonPositions <- function(methyl_obj,
                                                    output_file = NULL,
                                                    treatment_labels = c("0" = "Control", "1" = "Test")) {

  stopifnot(inherits(methyl_obj, "methylRawList"))

  sample_ids <- methylKit::getSampleID(methyl_obj)
  treatments <- methyl_obj@treatment
  n <- length(sample_ids)

  # Labels
  group_labels <- treatment_labels[as.character(treatments)]
  names(group_labels) <- sample_ids

  ctrl_ids <- sample_ids[treatments == 0]
  test_ids <- sample_ids[treatments == 1]

  message(sprintf("[mETHYLotest] %d Controls, %d Tests", length(ctrl_ids), length(test_ids)))

  # ══════════════════════════════════════════════════════
  # 1. Extract position sets (fast: just chr:start strings)
  # ══════════════════════════════════════════════════════
  message("[mETHYLotest] Extracting positions...")

  pos_list <- lapply(seq_len(n), function(i) {
    d <- methylKit::getData(methyl_obj[[i]])
    paste0(d$chr, ":", d$start)
  })
  names(pos_list) <- sample_ids

  n_per_sample <- sapply(pos_list, length)

  # ══════════════════════════════════════════════════════
  # 2. Common positions among ALL Controls
  # ══════════════════════════════════════════════════════
  message("[mETHYLotest] Computing Control common set...")

  ctrl_common <- pos_list[[ctrl_ids[1]]]
  for (cid in ctrl_ids[-1]) {
    ctrl_common <- intersect(ctrl_common, pos_list[[cid]])
  }
  n_ctrl_common <- length(ctrl_common)

  message(sprintf("[mETHYLotest]   Control common CpGs: %s",
                  format(n_ctrl_common, big.mark = ",")))

  # ══════════════════════════════════════════════════════
  # 3. For EVERY sample: overlap with control common set
  # ══════════════════════════════════════════════════════
  message("[mETHYLotest] Computing overlaps...")

  stats <- data.frame(
    Sample       = sample_ids,
    Group        = unname(group_labels),
    Total_CpGs   = n_per_sample,
    In_Ctrl_Common = sapply(pos_list, function(p) sum(p %in% ctrl_common)),
    stringsAsFactors = FALSE
  )
  stats$Pct_Coverage <- round(100 * stats$In_Ctrl_Common / max(n_ctrl_common, 1), 1)
  stats$Missing      <- n_ctrl_common - stats$In_Ctrl_Common

  # Sort: Controls first, then Tests by coverage
  stats <- stats[order(stats$Group != "Control", -stats$Pct_Coverage), ]

  # ══════════════════════════════════════════════════════
  # 4. What happens if we add each Test sample?
  # ══════════════════════════════════════════════════════
  # Common across all Controls + this Test sample
  stats$Common_With_Ctrls <- sapply(stats$Sample, function(sid) {
    if (group_labels[sid] == "Control") return(NA_integer_)
    length(intersect(ctrl_common, pos_list[[sid]]))
  })

  # ══════════════════════════════════════════════════════
  # 5. Single combined figure
  # ══════════════════════════════════════════════════════
  message("[mETHYLotest] Generating plot...")

  library(ggplot2)

  stats$Sample <- factor(stats$Sample, levels = rev(stats$Sample))

  col_ctrl <- "#3498db"
  col_test <- "#e74c3c"
  stats$Fill <- ifelse(stats$Group == "Control", col_ctrl, col_test)

  # ── Panel A: total CpGs per sample ──
  pA <- ggplot(stats, aes(Sample, Total_CpGs, fill = Group)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = c(Control = col_ctrl, Test = col_test)) +
    geom_hline(yintercept = n_ctrl_common, linetype = "dashed",
               colour = "black", linewidth = 0.6) +
    annotate("text", x = nrow(stats), y = n_ctrl_common,
             label = sprintf("Ctrl common: %s", format(n_ctrl_common, big.mark = ",")),
             vjust = -0.5, hjust = 0, size = 3) +
    coord_flip() +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(title = "A. Total CpG positions per sample",
         x = NULL, y = "CpG positions") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  # ── Panel B: % coverage of control common set ──
  pB <- ggplot(stats, aes(Sample, Pct_Coverage, fill = Group)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = c(Control = col_ctrl, Test = col_test)) +
    geom_hline(yintercept = 90, linetype = "dotted", colour = "orange") +
    geom_text(aes(label = sprintf("%.1f%%", Pct_Coverage)),
              hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 105)) +
    labs(title = "B. Coverage of Control common CpGs",
         subtitle = sprintf("Reference: %s CpGs common across all %d Controls",
                            format(n_ctrl_common, big.mark = ","), length(ctrl_ids)),
         x = NULL, y = "% of Control common CpGs covered") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  # ── Combine ──
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_final <- pA + pB +
      patchwork::plot_layout(ncol = 2, widths = c(1, 1.2)) +
      patchwork::plot_annotation(
        title    = "Common Positions Analysis",
        subtitle = sprintf("%d Controls (blue) | %d Tests (red)",
                           length(ctrl_ids), length(test_ids)),
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
  } else {
    # Fallback: just panel B (the most informative)
    p_final <- pB
  }

  # ── Save ──
  if (!is.null(output_file)) {
    d <- dirname(output_file)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)

    h <- max(5, 0.45 * n)
    ggsave(output_file, p_final, width = 14, height = h, dpi = 150)
    message("[mETHYLotest] Saved: ", output_file)

    # PDF too
    pdf_file <- sub("\\.png$", ".pdf", output_file)
    ggsave(pdf_file, p_final, width = 14, height = h)
  } else {
    print(p_final)
  }

  # ══════════════════════════════════════════════════════
  # 6. Print summary
  # ══════════════════════════════════════════════════════
  message("\n[mETHYLotest] === Summary ===")
  message(sprintf("  Control common CpGs : %s", format(n_ctrl_common, big.mark = ",")))

  test_stats <- stats[stats$Group == "Test", ]
  if (nrow(test_stats) > 0) {
    message(sprintf("  Test coverage range : %.1f%% - %.1f%%",
                    min(test_stats$Pct_Coverage), max(test_stats$Pct_Coverage)))

    low <- test_stats$Sample[test_stats$Pct_Coverage < 90]
    if (length(low) > 0) {
      message(sprintf("  WARNING - Test samples < 90%% coverage: %s",
                      paste(low, collapse = ", ")))
    } else {
      message("  All Test samples cover > 90% of Control common CpGs.")
    }
  }

  return(invisible(stats))
}
