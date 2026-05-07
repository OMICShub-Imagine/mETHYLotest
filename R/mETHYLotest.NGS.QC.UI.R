#' Interactive QC Dashboard (Iterative) - Refined UI
#'
#' @description
#' Launches a Shiny Dashboard for QC visualization and parameter refinement.
#' Includes controls for Coverage, Percentiles, Samples, and Chromosomes filtering.
#'
#' @param df_meth Dataframe. Chromosome methylation stats.
#' @param df_controls Dataframe. Control genome stats.
#' @param df_stats Dataframe. Global coverage/position stats.
#' @param df_global_meth Dataframe. Global methylation per sample (from QC function).
#' @param current_cov Integer. Current minimum coverage used (lo.count).
#' @param current_hi_perc Numeric. Current High Percentile filter (default 99.9).
#' @param current_lo_perc Numeric. Current Low Percentile filter (default NULL).
#' @param current_excluded Vector. List of already excluded sample IDs.
#' @param all_chrs Vector. List of all available chromosomes in the dataset.
#' @param current_chrs_kept Vector. List of currently kept chromosomes.
#' @param sample_treatments Named vector. Treatment per sample (names = sample IDs, values = 0/1).
#'
#' @return A list containing action="update"|"proceed" and all filtering parameters.
#' @export
#' @import shiny
#' @import shinydashboard
#' @import plotly
#' @import DT
#' @import ggplot2
mETHYLotest.NGS.QC.UI <- function(df_meth = NULL,
                                  df_controls = NULL,
                                  df_stats = NULL,
                                  df_global_meth = NULL,
                                  current_cov = 1,
                                  current_hi_perc = 99.9,
                                  current_lo_perc = NULL,
                                  current_excluded = NULL,
                                  all_chrs = NULL,
                                  current_chrs_kept = NULL,
                                  sample_treatments = NULL) {

  requireNamespace("shiny")
  requireNamespace("shinydashboard")
  requireNamespace("plotly")
  requireNamespace("DT")

  # --- PRE-PROCESSING ---
  available_samples <- character(0)
  if (!is.null(df_stats) && "Sample" %in% colnames(df_stats)) {
    available_samples <- unique(as.character(df_stats$Sample[df_stats$Type != "Common"]))
  } else if (!is.null(df_meth)) {
    available_samples <- unique(as.character(df_meth$Sample))
  }

  if (is.null(all_chrs) && !is.null(df_meth)) {
    all_chrs <- unique(df_meth$chr)
  }
  if (is.null(current_chrs_kept)) {
    current_chrs_kept <- all_chrs
  }

  # ── Build treatment lookup ──
  get_group_label <- function(sample_id) {
    if (is.null(sample_treatments)) return("Unknown")
    sid <- as.character(sample_id)
    val <- sample_treatments[sid]
    if (length(val) == 0 || is.na(val)) return("Unknown")
    if (val == 0) return("Control (0)")
    if (val == 1) return("Case (1)")
    return(paste0("Group ", val))
  }

  # ── Enrich dataframes with treatment group ──
  if (!is.null(df_stats) && !is.null(sample_treatments)) {
    df_stats$Group <- vapply(as.character(df_stats$Sample),
                             get_group_label, character(1))
  }
  if (!is.null(df_meth) && !is.null(sample_treatments)) {
    df_meth$Group <- vapply(as.character(df_meth$Sample),
                            get_group_label, character(1))
  }
  if (!is.null(df_controls) && !is.null(sample_treatments)) {
    df_controls$Group <- vapply(as.character(df_controls$Sample),
                                get_group_label, character(1))
  }

  # ── Treatment summary (active samples only) ──
  if (!is.null(sample_treatments)) {
    active_tx <- sample_treatments[names(sample_treatments) %in% available_samples]
    n_ctrl <- sum(active_tx == 0, na.rm = TRUE)
    n_case <- sum(active_tx == 1, na.rm = TRUE)
  } else {
    n_ctrl <- 0L
    n_case <- 0L
  }

  # ── Global methylation summary per group ──
  ctrl_meth <- NA; case_meth <- NA
  ctrl_sd   <- NA; case_sd   <- NA

  if (!is.null(df_global_meth) && "Group" %in% colnames(df_global_meth)) {
    ctrl_rows <- df_global_meth[df_global_meth$Group == "Control (0)", ]
    case_rows <- df_global_meth[df_global_meth$Group == "Case (1)", ]
    if (nrow(ctrl_rows) > 0) {
      ctrl_meth <- round(mean(ctrl_rows$Global_Meth_Pct, na.rm = TRUE), 1)
      ctrl_sd   <- round(sd(ctrl_rows$Global_Meth_Pct, na.rm = TRUE), 1)
    }
    if (nrow(case_rows) > 0) {
      case_meth <- round(mean(case_rows$Global_Meth_Pct, na.rm = TRUE), 1)
      case_sd   <- round(sd(case_rows$Global_Meth_Pct, na.rm = TRUE), 1)
    }
  }

  # --- UI DEFINITION ---
  ui <- shinydashboard::dashboardPage(
    skin = "blue",

    shinydashboard::dashboardHeader(title = "mETHYLotest QC Dashboard"),

    shinydashboard::dashboardSidebar(
      width = 300,

      div(style = "padding: 15px;",
          h5("Current Parameters", style = "color: #b8c7ce;"),
          div(class = "callout callout-info",
              style = "padding: 10px; border-left-color: #00c0ef; background-color: #222d32 !important; border: 1px solid #444;",
              p(strong("Min Cov: "), current_cov, "x"),
              p(strong("Hi. Perc: "), current_hi_perc, "%"),
              p(strong("Samples: "), length(available_samples)),
              if (!is.null(sample_treatments)) {
                tagList(
                  p(strong("Controls: "),
                    span(n_ctrl, style = "color: #27ae60; font-weight: bold;"),
                    if (!is.na(ctrl_meth))
                      span(sprintf(" (%.1f%%)", ctrl_meth),
                           style = "color: #27ae60; font-size: 0.9em;")),
                  p(strong("Cases: "),
                    span(n_case, style = "color: #e74c3c; font-weight: bold;"),
                    if (!is.na(case_meth))
                      span(sprintf(" (%.1f%%)", case_meth),
                           style = "color: #e74c3c; font-size: 0.9em;")),
                  if (!is.na(ctrl_meth) && !is.na(case_meth)) {
                    delta <- round(case_meth - ctrl_meth, 1)
                    delta_col <- if (abs(delta) > 5) "#e74c3c" else "#95a5a6"
                    p(strong("Delta: "),
                      span(sprintf("%+.1f%%", delta),
                           style = sprintf("color:%s;font-weight:bold;", delta_col)))
                  }
                )
              },
              p(strong("Chromosomes: "), length(current_chrs_kept))
          ),
          if (!is.null(current_excluded) && length(current_excluded) > 0) {
            div(style = "margin-top: 10px; color: #ffcccc; font-size: 0.9em;",
                icon("trash"),
                paste(length(current_excluded), "excluded:",
                      paste(current_excluded, collapse = ", "))
            )
          }
      ),

      shinydashboard::sidebarMenu(
        id = "tabs",
        shinydashboard::menuItem("1. Visualisation", tabName = "plots",
                                 icon = icon("chart-area")),
        shinydashboard::menuItem("2. Data Tables", tabName = "tables",
                                 icon = icon("table")),
        shinydashboard::menuItem("3. Refine & Validate", tabName = "actions",
                                 icon = icon("sliders-h"),
                                 badgeLabel = "Action", badgeColor = "orange")
      )
    ),

    shinydashboard::dashboardBody(
      tags$head(tags$style(HTML("
        .group-legend {
          display: inline-block; padding: 2px 8px; border-radius: 3px;
          font-size: 0.85em; font-weight: bold; margin: 2px;
        }
        .group-ctrl { background: #d5f5e3; color: #1e8449; border: 1px solid #27ae60; }
        .group-case { background: #fadbd8; color: #922b21; border: 1px solid #e74c3c; }
      "))),

      shinydashboard::tabItems(

        # --- TAB 1: VISUALISATION ---
        shinydashboard::tabItem(
          tabName = "plots",

          if (!is.null(df_stats)) {
            fluidRow(
              shinydashboard::infoBoxOutput("total_samples_box", width = 3),
              shinydashboard::infoBoxOutput("common_sites_box", width = 3),
              if (!is.null(sample_treatments)) {
                tagList(
                  shinydashboard::infoBoxOutput("group_balance_box", width = 3),
                  shinydashboard::infoBoxOutput("group_meth_box", width = 3)
                )
              }
            )
          },

          # ── Group legend ──
          if (!is.null(sample_treatments)) {
            fluidRow(
              column(12,
                     div(style = "text-align: center; margin-bottom: 10px;",
                         span(class = "group-legend group-ctrl",
                              icon("circle"), " Control (0)"),
                         span(class = "group-legend group-case",
                              icon("circle"), " Case (1)")
                     )
              )
            )
          },

          fluidRow(
            shinydashboard::tabBox(
              title = tagList(icon("chart-bar"), " QC Metrics"),
              id = "tabset_plots", width = 12,

              tabPanel("CpG Positions",
                       if (!is.null(df_stats)) {
                         tagList(
                           div(class = "alert alert-info",
                               style = "font-size: 12px; padding: 8px; margin: 5px 0;",
                               icon("info-circle"),
                               " This plot shows the number of ",
                               strong("CpG positions passing coverage filters"),
                               " per sample. ",
                               "'Individual' = positions in that sample alone; ",
                               "'Common' = positions shared across ALL samples."),
                           plotly::plotlyOutput("plot_stats", height = "550px")
                         )
                       } else "No data."
              ),

              tabPanel("Chromosomes",
                       if (!is.null(df_meth)) {
                         fluidRow(
                           column(3,
                                  selectInput("chr_sample_select",
                                              "Select Sample:",
                                              choices = unique(df_meth$Sample))),
                           column(9,
                                  plotly::plotlyOutput("plot_chr",
                                                       height = "500px")))
                       } else "No data."
              ),

              tabPanel("Controls",
                       if (!is.null(df_controls)) {
                         plotly::plotlyOutput("plot_controls",
                                              height = "550px")
                       } else "No controls."
              ),

              tabPanel("Sample Overview",
                       if (!is.null(df_global_meth)) {
                         tagList(
                           div(class = "alert alert-info",
                               style = "font-size: 12px; padding: 8px;",
                               icon("info-circle"),
                               " Global methylation per sample, colored by treatment group."),
                           plotly::plotlyOutput("plot_sample_overview",
                                                height = "400px"),
                           if ("Group" %in% colnames(df_global_meth)) {
                             tagList(
                               hr(),
                               plotly::plotlyOutput("plot_group_comparison",
                                                    height = "400px")
                             )
                           }
                         )
                       } else if (!is.null(df_meth) && !is.null(sample_treatments)) {
                         tagList(
                           div(class = "alert alert-info",
                               style = "font-size: 12px; padding: 8px;",
                               icon("info-circle"),
                               " Average methylation per sample, colored by treatment group."),
                           plotly::plotlyOutput("plot_sample_overview_fallback",
                                                height = "500px")
                         )
                       } else "No data or treatment info."
              )
            )
          )
        ),

        # --- TAB 2: DATA TABLES ---
        shinydashboard::tabItem(
          tabName = "tables",
          shinydashboard::box(
            title = "Raw Data", status = "primary",
            solidHeader = TRUE, width = 12,
            tabsetPanel(
              if (!is.null(df_stats))
                tabPanel("Global Stats", br(),
                         DT::DTOutput("table_stats")) else NULL,
              if (!is.null(df_meth))
                tabPanel("Chromosome Stats", br(),
                         DT::DTOutput("table_meth")) else NULL,
              if (!is.null(df_controls))
                tabPanel("Control Stats", br(),
                         DT::DTOutput("table_controls")) else NULL,
              if (!is.null(df_global_meth))
                tabPanel("Global Methylation", br(),
                         DT::DTOutput("table_global_meth")) else NULL
            )
          )
        ),

        # --- TAB 3: REFINE & VALIDATE ---
        shinydashboard::tabItem(
          tabName = "actions",
          fluidRow(
            shinydashboard::box(
              title = "Filter Settings", status = "warning",
              solidHeader = TRUE, width = 6,

              h4(icon("filter"), " Coverage Filtering"),

              sliderInput("new_min_cov",
                          "Min Coverage Count (lo.count):",
                          min = 1, max = 50,
                          value = current_cov, step = 1),

              h5(strong("Percentile Filtering (Advanced)")),
              splitLayout(
                numericInput("new_hi_perc", "High % (hi.perc):",
                             value = current_hi_perc,
                             min = 80, max = 100, step = 0.1),
                numericInput("new_lo_perc", "Low % (lo.perc):",
                             value = current_lo_perc,
                             min = 0, max = 20, step = 0.1)
              ),
              helpText("Use Hi % (e.g., 99.9) to remove PCR artifacts."),

              hr(),

              h4(icon("dna"), " Chromosome Filtering"),
              selectInput("chrs_to_keep",
                          "Select Chromosomes to KEEP:",
                          choices = all_chrs,
                          selected = current_chrs_kept,
                          multiple = TRUE),
              helpText("Remove sex chromosomes (chrX, chrY) or",
                       " unplaced contigs if needed."),

              hr(),

              h4(icon("trash-alt"), " Sample Exclusion"),
              if (length(available_samples) > 0) {
                sample_choices <- available_samples
                if (!is.null(sample_treatments)) {
                  sample_labels <- vapply(available_samples, function(s) {
                    g <- get_group_label(s)
                    paste0(s, " [", g, "]")
                  }, character(1))
                  names(sample_choices) <- sample_labels
                }
                selectInput("samples_to_remove",
                            "Select samples to EXCLUDE:",
                            choices = sample_choices,
                            multiple = TRUE)
              } else {
                p("No samples left.", style = "color:red")
              }
            ),

            shinydashboard::box(
              title = "Validation", status = "success",
              solidHeader = TRUE, width = 6,
              h4("Next Step"),
              br(),
              div(style = paste0("text-align: center; background-color: #f9f9f9;",
                                 " padding: 15px; border-radius: 5px;",
                                 " border: 1px solid #ddd;"),
                  strong("Option A: Apply Changes"),
                  p("Re-calculate QC metrics with these new settings."),
                  actionButton("btn_update", "Update / Re-Filter",
                               icon = icon("sync"),
                               class = "btn-warning btn-lg",
                               style = "width: 80%;")
              ),
              br(),
              div(style = paste0("text-align: center; background-color: #f9f9f9;",
                                 " padding: 15px; border-radius: 5px;",
                                 " border: 1px solid #ddd;"),
                  strong("Option B: Finish QC"),
                  p("Accept settings and proceed to analysis."),
                  actionButton("btn_proceed", "Validate & Proceed",
                               icon = icon("check-circle"),
                               class = "btn-success btn-lg",
                               style = "width: 80%;")
              )
            )
          ),

          if (!is.null(current_excluded) && length(current_excluded) > 0) {
            fluidRow(shinydashboard::box(
              title = "History of Exclusions", status = "danger",
              width = 12, collapsible = TRUE, collapsed = TRUE,
              tags$ul(lapply(current_excluded, function(ex) {
                g <- get_group_label(ex)
                tags$li(strong(ex),
                        span(paste0(" [", g, "]"),
                             style = if (grepl("Control", g))
                               "color: #27ae60;"
                             else if (grepl("Case", g))
                               "color: #e74c3c;"
                             else "color: #999;"))
              }))
            ))
          }
        )
      )
    )
  )

  # --- SERVER ---
  server <- function(input, output, session) {

    # ── Color palette for groups ──
    group_colors <- c("Control (0)" = "#27ae60",
                      "Case (1)"    = "#e74c3c",
                      "Unknown"     = "#95a5a6")

    # ══════════════════════════════════════════════════════════════
    # TAB 1: PLOTS
    # ══════════════════════════════════════════════════════════════

    if (!is.null(df_stats)) {

      output$plot_stats <- plotly::renderPlotly({
        plot_df <- df_stats

        if ("Group" %in% colnames(plot_df)) {
          indiv_df <- plot_df[plot_df$Type == "Individual", ]

          p <- ggplot2::ggplot(
            indiv_df,
            ggplot2::aes(x = reorder(Sample, -PositionCount),
                         y = PositionCount, fill = Group,
                         text = paste0("Sample: ", Sample,
                                       "\nGroup: ", Group,
                                       "\nCpG positions: ",
                                       format(PositionCount, big.mark = ",")))) +
            ggplot2::geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
            ggplot2::scale_fill_manual(values = group_colors) +
            ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
            ggplot2::labs(title = "CpG Positions per Sample (after coverage filter)",
                          y = "Number of CpG positions passing filters",
                          x = "", fill = "Group") +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

          common_val <- plot_df$PositionCount[plot_df$Type == "Common"]
          if (length(common_val) > 0 && common_val[1] > 0) {
            p <- p + ggplot2::geom_hline(yintercept = common_val[1],
                                         linetype = "dashed", color = "firebrick", linewidth = 0.7) +
              ggplot2::annotate("text", x = 1, y = common_val[1],
                                label = paste0("Common: ", format(common_val[1], big.mark = ",")),
                                hjust = 0, vjust = -0.5, color = "firebrick", size = 3.5)
          }
          plotly::ggplotly(p, tooltip = "text")
        } else {
          p <- ggplot2::ggplot(plot_df,
                               ggplot2::aes(x = Sample, y = PositionCount, fill = Type)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
            ggplot2::scale_fill_manual(values = c("Individual" = "steelblue", "Common" = "firebrick")) +
            ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
            ggplot2::labs(y = "Number of CpG positions passing filters", x = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          plotly::ggplotly(p)
        }
      })

      output$table_stats <- DT::renderDT({
        DT::datatable(df_stats, options = list(pageLength = 15, scrollX = TRUE))
      })

      output$total_samples_box <- shinydashboard::renderInfoBox({
        n <- length(unique(df_stats$Sample[df_stats$Type != "Common"]))
        shinydashboard::infoBox("Active Samples", n,
                                icon = icon("vial"), color = "blue", fill = TRUE)
      })

      output$common_sites_box <- shinydashboard::renderInfoBox({
        val <- df_stats$PositionCount[df_stats$Type == "Common"]
        shinydashboard::infoBox("Common CpG Positions",
                                if (length(val) > 0) format(val, big.mark = ",") else "0",
                                icon = icon("compress-arrows-alt"), color = "red", fill = TRUE)
      })
    }

    if (!is.null(sample_treatments)) {
      output$group_balance_box <- shinydashboard::renderInfoBox({
        shinydashboard::infoBox("Group Balance",
                                paste0(n_ctrl, " ctrl / ", n_case, " case"),
                                icon = icon("balance-scale"),
                                color = if (n_ctrl == 0 || n_case == 0) "red" else "green",
                                fill = TRUE)
      })

      output$group_meth_box <- shinydashboard::renderInfoBox({
        if (!is.na(ctrl_meth) && !is.na(case_meth)) {
          delta <- round(case_meth - ctrl_meth, 1)
          shinydashboard::infoBox(
            "Global Methylation",
            sprintf("CTL %.1f%% / Case %.1f%%", ctrl_meth, case_meth),
            subtitle = sprintf("Delta: %+.1f%%", delta),
            icon = icon("dna"),
            color = if (abs(delta) > 10) "red"
            else if (abs(delta) > 5) "yellow" else "blue",
            fill = TRUE)
        } else {
          shinydashboard::infoBox("Global Methylation", "N/A",
                                  icon = icon("dna"), color = "gray", fill = TRUE)
        }
      })
    }

    # ── Chromosome plot ──
    if (!is.null(df_meth)) {
      filtered_meth <- shiny::reactive({
        shiny::req(input$chr_sample_select)
        df_meth[df_meth$Sample == input$chr_sample_select, ]
      })

      output$plot_chr <- plotly::renderPlotly({
        d <- filtered_meth()
        max_y <- if (nrow(d) > 0) max(df_meth$avg_methylation, na.rm = TRUE) * 1.05 else 100
        sample_group <- get_group_label(input$chr_sample_select)
        bar_color <- group_colors[sample_group]
        if (is.na(bar_color)) bar_color <- "#3498db"

        p <- ggplot2::ggplot(d,
                             ggplot2::aes(x = chr, y = avg_methylation,
                                          text = paste0("Chr: ", chr, "\nMeth: ",
                                                        round(avg_methylation, 1), "%"))) +
          ggplot2::geom_bar(stat = "identity", fill = bar_color, color = "white", linewidth = 0.3) +
          ggplot2::scale_y_continuous(limits = c(0, max_y)) +
          ggplot2::labs(title = paste0(input$chr_sample_select, " [", sample_group, "]"),
                        y = "Avg Methylation %", x = "Chromosome") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                         legend.position = "none")
        plotly::ggplotly(p, tooltip = "text")
      })

      output$table_meth <- DT::renderDT({
        DT::datatable(df_meth, options = list(pageLength = 15, scrollX = TRUE))
      })
    }

    # ── Controls plot ──
    if (!is.null(df_controls)) {
      output$plot_controls <- plotly::renderPlotly({
        if ("Group" %in% colnames(df_controls)) {
          df_controls$SampleLabel <- paste0(df_controls$Sample, " [",
                                            substr(df_controls$Group, 1, 4), "]")
          p <- ggplot2::ggplot(df_controls,
                               ggplot2::aes(x = SampleLabel, y = MethylationPercentage,
                                            fill = Chromosome,
                                            text = paste0("Sample: ", Sample,
                                                          "\nGroup: ", Group,
                                                          "\nControl: ", Chromosome,
                                                          "\nMeth: ", round(MethylationPercentage, 2), "%"))) +
            ggplot2::geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.3) +
            ggplot2::scale_fill_brewer(palette = "Set2") +
            ggplot2::geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 0.5) +
            ggplot2::geom_hline(yintercept = 95, linetype = "dashed", color = "darkgreen", linewidth = 0.5) +
            ggplot2::labs(y = "Methylation %", x = "", fill = "Control") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          plotly::ggplotly(p, tooltip = "text")
        } else {
          p <- ggplot2::ggplot(df_controls,
                               ggplot2::aes(x = Sample, y = MethylationPercentage, fill = Chromosome)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
            ggplot2::labs(y = "Methylation %", x = "Sample") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          plotly::ggplotly(p)
        }
      })

      output$table_controls <- DT::renderDT({
        DT::datatable(df_controls, options = list(pageLength = 15, scrollX = TRUE))
      })
    }

    # ── Global Methylation table ──
    if (!is.null(df_global_meth)) {
      output$table_global_meth <- DT::renderDT({
        DT::datatable(df_global_meth,
                      options = list(pageLength = 20, scrollX = TRUE)) |>
          DT::formatRound("Global_Meth_Pct", digits = 2)
      })
    }

    # ══════════════════════════════════════════════════════════════
    # Sample Overview: barplot per sample
    # ══════════════════════════════════════════════════════════════
    if (!is.null(df_global_meth)) {
      output$plot_sample_overview <- plotly::renderPlotly({
        gm <- df_global_meth

        if ("Group" %in% colnames(gm)) {
          p <- ggplot2::ggplot(
            gm,
            ggplot2::aes(x = reorder(Sample, Global_Meth_Pct),
                         y = Global_Meth_Pct, fill = Group,
                         text = paste0("Sample: ", Sample,
                                       "\nGroup: ", Group,
                                       "\nMeth: ", Global_Meth_Pct, "%",
                                       "\nPositions: ", format(N_Positions, big.mark = ",")))) +
            ggplot2::geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
            ggplot2::scale_fill_manual(values = group_colors) +
            ggplot2::coord_flip() +
            ggplot2::labs(title = "Global Methylation per Sample",
                          x = "", y = "Global Methylation (%)", fill = "Group") +
            ggplot2::theme_minimal(base_size = 12)
        } else {
          p <- ggplot2::ggplot(
            gm,
            ggplot2::aes(x = reorder(Sample, Global_Meth_Pct),
                         y = Global_Meth_Pct)) +
            ggplot2::geom_bar(stat = "identity", fill = "#3498db") +
            ggplot2::coord_flip() +
            ggplot2::labs(title = "Global Methylation per Sample",
                          x = "", y = "Global Methylation (%)") +
            ggplot2::theme_minimal(base_size = 12)
        }
        plotly::ggplotly(p, tooltip = "text")
      })

      # ── Group comparison: boxplot + points ──
      if ("Group" %in% colnames(df_global_meth)) {
        output$plot_group_comparison <- plotly::renderPlotly({
          gm <- df_global_meth[df_global_meth$Group != "Unknown", ]
          if (nrow(gm) == 0) return(plotly::plotly_empty())

          p <- ggplot2::ggplot(
            gm,
            ggplot2::aes(x = Group, y = Global_Meth_Pct, fill = Group,
                         text = paste0("Sample: ", Sample,
                                       "\nGroup: ", Group,
                                       "\nMeth: ", Global_Meth_Pct, "%"))) +
            ggplot2::geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
            ggplot2::geom_jitter(width = 0.15, size = 3, alpha = 0.8,
                                 ggplot2::aes(colour = Group)) +
            ggplot2::scale_fill_manual(values = group_colors) +
            ggplot2::scale_colour_manual(values = group_colors) +
            ggplot2::labs(
              title = "Global Methylation: Control vs Case",
              subtitle = if (!is.na(ctrl_meth) && !is.na(case_meth))
                sprintf("CTL: %.1f%% +/- %.1f | Case: %.1f%% +/- %.1f | Delta: %+.1f%%",
                        ctrl_meth, ctrl_sd, case_meth, case_sd,
                        case_meth - ctrl_meth) else "",
              x = "", y = "Global Methylation (%)", fill = "Group") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(legend.position = "none")

          plotly::ggplotly(p, tooltip = "text")
        })
      }
    }

    # ── Fallback: use df_meth if no df_global_meth ──
    if (is.null(df_global_meth) && !is.null(df_meth) && !is.null(sample_treatments)) {
      output$plot_sample_overview_fallback <- plotly::renderPlotly({
        avg_df <- stats::aggregate(avg_methylation ~ Sample,
                                   data = df_meth, FUN = mean, na.rm = TRUE)
        avg_df$Group <- vapply(as.character(avg_df$Sample),
                               get_group_label, character(1))

        p <- ggplot2::ggplot(
          avg_df,
          ggplot2::aes(x = reorder(Sample, avg_methylation),
                       y = avg_methylation, fill = Group,
                       text = paste0("Sample: ", Sample,
                                     "\nGroup: ", Group,
                                     "\nAvg Meth: ", round(avg_methylation, 1), "%"))) +
          ggplot2::geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
          ggplot2::scale_fill_manual(values = group_colors) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = "Average Methylation per Sample",
                        x = "", y = "Average Methylation (%)", fill = "Group") +
          ggplot2::theme_minimal(base_size = 12)

        plotly::ggplotly(p, tooltip = "text")
      })
    }

    # ══════════════════════════════════════════════════════════════
    # TAB 3: ACTIONS
    # ══════════════════════════════════════════════════════════════

    observeEvent(input$btn_update, {
      hi_p <- if (is.na(input$new_hi_perc)) 99.9 else input$new_hi_perc
      lo_p <- if (is.na(input$new_lo_perc)) NULL else input$new_lo_perc
      stopApp(list(action = "update", min_coverage = input$new_min_cov,
                   hi_perc = hi_p, lo_perc = lo_p,
                   samples_to_exclude = input$samples_to_remove,
                   chrs_to_keep = input$chrs_to_keep))
    })

    observeEvent(input$btn_proceed, {
      if (!is.null(input$samples_to_remove) && length(input$samples_to_remove) > 0) {
        showModal(modalDialog(
          title = "Warning",
          paste0("You selected ", length(input$samples_to_remove),
                 " sample(s) to remove but clicked 'Proceed'. ",
                 "They will NOT be removed unless you click 'Update' first. Proceed anyway?"),
          footer = tagList(modalButton("Cancel"),
                           actionButton("btn_proceed_confirm", "Yes, Proceed", class = "btn-danger"))
        ))
      } else {
        hi_p <- if (is.na(input$new_hi_perc)) 99.9 else input$new_hi_perc
        lo_p <- if (is.na(input$new_lo_perc)) NULL else input$new_lo_perc
        stopApp(list(action = "proceed", min_coverage = input$new_min_cov,
                     hi_perc = hi_p, lo_perc = lo_p,
                     samples_to_exclude = NULL, chrs_to_keep = input$chrs_to_keep))
      }
    })

    observeEvent(input$btn_proceed_confirm, {
      hi_p <- if (is.na(input$new_hi_perc)) 99.9 else input$new_hi_perc
      lo_p <- if (is.na(input$new_lo_perc)) NULL else input$new_lo_perc
      stopApp(list(action = "proceed", min_coverage = input$new_min_cov,
                   hi_perc = hi_p, lo_perc = lo_p,
                   samples_to_exclude = NULL, chrs_to_keep = input$chrs_to_keep))
    })
  }

  runApp(shinyApp(ui, server))
}
