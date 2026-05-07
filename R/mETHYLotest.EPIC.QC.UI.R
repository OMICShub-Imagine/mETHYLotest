#' Interactive EPIC QC Dashboard
#'
#' @description
#' Shiny Dashboard for EPIC QC visualization and sample removal.
#' Style matches \code{mETHYLotest.EPIC.ProjectUI}: \code{skin = "purple"},
#' \code{shinydashboard} layout. Provides seven visualization tabs:
#' PCA, t-SNE, SVD, Beta Distribution, Sample Correlation,
#' Detection P-value, and Sex Prediction.
#'
#' @param myLoad        ChAMP object with \code{$beta} and \code{$pd}.
#' @param pca_result    Output of \code{stats::prcomp}. \code{NULL} disables tab.
#' @param tsne_result   Data frame (TSNE1, TSNE2), row-named by sample.
#' @param svd_plot      ggplot from \code{mETHYLotest.EPIC.utils.SVDCorrelation}.
#' @param beta_density  ggplot from \code{mETHYLotest.EPIC.utils.PlotBetaDensity}.
#' @param cor_heatmap   ggplot from \code{mETHYLotest.EPIC.utils.PlotCorHeatmap}.
#' @param detp_summary  ggplot from \code{mETHYLotest.EPIC.utils.PlotDetPSummary}.
#' @param sex_pred      List from \code{mETHYLotest.EPIC.utils.PlotSexPrediction}
#'   with elements \code{$plot}, \code{$data}, \code{$mismatches}.
#'
#' @return Character vector of sample names flagged for removal (invisibly).
#'
#' @export
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import plotly
#' @import ggplot2
#' @importFrom DT DTOutput renderDT datatable dataTableProxy selectRows
mETHYLotest.EPIC.QC.UI <- function(myLoad,
                                  pca_result   = NULL,
                                  tsne_result  = NULL,
                                  svd_plot     = NULL,
                                  beta_density = NULL,
                                  cor_heatmap  = NULL,
                                  detp_summary = NULL,
                                  sex_pred     = NULL) {

  for (pkg in c("shiny", "shinydashboard", "shinyjs", "plotly", "DT",
                "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", pkg))
  }

  # ── Pre-processing ──────────────────────────────────────────────────────────
  sample_names  <- rownames(myLoad$pd)
  n_samples     <- length(sample_names)
  n_probes      <- nrow(myLoad$beta)
  pd_cols       <- colnames(myLoad$pd)

  default_color <- if ("Treatment"     %in% pd_cols) "Treatment"
  else if ("Sample_Status"  %in% pd_cols) "Sample_Status"
  else pd_cols[1L]

  pc_choices <- if (!is.null(pca_result))
    paste0("PC", seq_len(min(10L, ncol(pca_result$x))))
  else character(0)

  array_label <- if ("Sample_Plate" %in% pd_cols)
    paste(unique(myLoad$pd[["Sample_Plate"]]), collapse = ", ")
  else "Unknown"

  # Extract sex_pred components safely
  sex_plot       <- if (is.list(sex_pred) && !is.null(sex_pred$plot))
    sex_pred$plot   else NULL
  sex_data       <- if (is.list(sex_pred) && !is.null(sex_pred$data))
    sex_pred$data   else NULL
  sex_mismatches <- if (is.list(sex_pred) && !is.null(sex_pred$mismatches))
    sex_pred$mismatches else character(0)

  # Outlier samples from correlation heatmap (stored as attribute)
  cor_outliers   <- if (!is.null(cor_heatmap))
    attr(cor_heatmap, "outlier_samples", exact = TRUE)
  else character(0)
  if (is.null(cor_outliers)) cor_outliers <- character(0)

  # Helper: make a tab unavailable warning box
  unavailable_box <- function(label) {
    div(class = "alert alert-warning", style = "margin:15px;",
        icon("exclamation-triangle"),
        paste0(" ", label, " not available."))
  }

  # --------------------------------------------------------------------------
  # UI
  # --------------------------------------------------------------------------
  ui <- shinydashboard::dashboardPage(
    skin = "purple",

    shinydashboard::dashboardHeader(title = "mETHYLotest EPIC QC"),

    # ── Sidebar ────────────────────────────────────────────────────────────────
    shinydashboard::dashboardSidebar(
      width = 300,

      div(style = "padding:15px;",

          h5("Dataset Summary", style = "color:#b8c7ce;"),
          div(class = "callout callout-info",
              style = paste0("padding:10px; border-left-color:#00c0ef;",
                             " background-color:#222d32 !important;",
                             " border:1px solid #444;"),
              p(strong("Samples: "),  n_samples),
              p(strong("Probes: "),   format(n_probes, big.mark = ",")),
              p(strong("Array(s): "), array_label)
          ),

          br(),
          uiOutput("removal_summary_sidebar")
      ),

      shinydashboard::sidebarMenu(
        id = "tabs",
        shinydashboard::menuItem("1. QC Plots",       tabName = "plots",
                                 icon = icon("chart-bar")),
        shinydashboard::menuItem("2. Data Tables",    tabName = "tables",
                                 icon = icon("table")),
        shinydashboard::menuItem("3. Remove Samples", tabName = "samples",
                                 icon = icon("filter"),
                                 badgeLabel = "Action",
                                 badgeColor = "orange")
      )
    ),

    # ── Body ───────────────────────────────────────────────────────────────────
    shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),

      shinydashboard::tabItems(

        # ── Tab 1 : QC Plots ─────────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "plots",

                                fluidRow(
                                  shinydashboard::infoBoxOutput("ibox_samples",  width = 4),
                                  shinydashboard::infoBoxOutput("ibox_probes",   width = 4),
                                  shinydashboard::infoBoxOutput("ibox_flagged",  width = 4)
                                ),

                                fluidRow(
                                  shinydashboard::tabBox(
                                    title = tagList(icon("chart-bar"), " QC Visualizations"),
                                    id    = "tabset_plots",
                                    width = 12,

                                    # ---- PCA ---------------------------------------------------
                                    shiny::tabPanel("PCA",
                                                    if (is.null(pca_result)) {
                                                      unavailable_box("PCA result")
                                                    } else {
                                                      tagList(
                                                        br(),
                                                        fluidRow(
                                                          column(3, selectInput("pca_x", "X Axis:",
                                                                                choices  = pc_choices,
                                                                                selected = pc_choices[1L])),
                                                          column(3, selectInput("pca_y", "Y Axis:",
                                                                                choices  = pc_choices,
                                                                                selected = pc_choices[2L])),
                                                          column(3, selectInput("pca_color", "Colour by:",
                                                                                choices  = pd_cols,
                                                                                selected = default_color)),
                                                          column(3, br(),
                                                                 checkboxInput("pca_labels",
                                                                               "Show sample labels",
                                                                               value = TRUE))
                                                        ),
                                                        plotly::plotlyOutput("plot_pca", height = "500px")
                                                      )
                                                    }
                                    ),

                                    # ---- t-SNE ------------------------------------------------
                                    shiny::tabPanel("t-SNE",
                                                    if (is.null(tsne_result)) {
                                                      unavailable_box("t-SNE result")
                                                    } else {
                                                      tagList(
                                                        br(),
                                                        fluidRow(
                                                          column(4, selectInput("tsne_color", "Colour by:",
                                                                                choices  = pd_cols,
                                                                                selected = default_color)),
                                                          column(4, br(),
                                                                 checkboxInput("tsne_labels",
                                                                               "Show sample labels",
                                                                               value = TRUE))
                                                        ),
                                                        plotly::plotlyOutput("plot_tsne", height = "500px")
                                                      )
                                                    }
                                    ),

                                    # ---- SVD --------------------------------------------------
                                    shiny::tabPanel("SVD",
                                                    if (is.null(svd_plot)) {
                                                      unavailable_box("SVD result")
                                                    } else {
                                                      tagList(br(),
                                                              shiny::plotOutput("plot_svd", height = "500px"))
                                                    }
                                    ),

                                    # ---- Beta Distribution --------------------------------------------------
                                    shiny::tabPanel("Beta Distribution",
                                                    if (is.null(beta_density)) {
                                                      unavailable_box("Beta density")
                                                    } else {
                                                      tagList(
                                                        br(),
                                                        plotly::plotlyOutput("plot_beta", height = "500px")  # changed
                                                      )
                                                    }
                                    ),

                                    # ---- Sample Correlation -------------------------------------------------
                                    shiny::tabPanel("Sample Correlation",
                                                    if (is.null(cor_heatmap)) {
                                                      unavailable_box("Correlation heatmap")
                                                    } else {
                                                      tagList(
                                                        br(),
                                                        if (length(cor_outliers) > 0L)
                                                          div(class = "alert alert-warning",
                                                              style = "margin:0 0 10px 0; font-size:12px;",
                                                              icon("exclamation-triangle"),
                                                              strong(" Potential outlier(s): "),
                                                              paste(cor_outliers, collapse = ", "))
                                                        else
                                                          div(class = "alert alert-success",
                                                              style = "margin:0 0 10px 0; font-size:12px;",
                                                              icon("check"),
                                                              " No outliers detected."),
                                                        plotly::plotlyOutput("plot_cor", height = "560px")  # changed
                                                      )
                                                    }
                                    ),

                                    # ---- Detection P-value -----------------------------------
                                    shiny::tabPanel("Detection P-value",
                                                    if (is.null(detp_summary)) {
                                                      unavailable_box("Detection p-value summary")
                                                    } else {
                                                      tagList(br(),
                                                              plotly::plotlyOutput("plot_detp",
                                                                                   height = "500px"))
                                                    }
                                    ),

                                    # ---- Sex Prediction --------------------------------------
                                    shiny::tabPanel("Sex Prediction",
                                                    if (is.null(sex_plot)) {
                                                      unavailable_box("Sex prediction")
                                                    } else {
                                                      tagList(
                                                        br(),
                                                        if (length(sex_mismatches) > 0)
                                                          div(class = "alert alert-danger",
                                                              style = "margin:0 0 10px 0; font-size:12px;",
                                                              icon("exclamation-triangle"),
                                                              strong(" Sex mismatch(es) detected: "),
                                                              paste(sex_mismatches, collapse = ", "),
                                                              " -- verify in the table below.")
                                                        else if (!is.null(sex_data) &&
                                                                 "Match" %in% colnames(sex_data))
                                                          div(class = "alert alert-success",
                                                              style = "margin:0 0 10px 0; font-size:12px;",
                                                              icon("check"),
                                                              " Declared and predicted sex match for all samples."),
                                                        plotly::plotlyOutput("plot_sex", height = "460px"),
                                                        if (!is.null(sex_data)) {
                                                          tagList(
                                                            hr(),
                                                            h5("Sex Prediction Details:"),
                                                            DT::DTOutput("table_sex_pred")
                                                          )
                                                        }
                                                      )
                                                    }
                                    )
                                  )
                                )
        ),

        # ── Tab 2 : Data Tables ───────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tables",

                                shinydashboard::box(
                                  title = "Phenotype Data", status = "primary",
                                  solidHeader = TRUE, width = 12,
                                  collapsible = TRUE, collapsed = FALSE,
                                  DT::DTOutput("table_pd")
                                ),

                                if (!is.null(pca_result))
                                  shinydashboard::box(
                                    title = "PCA Scores", status = "info",
                                    solidHeader = TRUE, width = 12,
                                    collapsible = TRUE, collapsed = TRUE,
                                    DT::DTOutput("table_pca")
                                  ),

                                if (!is.null(tsne_result))
                                  shinydashboard::box(
                                    title = "t-SNE Coordinates", status = "info",
                                    solidHeader = TRUE, width = 12,
                                    collapsible = TRUE, collapsed = TRUE,
                                    DT::DTOutput("table_tsne")
                                  )
        ),

        # ── Tab 3 : Remove Samples ────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "samples",

                                # Warnings row
                                if (length(sex_mismatches) > 0 || length(cor_outliers) > 0) {
                                  fluidRow(
                                    shinydashboard::box(
                                      title = "Automated Warnings", status = "danger",
                                      solidHeader = TRUE, width = 12,
                                      collapsible = TRUE, collapsed = FALSE,

                                      if (length(sex_mismatches) > 0)
                                        div(class = "alert alert-danger",
                                            style = "margin-bottom:8px;",
                                            icon("venus-mars"),
                                            strong(" Sex mismatch: "),
                                            paste(sex_mismatches, collapse = ", "),
                                            tags$small(" -- See 'Sex Prediction' tab for details.")),

                                      if (length(cor_outliers) > 0)
                                        div(class = "alert alert-warning",
                                            icon("project-diagram"),
                                            strong(" Correlation outlier(s): "),
                                            paste(cor_outliers, collapse = ", "),
                                            tags$small(" -- See 'Sample Correlation' tab for details."))
                                    )
                                  )
                                },

                                fluidRow(
                                  # ---- Sample selection -----------------------------------------
                                  shinydashboard::box(
                                    title = "Sample Exclusion", status = "warning",
                                    solidHeader = TRUE, width = 8,

                                    div(class = "alert alert-info",
                                        icon("info-circle"),
                                        " Click rows to flag samples for ",
                                        strong("REMOVAL"), ". ",
                                        "Unselected rows are kept."),

                                    DT::DTOutput("sample_table"),
                                    br(),

                                    fluidRow(
                                      column(4,
                                             actionButton("btn_deselect_all", "Deselect All",
                                                          icon  = icon("times-circle"),
                                                          class = "btn-warning btn-block")
                                      ),
                                      column(4,
                                             actionButton("btn_select_all", "Select All",
                                                          icon  = icon("check-square"),
                                                          class = "btn-danger btn-block")
                                      )
                                    )
                                  ),

                                  # ---- Validation -----------------------------------------------
                                  shinydashboard::box(
                                    title = "Validation", status = "success",
                                    solidHeader = TRUE, width = 4,

                                    h4("Removal Summary"),
                                    uiOutput("removal_summary_main"),

                                    hr(),

                                    div(style = paste0("text-align:center;",
                                                       " background-color:#f9f9f9;",
                                                       " padding:15px;",
                                                       " border-radius:5px;",
                                                       " border:1px solid #ddd;"),
                                        strong("Confirm and proceed"),
                                        p("Accept selection and return to pipeline."),
                                        actionButton("btn_proceed", " Confirm & Proceed",
                                                     icon  = icon("check"),
                                                     class = "btn-success btn-lg",
                                                     style = "width:80%;")
                                    )
                                  )
                                )
        )
      )
    )
  )

  # --------------------------------------------------------------------------
  # Server
  # --------------------------------------------------------------------------
  server <- function(input, output, session) {

    flagged <- reactiveVal(character(0))

    # ── Info boxes ────────────────────────────────────────────────────────────
    output$ibox_samples <- shinydashboard::renderInfoBox(
      shinydashboard::infoBox("Total Samples", n_samples,
                              icon = icon("vial"),
                              color = "purple", fill = TRUE)
    )
    output$ibox_probes <- shinydashboard::renderInfoBox(
      shinydashboard::infoBox("CpG Probes",
                              format(n_probes, big.mark = ","),
                              icon = icon("dna"),
                              color = "blue", fill = TRUE)
    )
    output$ibox_flagged <- shinydashboard::renderInfoBox({
      rem <- length(flagged())
      shinydashboard::infoBox("Flagged for Removal", rem,
                              icon  = icon("trash"),
                              color = if (rem > 0L) "red" else "green",
                              fill  = TRUE)
    })

    # ── Sidebar summary ───────────────────────────────────────────────────────
    output$removal_summary_sidebar <- renderUI({
      rem <- flagged()
      if (length(rem) > 0L) {
        div(class = "callout callout-warning",
            style = paste0("padding:10px; border-left-color:#f39c12;",
                           " background-color:#222d32 !important;",
                           " border:1px solid #444;"),
            p(style = "color:#f39c12;",
              strong("Flagged: "), length(rem), " sample(s)"),
            p(style = "color:#f39c12;",
              strong("Kept: "),    n_samples - length(rem)),
            tags$small(style = "color:#f39c12; word-break:break-all;",
                       paste(rem, collapse = ", "))
        )
      } else {
        div(class = "callout callout-success",
            style = paste0("padding:10px; border-left-color:#00a65a;",
                           " background-color:#222d32 !important;",
                           " border:1px solid #444;"),
            p(style = "color:#00a65a;",
              icon("check"), " No samples flagged.")
        )
      }
    })

    # ── Main removal summary ──────────────────────────────────────────────────
    output$removal_summary_main <- renderUI({
      rem  <- flagged()
      kept <- setdiff(sample_names, rem)
      tagList(
        div(class = "alert alert-success", style = "padding:8px;",
            icon("check"),
            strong(paste0(" Kept: ", length(kept)))),
        div(class = "alert alert-danger", style = "padding:8px;",
            icon("trash"),
            strong(paste0(" Removed: ", length(rem))),
            if (length(rem) > 0L)
              tags$small(br(), paste(rem, collapse = ", ")))
      )
    })

    # ── Scatter helper ────────────────────────────────────────────────────────
    make_scatter <- function(df, x_col, y_col, x_lab, y_lab,
                             color_col, show_labels) {
      df$Status <- ifelse(df$Sample %in% flagged(), "Flagged", "Keep")
      mode_plot <- if (isTRUE(show_labels)) "markers+text" else "markers"
      plotly::plot_ly(
        data         = df,
        x            = df[[x_col]],
        y            = df[[y_col]],
        color        = if (color_col %in% colnames(df)) df[[color_col]]
        else rep("Sample", nrow(df)),
        symbol       = ~Status,
        symbols      = c("Keep" = "circle", "Flagged" = "x"),
        text         = ~Sample,
        type         = "scatter",
        mode         = mode_plot,
        textposition = "top center",
        marker       = list(size = 10)
      ) %>%
        plotly::layout(xaxis = list(title = x_lab),
                       yaxis = list(title = y_lab))
    }

    # ── PCA ───────────────────────────────────────────────────────────────────
    output$plot_pca <- plotly::renderPlotly({
      req(pca_result, input$pca_x, input$pca_y, input$pca_color)
      df           <- as.data.frame(pca_result$x)
      df$Sample    <- rownames(df)
      if (input$pca_color %in% pd_cols)
        df[[input$pca_color]] <- myLoad$pd[df$Sample, input$pca_color]
      var_exp <- summary(pca_result)$importance[2L, ] * 100
      make_scatter(df, input$pca_x, input$pca_y,
                   sprintf("%s (%.1f%%)", input$pca_x,
                           var_exp[[input$pca_x]]),
                   sprintf("%s (%.1f%%)", input$pca_y,
                           var_exp[[input$pca_y]]),
                   input$pca_color, input$pca_labels)
    })

    # ── t-SNE ─────────────────────────────────────────────────────────────────
    output$plot_tsne <- plotly::renderPlotly({
      req(tsne_result, input$tsne_color)
      df        <- tsne_result
      df$Sample <- rownames(df)
      if (input$tsne_color %in% pd_cols)
        df[[input$tsne_color]] <- myLoad$pd[df$Sample, input$tsne_color]
      make_scatter(df, "TSNE1", "TSNE2", "t-SNE 1", "t-SNE 2",
                   input$tsne_color, input$tsne_labels)
    })

    # ── SVD ───────────────────────────────────────────────────────────────────
    output$plot_svd <- shiny::renderPlot({
      req(svd_plot)
      print(svd_plot)
    })

    # ── Beta density — renderPlotly (plotly object natif) ─────────────────────
  output$plot_beta <- plotly::renderPlotly({
    req(beta_density)
    beta_density   # already a plotly object — no ggplotly() needed
  })

  # ── Correlation heatmap — renderPlotly (plotly object natif) ──────────────
  output$plot_cor <- plotly::renderPlotly({
    req(cor_heatmap)
    cor_heatmap    # already a plotly object — no ggplotly() needed
  })

    # ── Detection p-value (plotly) ────────────────────────────────────────────
    output$plot_detp <- plotly::renderPlotly({
      req(detp_summary)
      plotly::ggplotly(detp_summary)
    })

    # ── Sex prediction (plotly) ───────────────────────────────────────────────
    output$plot_sex <- plotly::renderPlotly({
      req(sex_plot)
      plotly::ggplotly(sex_plot) %>%
        plotly::layout(legend = list(orientation = "v"))
    })

    # ── Sex prediction table ──────────────────────────────────────────────────
    output$table_sex_pred <- DT::renderDT({
      req(sex_data)
      DT::datatable(
        sex_data,
        options  = list(dom = "tp", pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      ) %>%
        DT::formatStyle(
          if ("Match" %in% colnames(sex_data)) "Match" else character(0),
          backgroundColor = DT::styleEqual(
            c("Match", "Mismatch"),
            c("#d4edda",  "#f8d7da")
          )
        )
    })

    # ── Data tables ───────────────────────────────────────────────────────────
    output$table_pd <- DT::renderDT(
      DT::datatable(myLoad$pd,
                    options = list(dom = "ftp", pageLength = 10,
                                   scrollX = TRUE))
    )
    output$table_pca <- DT::renderDT({
      req(pca_result)
      DT::datatable(round(as.data.frame(pca_result$x), 4L),
                    options = list(dom = "ftp", pageLength = 10,
                                   scrollX = TRUE))
    })
    output$table_tsne <- DT::renderDT({
      req(tsne_result)
      DT::datatable(round(tsne_result, 4L),
                    options = list(dom = "ftp", pageLength = 10,
                                   scrollX = TRUE))
    })

    # ── Sample table (Tab 3) ──────────────────────────────────────────────────
    output$sample_table <- DT::renderDT(
      DT::datatable(
        myLoad$pd,
        selection = list(mode = "multiple", selected = NULL, target = "row"),
        options   = list(dom = "ftp", pageLength = 15, scrollX = TRUE)
      )
    )

    proxy_table <- DT::dataTableProxy("sample_table")

    observeEvent(input$sample_table_rows_selected, {
      sel <- input$sample_table_rows_selected
      flagged(if (length(sel) > 0L) rownames(myLoad$pd)[sel]
              else character(0))
    }, ignoreNULL = FALSE)

    observeEvent(input$btn_select_all, {
      DT::selectRows(proxy_table, seq_len(nrow(myLoad$pd)))
    })
    observeEvent(input$btn_deselect_all, {
      DT::selectRows(proxy_table, NULL)
      flagged(character(0))
    })

    # ── Proceed ───────────────────────────────────────────────────────────────
    observeEvent(input$btn_proceed, {
      rem <- flagged()
      message("[mETHYLotest EPIC QC] Confirmed: ",
              length(rem), " sample(s) removed, ",
              n_samples - length(rem), " kept.")
      shiny::stopApp(returnValue = rem)
    })

    session$onSessionEnded(function() {
      shiny::stopApp(returnValue = isolate(flagged()))
    })
  }

  message("[mETHYLotest] Starting EPIC QC UI...")
  invisible(shiny::runApp(shiny::shinyApp(ui, server),
                          host = "0.0.0.0",
                          launch.browser = TRUE))
}
