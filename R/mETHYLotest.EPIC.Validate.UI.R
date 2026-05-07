#' mETHYLotest Signature Validator UI
#'
#' @description
#' Shiny Dashboard for interactive signature validation exploration.
#' Displays SVM performance, PCA, Silhouette, Heatmap and variable importance.
#' Style matches the rest of mETHYLotest (skin = "purple").
#'
#' @param validation_results Named list returned by
#'   \code{mETHYLotest.EPIC.validate()}.
#'
#' @return NULL (invisibly). The app blocks until closed.
#' @export
#' @import shiny shinydashboard
#' @import plotly ggplot2
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom heatmaply heatmaply
#' @importFrom caret varImp confusionMatrix
#' @importFrom pROC roc
#' @importFrom cluster silhouette
mETHYLotest.EPIC.Validate.UI <- function(validation_results) {

  for (pkg in c("shiny", "shinydashboard", "plotly", "ggplot2",
                "DT", "caret", "heatmaply", "pROC"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", pkg))

  if (length(validation_results) == 0L)
    stop("[Validate UI] No validation results provided.")

  # ========================================================================
  # UI
  # ========================================================================
  ui <- shinydashboard::dashboardPage(
    skin = "purple",

    shinydashboard::dashboardHeader(
      title = "mETHYLotest Signature Validator"
    ),

    shinydashboard::dashboardSidebar(
      width = 280,
      div(style = "padding:15px;",

          h4("Signature Selection",
             style = "color:#b8c7ce;"),

          selectInput("selected_sig", "Signature:",
                      choices = names(validation_results)),

          hr(style = "border-color:#444;"),

          h5("Details", style = "color:#b8c7ce;"),
          uiOutput("sig_details_ui"),

          hr(style = "border-color:#444;"),

          div(class = "callout callout-info",
              style = paste0("padding:8px; border-left-color:#00c0ef;",
                             " background-color:#222d32 !important;",
                             " border:1px solid #444;"),
              tags$small(
                style = "color:#b8c7ce;",
                icon("info-circle"),
                " Heatmap & PCA are available even if SVM fails",
                " (e.g. low sample size)."))
      ),

      shinydashboard::sidebarMenu(
        id = "val_tabs",
        shinydashboard::menuItem("Metrics Overview",
                                 tabName = "tab_metrics",
                                 icon = icon("tachometer-alt")),
        shinydashboard::menuItem("PCA",
                                 tabName = "tab_pca",
                                 icon = icon("braille")),
        shinydashboard::menuItem("Heatmap",
                                 tabName = "tab_heatmap",
                                 icon = icon("border-all")),
        shinydashboard::menuItem("SVM Performance",
                                 tabName = "tab_svm",
                                 icon = icon("chart-line")),
        shinydashboard::menuItem("Top Contributors",
                                 tabName = "tab_varimp",
                                 icon = icon("sort-amount-down"))
      ),

      div(style = "padding:15px; margin-top:10px;",
          actionButton("close_val_btn", " Close & Return",
                       icon  = icon("power-off"),
                       class = "btn-danger btn-block"))
    ),

    shinydashboard::dashboardBody(

      shinydashboard::tabItems(

        # ── Metrics Overview ─────────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_metrics",
          fluidRow(
            shinydashboard::valueBoxOutput("vbox_metric1", width = 4),
            shinydashboard::valueBoxOutput("vbox_metric2", width = 4),
            shinydashboard::valueBoxOutput("vbox_sil",     width = 4)
          ),
          fluidRow(
            shinydashboard::box(
              title = "Signature Summary", status = "primary",
              solidHeader = TRUE, width = 12,
              DT::DTOutput("benchmark_table")
            )
          )
        ),

        # ── PCA ──────────────────────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_pca",
          fluidRow(
            shinydashboard::box(
              title = "PCA (Unsupervised Separation)",
              status = "primary", solidHeader = TRUE, width = 12,
              p("Does the signature separate biological groups",
                " without supervision?"),
              plotly::plotlyOutput("pcaPlot", height = "550px")
            )
          )
        ),

        # ── Heatmap ─────────────────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_heatmap",
          fluidRow(
            shinydashboard::box(
              title = "Hierarchical Clustering",
              status = "primary", solidHeader = TRUE, width = 12,
              p("Clustering of samples based on the signature CpGs."),
              plotly::plotlyOutput("heatmapPlot", height = "700px")
            )
          )
        ),

        # ── SVM Performance ──────────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_svm",
          fluidRow(
            shinydashboard::box(
              title = "SVM Classification Performance",
              status = "success", solidHeader = TRUE, width = 12,
              uiOutput("svm_performance_ui")
            )
          )
        ),

        # ── Variable Importance ──────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_varimp",
          fluidRow(
            shinydashboard::box(
              title = "Top Discriminative CpGs (SVM Weights)",
              status = "warning", solidHeader = TRUE, width = 12,
              uiOutput("svm_importance_ui")
            )
          )
        )
      )
    )
  )

  # ========================================================================
  # Server
  # ========================================================================
  server <- function(input, output, session) {

    curr <- reactive({
      req(input$selected_sig)
      validation_results[[input$selected_sig]]
    })

    # ── Sidebar details ──────────────────────────────────────────────────
    output$sig_details_ui <- renderUI({
      res <- curr()
      metric_lbl <- if (is.null(res$svm_model)) "None (SVM failed)"
      else res$metric_used

      perf_val <- if (!is.null(res$svm_model))
        round(max(res$svm_model$results[[res$metric_used]],
                  na.rm = TRUE), 3)
      else "N/A"

      tagList(
        p(style = "color:#b8c7ce;",
          strong("CpGs: "), res$n_cpgs),
        p(style = "color:#b8c7ce;",
          strong("Metric: "), metric_lbl),
        p(style = "color:#b8c7ce;",
          strong("Performance: "), perf_val),
        p(style = "color:#b8c7ce;",
          strong("Silhouette: "), round(res$silhouette, 3))
      )
    })

    # ── Value boxes ──────────────────────────────────────────────────────
    output$vbox_metric1 <- shinydashboard::renderValueBox({
      res <- curr()
      if (is.null(res$svm_model)) {
        shinydashboard::valueBox(
          "FAILED", "SVM Result",
          icon = icon("times-circle"), color = "red")
      } else {
        lbl <- if (res$metric_used == "ROC") "ROC (AUC)" else "Accuracy"
        val <- round(max(res$svm_model$results[[res$metric_used]],
                         na.rm = TRUE), 3)
        shinydashboard::valueBox(
          val, lbl,
          icon = icon("bullseye"), color = "blue")
      }
    })

    output$vbox_metric2 <- shinydashboard::renderValueBox({
      res <- curr()
      if (is.null(res$svm_model)) {
        shinydashboard::valueBox(
          "N/A", "Kappa / Sensitivity",
          icon = icon("question"), color = "yellow")
      } else {
        if (res$metric_used == "ROC") {
          val <- round(max(res$svm_model$results$Sens,
                           na.rm = TRUE), 3)
          shinydashboard::valueBox(
            val, "Sensitivity",
            icon = icon("heartbeat"), color = "green")
        } else {
          val <- round(max(res$svm_model$results$Kappa,
                           na.rm = TRUE), 3)
          shinydashboard::valueBox(
            val, "Kappa",
            icon = icon("chart-bar"), color = "green")
        }
      }
    })

    output$vbox_sil <- shinydashboard::renderValueBox({
      val <- round(curr()$silhouette, 3)
      col <- if (val > 0.5) "green"
      else if (val > 0.25) "yellow"
      else "red"
      shinydashboard::valueBox(
        val, "Silhouette",
        icon = icon("layer-group"), color = col)
    })

    # ── Benchmark table (all signatures) ─────────────────────────────────
    output$benchmark_table <- DT::renderDT({
      rows <- lapply(names(validation_results), function(nm) {
        r <- validation_results[[nm]]
        perf <- if (!is.null(r$svm_model))
          round(max(r$svm_model$results[[r$metric_used]],
                    na.rm = TRUE), 3)
        else NA
        data.frame(
          Signature  = nm,
          CpGs       = r$n_cpgs,
          Metric     = if (!is.null(r$svm_model)) r$metric_used
          else "N/A",
          Performance = perf,
          Silhouette = round(r$silhouette, 3),
          stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, rows)
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 20,
                                   order = list(list(3, "desc"))))
    })

    # ── PCA ──────────────────────────────────────────────────────────────
    output$pcaPlot <- plotly::renderPlotly({
      res <- curr()
      if (is.null(res$pca))
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "PCA unavailable."))

      pca <- res$pca
      df  <- data.frame(
        PC1    = pca$x[, 1L],
        PC2    = pca$x[, 2L],
        Group  = res$raw_data$Class,
        Sample = rownames(res$raw_data)
      )
      ve <- round(summary(pca)$importance[2, 1:2] * 100, 1)

      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = PC1, y = PC2, color = Group, text = Sample)) +
        ggplot2::geom_point(size = 3, alpha = 0.8) +
        ggplot2::stat_ellipse(level = 0.95, linetype = "dashed",
                              alpha = 0.4) +
        ggplot2::labs(
          x = paste0("PC1 (", ve[1], "%)"),
          y = paste0("PC2 (", ve[2], "%)"),
          title = paste0("PCA — ", input$selected_sig)) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_brewer(palette = "Set1")

      plotly::ggplotly(p, tooltip = "text")
    })

    # ── Heatmap ──────────────────────────────────────────────────────────
    output$heatmapPlot <- plotly::renderPlotly({
      res <- curr()
      annot_col <- data.frame(Group = res$raw_data$Class)
      rownames(annot_col) <- rownames(res$raw_data)

      heatmaply::heatmaply(
        res$beta_matrix,
        limits         = c(0, 1),
        colors         = grDevices::colorRampPalette(
          c("navy", "white", "firebrick"))(100),
        col_side_colors = annot_col,
        row_dend_left  = TRUE,
        col_dend_top   = TRUE,
        branches_lwd   = 0.2,
        labRow = if (nrow(res$beta_matrix) > 50) NULL
        else rownames(res$beta_matrix),
        main   = paste0("Clustering — ", res$n_cpgs, " CpGs")
      )
    })

    # ── SVM Performance ──────────────────────────────────────────────────
    output$svm_performance_ui <- renderUI({
      res <- curr()
      if (is.null(res$svm_model)) {
        return(div(
          class = "alert alert-danger",
          style = "padding:15px;",
          icon("exclamation-triangle"),
          strong(" SVM analysis failed."),
          br(), br(),
          "Sample size too small or classes too unbalanced.",
          br(),
          "Please rely on the PCA and Heatmap tabs."))
      }

      fluidRow(
        column(6,
               h4("Confusion Matrix"),
               plotOutput("confMatPlot", height = "400px")),
        column(6, uiOutput("roc_panel"))
      )
    })

    output$confMatPlot <- renderPlot({
      res <- curr()
      if (is.null(res$svm_model)) return(NULL)

      preds <- res$svm_model$pred
      preds <- preds[preds$C == res$svm_model$bestTune$C, ]
      cm    <- caret::confusionMatrix(preds$pred, preds$obs)
      df_cm <- as.data.frame(cm$table)

      ggplot2::ggplot(df_cm,
                      ggplot2::aes(x = Reference, y = Prediction,
                                   fill = Freq)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = Freq),
                           size = 8, color = "white") +
        ggplot2::scale_fill_gradient(low = "#868e96",
                                     high = "#0275d8") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Confusion Matrix")
    })

    output$roc_panel <- renderUI({
      res <- curr()
      if (is.null(res$svm_model)) return(NULL)
      if (res$metric_used == "ROC") {
        tagList(
          h4("ROC Curve"),
          plotOutput("rocPlot", height = "400px"))
      } else {
        div(class = "alert alert-warning",
            style = "padding:10px; margin-top:30px;",
            icon("info-circle"),
            strong(" ROC unavailable."), br(),
            "Accuracy mode was used (sample size too small",
            " for class probabilities).")
      }
    })

    output$rocPlot <- renderPlot({
      res <- curr()
      if (res$metric_used != "ROC" || is.null(res$svm_model))
        return(NULL)

      preds   <- res$svm_model$pred
      preds   <- preds[preds$C == res$svm_model$bestTune$C, ]
      roc_obj <- pROC::roc(preds$obs,
                           preds[[levels(preds$obs)[1L]]])
      plot(roc_obj,
           print.auc         = TRUE,
           auc.polygon       = TRUE,
           grid               = c(0.1, 0.2),
           grid.col           = c("green", "red"),
           max.auc.polygon    = TRUE,
           auc.polygon.col    = "lightblue",
           print.thres        = TRUE,
           legacy.axes        = TRUE)
    })

    # ── Variable Importance ──────────────────────────────────────────────
    output$svm_importance_ui <- renderUI({
      res <- curr()
      if (is.null(res$svm_model))
        return(div(
          class = "alert alert-danger",
          style = "padding:15px;",
          icon("exclamation-triangle"),
          strong(" Requires a valid SVM model."),
          br(), "SVM failed for this signature."))

      plotly::plotlyOutput("varImpPlot", height = "600px")
    })

    output$varImpPlot <- plotly::renderPlotly({
      res <- curr()
      if (is.null(res$svm_model)) return(NULL)

      imp_obj <- tryCatch(
        caret::varImp(res$svm_model, scale = FALSE),
        error = function(e) NULL)
      if (is.null(imp_obj)) return(NULL)

      imp <- imp_obj$importance
      if (nrow(imp) == 0L) return(NULL)

      # Score: Overall or max across classes
      imp$Score <- if ("Overall" %in% colnames(imp)) imp$Overall
      else apply(imp, 1, max)

      # CpG name recovery (caret renames '-' to '.')
      imp$RawName <- rownames(imp)
      known <- names(res$genes_map)

      recover_name <- function(x) {
        if (x %in% known) return(x)
        x2 <- gsub("\\.", "-", x)
        if (x2 %in% known) return(x2)
        x3 <- gsub("^X", "", x)
        if (x3 %in% known) return(x3)
        x
      }
      imp$CleanCpG <- vapply(imp$RawName, recover_name, character(1))

      # Gene mapping
      imp$Gene <- res$genes_map[imp$CleanCpG]
      imp$Gene[is.na(imp$Gene) | imp$Gene == ""] <- "Unknown"
      imp$Label <- paste0(imp$CleanCpG, " (", imp$Gene, ")")

      # Top 20
      imp <- imp[order(imp$Score, decreasing = TRUE), ]
      plot_data <- head(imp, 20L)
      plot_data$Label <- factor(plot_data$Label,
                                levels = rev(plot_data$Label))

      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = Label, y = Score)) +
        ggplot2::geom_bar(stat = "identity", fill = "#d9534f",
                          alpha = 0.8) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Top 20 Discriminative CpGs",
                      x = "", y = "Variable Importance") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))

      plotly::ggplotly(p)
    })

    # ── Close ────────────────────────────────────────────────────────────
    observeEvent(input$close_val_btn, shiny::stopApp())
    session$onSessionEnded(function() shiny::stopApp())
  }

  message("[mETHYLotest] Starting Signature Validator UI...")
  invisible(shiny::runApp(shiny::shinyApp(ui, server),
                          host = "0.0.0.0",
                          launch.browser = TRUE))
}
