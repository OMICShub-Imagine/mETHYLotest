#' mETHYLotest NGS Signature Validator UI
#'
#' @description
#' Shiny Dashboard for interactive signature validation exploration.
#' Adapted for NGS Long-Read data from methylKit objects.
#' Displays SVM performance, PCA, Silhouette, Heatmap and variable importance.
#'
#' @param validation_results Named list returned by
#'   \code{mETHYLotest.NGS.validate()}.
#'
#' @return NULL (invisibly). The app blocks until closed.
#' @export
mETHYLotest.NGS.Validate.UI <- function(validation_results) {

  for (pkg in c("shiny", "shinydashboard", "plotly", "ggplot2",
                "DT", "caret", "heatmaply", "pROC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("[Validate UI] Package '%s' is required.", pkg))
    }
  }

  if (length(validation_results) == 0L) {
    stop("[Validate UI] No validation results provided.")
  }

  # ========================================================================
  # UI
  # ========================================================================
  ui <- shinydashboard::dashboardPage(
    skin = "blue",

    shinydashboard::dashboardHeader(
      title = "mETHYLotest NGS Validator"
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

          div(style = paste0("padding:8px; border-left:3px solid #00c0ef;",
                             " background-color:#222d32;",
                             " border:1px solid #444; border-radius:4px;"),
              tags$small(
                style = "color:#b8c7ce;",
                icon("info-circle"),
                " NGS signatures use genomic loci (chr:start-end).",
                " Heatmap & PCA are available even if SVM fails."))
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
                                 icon = icon("sort-amount-down")),
        shinydashboard::menuItem("CpG Details",
                                 tabName = "tab_cpgs",
                                 icon = icon("dna"))
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
            shinydashboard::valueBoxOutput("vbox_ncpgs",   width = 3),
            shinydashboard::valueBoxOutput("vbox_metric1", width = 3),
            shinydashboard::valueBoxOutput("vbox_metric2", width = 3),
            shinydashboard::valueBoxOutput("vbox_sil",     width = 3)
          ),
          fluidRow(
            shinydashboard::box(
              title = "All Signatures Comparison", status = "primary",
              solidHeader = TRUE, width = 12,
              DT::DTOutput("benchmark_table")
            )
          ),
          fluidRow(
            shinydashboard::box(
              title = "Silhouette Interpretation", status = "info",
              solidHeader = TRUE, width = 6, collapsible = TRUE,
              uiOutput("sil_interpretation_ui")
            ),
            shinydashboard::box(
              title = "Confusion Matrix", status = "success",
              solidHeader = TRUE, width = 6, collapsible = TRUE,
              uiOutput("conf_matrix_ui")
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
              p("Does the signature separate Control vs Test without supervision?"),
              plotly::plotlyOutput("pcaPlot", height = "550px")
            )
          ),
          fluidRow(
            shinydashboard::box(
              title = "Scree Plot", status = "info",
              solidHeader = TRUE, width = 6,
              plotOutput("screePlot", height = "300px")
            ),
            shinydashboard::box(
              title = "PCA Variance", status = "info",
              solidHeader = TRUE, width = 6,
              DT::DTOutput("pca_var_table")
            )
          )
        ),

        # ── Heatmap ─────────────────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_heatmap",
          fluidRow(
            shinydashboard::box(
              title = "Hierarchical Clustering (Methylation Heatmap)",
              status = "primary", solidHeader = TRUE, width = 12,
              p("Clustering of samples based on signature CpG methylation levels."),
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
        ),

        # ── CpG Details ──────────────────────────────────────────────────
        shinydashboard::tabItem(
          tabName = "tab_cpgs",
          fluidRow(
            shinydashboard::box(
              title = "CpG Loci in Signature",
              status = "primary", solidHeader = TRUE, width = 12,
              p("Genomic loci in this signature with their mean methylation per group."),
              DT::DTOutput("cpg_detail_table")
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

      perf_val <- "N/A"
      if (!is.null(res$svm_model)) {
        val <- res$svm_model$results[[res$metric_used]]
        perf_val <- round(max(val, na.rm = TRUE), 3)
      }

      acc_str <- if (!is.null(res$accuracy) && !is.na(res$accuracy)) {
        sprintf("%.1f%%", res$accuracy * 100)
      } else { "N/A" }

      tagList(
        p(style = "color:#b8c7ce;",
          strong("CpGs: "), res$n_cpgs),
        p(style = "color:#b8c7ce;",
          strong("Metric: "), metric_lbl),
        p(style = "color:#b8c7ce;",
          strong("Performance: "), perf_val),
        p(style = "color:#b8c7ce;",
          strong("Accuracy: "), acc_str),
        p(style = "color:#b8c7ce;",
          strong("Silhouette: "), round(res$silhouette, 3))
      )
    })

    # ── Value boxes ──────────────────────────────────────────────────────
    output$vbox_ncpgs <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(
        curr()$n_cpgs, "CpGs in Signature",
        icon = icon("dna"), color = "purple")
    })

    output$vbox_metric1 <- shinydashboard::renderValueBox({
      res <- curr()
      if (is.null(res$svm_model)) {
        shinydashboard::valueBox(
          "FAILED", "SVM Result",
          icon = icon("times-circle"), color = "red")
      } else {
        lbl <- if (res$metric_used == "ROC") "AUC (ROC)" else "Accuracy"
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
          "N/A", "Secondary Metric",
          icon = icon("question"), color = "yellow")
      } else {
        if (!is.null(res$sensitivity) && !is.na(res$sensitivity)) {
          shinydashboard::valueBox(
            sprintf("%.1f%%", res$sensitivity * 100), "Sensitivity",
            icon = icon("heartbeat"), color = "green")
        } else if ("Kappa" %in% colnames(res$svm_model$results)) {
          val <- round(max(res$svm_model$results$Kappa, na.rm = TRUE), 3)
          shinydashboard::valueBox(
            val, "Kappa",
            icon = icon("chart-bar"), color = "green")
        } else {
          shinydashboard::valueBox(
            "N/A", "Secondary",
            icon = icon("question"), color = "yellow")
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

    # ── Silhouette interpretation ────────────────────────────────────────
    output$sil_interpretation_ui <- renderUI({
      sil <- curr()$silhouette
      interp <- if (sil > 0.7) list("Strong structure", "green", "check-circle")
      else if (sil > 0.5) list("Reasonable structure", "blue", "info-circle")
      else if (sil > 0.25) list("Weak structure", "orange", "exclamation-triangle")
      else list("No substantial structure", "red", "times-circle")

      div(style = sprintf("padding:15px; border-left:4px solid %s; background:#f9f9f9; border-radius:4px;",
                          interp[[2]]),
          icon(interp[[3]], style = sprintf("color:%s;", interp[[2]])),
          strong(sprintf(" Silhouette = %.3f", sil)),
          br(), br(),
          p(interp[[1]]),
          p(tags$small("Range: -1 to 1. Values > 0.5 indicate good cluster separation.",
                       "Values > 0.7 indicate strong separation."))
      )
    })

    # ── Confusion matrix ─────────────────────────────────────────────────
    output$conf_matrix_ui <- renderUI({
      res <- curr()
      if (!is.null(res$confusion)) {
        tagList(
          tableOutput("conf_table"),
          if (!is.null(res$accuracy) && !is.na(res$accuracy)) {
            p(strong("Accuracy: "), sprintf("%.1f%%", res$accuracy * 100))
          },
          if (!is.null(res$sensitivity) && !is.na(res$sensitivity)) {
            p(strong("Sensitivity: "), sprintf("%.1f%%", res$sensitivity * 100))
          },
          if (!is.null(res$specificity) && !is.na(res$specificity)) {
            p(strong("Specificity: "), sprintf("%.1f%%", res$specificity * 100))
          }
        )
      } else {
        div(class = "alert alert-info",
            icon("info-circle"), " Confusion matrix not available.")
      }
    })

    output$conf_table <- renderTable({
      res <- curr()
      if (!is.null(res$confusion)) {
        as.data.frame.matrix(res$confusion)
      }
    }, rownames = TRUE)

    # ── Benchmark table ──────────────────────────────────────────────────
    output$benchmark_table <- DT::renderDT({
      rows <- lapply(names(validation_results), function(nm) {
        r <- validation_results[[nm]]
        perf <- if (!is.null(r$svm_model)) {
          round(max(r$svm_model$results[[r$metric_used]], na.rm = TRUE), 3)
        } else { NA }

        acc <- if (!is.null(r$accuracy) && !is.na(r$accuracy)) {
          sprintf("%.1f%%", r$accuracy * 100)
        } else { "N/A" }

        data.frame(
          Signature   = nm,
          CpGs        = r$n_cpgs,
          Metric      = if (!is.null(r$svm_model)) r$metric_used else "N/A",
          Performance = perf,
          Accuracy    = acc,
          Silhouette  = round(r$silhouette, 3),
          stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, rows)
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 20,
                                   order = list(list(3, "desc")))) |>
        DT::formatStyle("Silhouette",
                        backgroundColor = DT::styleInterval(
                          c(0.25, 0.5),
                          c("#fadbd8", "#fef9e7", "#d5f5e3")))
    })

    # ── PCA ──────────────────────────────────────────────────────────────
    output$pcaPlot <- plotly::renderPlotly({
      res <- curr()
      if (is.null(res$pca)) {
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "PCA unavailable."))
      }

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
        ggplot2::labs(
          x = paste0("PC1 (", ve[1], "%)"),
          y = paste0("PC2 (", ve[2], "%)"),
          title = paste0("PCA - ", input$selected_sig)) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_manual(values = c("CTL" = "#3498db",
                                               "Test" = "#e74c3c"))

      # Only add ellipses for groups with >= 4 samples
      grp_n <- table(df$Group)
      ok_grps <- names(grp_n[grp_n >= 4])
      if (length(ok_grps) > 0) {
        tryCatch({
          p <- p + ggplot2::stat_ellipse(
            level = 0.95, linetype = "dashed", alpha = 0.4,
            data = df[df$Group %in% ok_grps, ])
        }, error = function(e) NULL)
      }

      plotly::ggplotly(p, tooltip = "text")
    })

    output$screePlot <- renderPlot({
      res <- curr()
      if (is.null(res$pca)) return(NULL)

      ve <- summary(res$pca)$importance[2, ] * 100
      n_show <- min(10, length(ve))
      df <- data.frame(PC = seq_len(n_show), Var = ve[seq_len(n_show)])

      ggplot2::ggplot(df, ggplot2::aes(PC, Var)) +
        ggplot2::geom_col(fill = "#2980b9", alpha = 0.8) +
        ggplot2::geom_line(colour = "#e74c3c", linewidth = 1) +
        ggplot2::geom_point(colour = "#e74c3c", size = 2) +
        ggplot2::scale_x_continuous(breaks = df$PC) +
        ggplot2::labs(title = "Scree Plot", x = "PC", y = "Variance (%)") +
        ggplot2::theme_minimal(base_size = 12)
    })

    output$pca_var_table <- DT::renderDT({
      res <- curr()
      if (is.null(res$pca)) return(NULL)

      imp <- summary(res$pca)$importance
      n_show <- min(10, ncol(imp))
      df <- data.frame(
        PC = paste0("PC", seq_len(n_show)),
        StdDev = round(imp[1, seq_len(n_show)], 3),
        PropVar = sprintf("%.1f%%", imp[2, seq_len(n_show)] * 100),
        CumVar = sprintf("%.1f%%", imp[3, seq_len(n_show)] * 100),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 10, dom = "t"))
    })

    # ── Heatmap ──────────────────────────────────────────────────────────
    output$heatmapPlot <- plotly::renderPlotly({
      res <- curr()
      if (is.null(res$beta_matrix)) {
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "No data."))
      }

      bm <- res$beta_matrix
      if (nrow(bm) > 200) {
        rv <- apply(bm, 1, var, na.rm = TRUE)
        bm <- bm[order(rv, decreasing = TRUE)[1:200], ]
      }

      annot_col <- data.frame(Group = res$raw_data$Class)
      rownames(annot_col) <- rownames(res$raw_data)

      heatmaply::heatmaply(
        bm,
        limits          = c(0, 1),
        colors          = grDevices::colorRampPalette(
          c("#2980b9", "white", "#c0392b"))(100),
        col_side_colors = annot_col,
        row_dend_left   = TRUE,
        col_dend_top    = TRUE,
        branches_lwd    = 0.2,
        labRow = if (nrow(bm) > 50) NULL else rownames(bm),
        main   = paste0("Methylation - ", res$n_cpgs, " CpGs"),
        xlab   = "Samples",
        ylab   = "CpG Loci"
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
          "Please rely on PCA and Heatmap tabs."))
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
      if (is.null(preds)) return(NULL)

      if ("C" %in% colnames(preds) && !is.null(res$svm_model$bestTune$C)) {
        preds <- preds[preds$C == res$svm_model$bestTune$C, ]
      }

      if (nrow(preds) == 0) return(NULL)

      cm    <- caret::confusionMatrix(preds$pred, preds$obs)
      df_cm <- as.data.frame(cm$table)

      ggplot2::ggplot(df_cm,
                      ggplot2::aes(x = Reference, y = Prediction,
                                   fill = Freq)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = Freq),
                           size = 8, color = "white") +
        ggplot2::scale_fill_gradient(low = "#868e96",
                                     high = "#2980b9") +
        ggplot2::theme_minimal(base_size = 14) +
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
            "Accuracy mode was used (insufficient samples for class probabilities).")
      }
    })

    output$rocPlot <- renderPlot({
      res <- curr()
      if (res$metric_used != "ROC" || is.null(res$svm_model)) return(NULL)

      preds <- res$svm_model$pred
      if (is.null(preds)) return(NULL)

      if ("C" %in% colnames(preds) && !is.null(res$svm_model$bestTune$C)) {
        preds <- preds[preds$C == res$svm_model$bestTune$C, ]
      }

      lvls <- levels(preds$obs)
      if (length(lvls) < 2 || !lvls[1] %in% colnames(preds)) return(NULL)

      roc_obj <- tryCatch(
        pROC::roc(preds$obs, preds[[lvls[1]]]),
        error = function(e) NULL)

      if (is.null(roc_obj)) return(NULL)

      plot(roc_obj,
           print.auc       = TRUE,
           auc.polygon     = TRUE,
           grid             = c(0.1, 0.2),
           grid.col         = c("green", "red"),
           max.auc.polygon  = TRUE,
           auc.polygon.col  = "lightblue",
           print.thres      = TRUE,
           legacy.axes      = TRUE)
    })

    # ── Variable Importance ──────────────────────────────────────────────
    output$svm_importance_ui <- renderUI({
      res <- curr()
      if (is.null(res$svm_model)) {
        return(div(
          class = "alert alert-danger",
          style = "padding:15px;",
          icon("exclamation-triangle"),
          strong(" Requires a valid SVM model."),
          br(), "SVM failed for this signature."))
      }
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

      imp$Score <- if ("Overall" %in% colnames(imp)) imp$Overall
      else apply(imp, 1, max)

      imp$RawName <- rownames(imp)
      known <- names(res$genes_map)

      recover_name <- function(x) {
        if (x %in% known) return(x)
        x2 <- gsub("\\.", ":", x, fixed = FALSE)
        x2 <- sub(":([0-9]+):([0-9]+)$", ":\\1-\\2", x2)
        if (x2 %in% known) return(x2)
        x3 <- gsub("^X", "", x)
        x3 <- gsub("\\.", ":", x3)
        x3 <- sub(":([0-9]+):([0-9]+)$", ":\\1-\\2", x3)
        if (x3 %in% known) return(x3)
        return(x)
      }

      imp$Locus <- vapply(imp$RawName, recover_name, character(1))

      imp <- imp[order(imp$Score, decreasing = TRUE), ]
      plot_data <- head(imp, 20L)
      plot_data$Locus <- factor(plot_data$Locus,
                                levels = rev(plot_data$Locus))

      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = Locus, y = Score)) +
        ggplot2::geom_bar(stat = "identity", fill = "#c0392b",
                          alpha = 0.8) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Top 20 Discriminative Loci",
                      x = "", y = "Variable Importance") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

      plotly::ggplotly(p)
    })

    # ── CpG Detail table ─────────────────────────────────────────────────
    output$cpg_detail_table <- DT::renderDT({
      res <- curr()
      if (is.null(res$beta_matrix)) return(NULL)

      bm <- res$beta_matrix
      groups <- res$raw_data$Class

      grp_levels <- levels(groups)
      mean_df <- data.frame(
        Locus = rownames(bm),
        stringsAsFactors = FALSE
      )

      for (g in grp_levels) {
        cols <- which(groups == g)
        if (length(cols) > 1) {
          mean_df[[paste0("Mean_", g)]] <- round(rowMeans(bm[, cols, drop = FALSE], na.rm = TRUE), 4)
          mean_df[[paste0("SD_", g)]]   <- round(apply(bm[, cols, drop = FALSE], 1, sd, na.rm = TRUE), 4)
        } else {
          mean_df[[paste0("Mean_", g)]] <- round(bm[, cols], 4)
          mean_df[[paste0("SD_", g)]]   <- NA
        }
      }

      if (length(grp_levels) == 2) {
        mean_df$Delta <- round(
          mean_df[[paste0("Mean_", grp_levels[2])]] -
            mean_df[[paste0("Mean_", grp_levels[1])]],
          4)
      }

      parts <- strsplit(mean_df$Locus, "[:|-]")
      mean_df$Chr   <- sapply(parts, function(x) x[1])
      mean_df$Start <- as.integer(sapply(parts, function(x) if (length(x) >= 2) x[2] else NA))

      first_cols <- c("Locus", "Chr", "Start")
      other_cols <- setdiff(colnames(mean_df), first_cols)
      mean_df <- mean_df[, c(first_cols, other_cols)]

      DT::datatable(mean_df, rownames = FALSE,
                    options = list(pageLength = 20, scrollX = TRUE,
                                   order = list(list(which(colnames(mean_df) == "Delta") - 1, "desc")))) |>
        DT::formatStyle("Delta",
                        backgroundColor = DT::styleInterval(
                          c(-0.1, 0.1),
                          c("#aed6f1", "white", "#f5b7b1")))
    })

    # ── Close ────────────────────────────────────────────────────────────
    observeEvent(input$close_val_btn, shiny::stopApp())
    session$onSessionEnded(function() shiny::stopApp())
  }

  message("[mETHYLotest] Starting NGS Signature Validator UI...")
  invisible(shiny::runApp(shiny::shinyApp(ui, server),
                          host = "0.0.0.0",
                          launch.browser = TRUE))
}
