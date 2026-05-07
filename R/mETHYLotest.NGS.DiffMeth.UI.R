#' Interactive Differential Methylation Results UI
#'
#' @param diff_list Named list of methylDiff objects.
#' @param output_dir Directory to save filtered results.
#' @return NULL (invisibly). Blocks until closed.
#' @export
mETHYLotest.NGS.DiffMeth.UI <- function(diff_list, output_dir) {

  for (pkg in c("shiny", "shinydashboard", "plotly", "DT", "ggplot2", "writexl")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("[mETHYLotest] Package '%s' is required.", pkg))
    }
  }

  ui <- shinydashboard::dashboardPage(
    skin = "blue",
    shinydashboard::dashboardHeader(title = "mETHYLotest DiffMeth Explorer"),

    shinydashboard::dashboardSidebar(
      width = 300,
      div(style = "padding:15px;",
          h4(icon("layer-group"), " Model Selection"),
          selectInput("model_select", "Scenario:",
                      choices = names(diff_list)),
          hr(),
          h4(icon("filter"), " Filters"),
          sliderInput("qvalue_cutoff", "Max q-value (FDR):",
                      min = 0.001, max = 0.1,
                      value = 0.05, step = 0.001),
          sliderInput("diff_cutoff", "Min Meth Diff (%):",
                      min = 0, max = 100,
                      value = 10, step = 5),
          selectInput("type_select", "Type:",
                      choices = c("All" = "all",
                                  "Hyper" = "hyper",
                                  "Hypo" = "hypo")),
          hr(),
          h4(icon("pen"), " Create Signature"),
          textInput("sig_name", "Signature name:", value = ""),
          numericInput("sig_max_cpgs", "Max CpGs:", value = 500,
                       min = 10, max = 5000, step = 50),
          actionButton("btn_save_sig", "Save as Signature",
                       icon = icon("file-signature"),
                       class = "btn-warning",
                       style = "width:100%;"),
          hr(),
          h4(icon("save"), " Export"),
          p(style = "font-size:0.8em; color:#ccc;",
            paste("Dest:", basename(output_dir))),
          actionButton("btn_save_filtered", "Save Filtered Results",
                       icon = icon("hdd"),
                       class = "btn-success",
                       style = "width:100%;"),
          hr(),
          actionButton("btn_close", "Close & Continue",
                       icon = icon("power-off"),
                       class = "btn-danger",
                       style = "width:100%;")
      )
    ),

    shinydashboard::dashboardBody(
      fluidRow(
        shinydashboard::valueBoxOutput("box_tested", width = 3),
        shinydashboard::valueBoxOutput("box_total",  width = 3),
        shinydashboard::valueBoxOutput("box_hyper",  width = 3),
        shinydashboard::valueBoxOutput("box_hypo",   width = 3)
      ),
      fluidRow(
        shinydashboard::tabBox(
          width = 12, title = "Results", id = "tabs",
          tabPanel("Volcano", icon = icon("chart-area"),
                   plotly::plotlyOutput("volcano", height = "600px")),
          tabPanel("Distribution", icon = icon("chart-bar"),
                   plotly::plotlyOutput("chr_plot", height = "400px"),
                   plotly::plotlyOutput("diff_hist", height = "400px")),
          tabPanel("Table", icon = icon("table"),
                   uiOutput("table_status"),
                   DT::DTOutput("table_dmc"))
        )
      )
    )
  )

  server <- function(input, output, session) {

    current_obj <- reactive({
      req(input$model_select)
      diff_list[[input$model_select]]
    })

    # ── Full data (all positions, no filter) ──
    full_data <- reactive({
      obj <- current_obj()
      req(obj)
      df <- methylKit::getData(obj)
      df <- df[!is.na(df$qvalue) & !is.na(df$meth.diff), ]
      df
    })

    # ── Filtered data (user thresholds) ──
    filtered_data <- reactive({
      obj <- current_obj()
      req(obj)

      res <- tryCatch(
        methylKit::getMethylDiff(
          obj,
          difference = input$diff_cutoff,
          qvalue     = input$qvalue_cutoff,
          type       = input$type_select),
        error = function(e) NULL)

      if (is.null(res)) return(data.frame())
      df <- methylKit::getData(res)
      if (is.null(df) || nrow(df) == 0L) return(data.frame())
      df
    })

    # ── Value boxes ──
    output$box_tested <- shinydashboard::renderValueBox({
      n <- nrow(full_data())
      shinydashboard::valueBox(
        format(n, big.mark = ","), "Total Tested",
        icon = icon("microscope"), color = "light-blue")
    })

    output$box_total <- shinydashboard::renderValueBox({
      n <- nrow(filtered_data())
      shinydashboard::valueBox(
        format(n, big.mark = ","), "DMCs (filtered)",
        icon = icon("dna"),
        color = if (n > 0) "purple" else "black")
    })

    output$box_hyper <- shinydashboard::renderValueBox({
      df <- filtered_data()
      n <- if (nrow(df) > 0L) sum(df$meth.diff > 0, na.rm = TRUE) else 0L
      shinydashboard::valueBox(
        format(n, big.mark = ","), "Hyper-methylated",
        icon = icon("arrow-up"), color = "red")
    })

    output$box_hypo <- shinydashboard::renderValueBox({
      df <- filtered_data()
      n <- if (nrow(df) > 0L) sum(df$meth.diff < 0, na.rm = TRUE) else 0L
      shinydashboard::valueBox(
        format(n, big.mark = ","), "Hypo-methylated",
        icon = icon("arrow-down"), color = "blue")
    })

    # ── Volcano ──
    output$volcano <- plotly::renderPlotly({
      fd <- full_data()
      if (nrow(fd) == 0L) {
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "No data available."))
      }

      # Handle q=0
      min_q <- min(fd$qvalue[fd$qvalue > 0], na.rm = TRUE)
      if (is.infinite(min_q)) min_q <- 1e-300
      fd$qvalue[fd$qvalue == 0] <- min_q

      # Status
      fd$Status <- "NS"
      sig <- fd$qvalue < input$qvalue_cutoff &
        abs(fd$meth.diff) >= input$diff_cutoff
      fd$Status[sig & fd$meth.diff > 0] <- "Hyper"
      fd$Status[sig & fd$meth.diff < 0] <- "Hypo"

      fd$logQ <- -log10(fd$qvalue)

      # Downsample NS for speed
      df_sig <- fd[fd$Status != "NS", ]
      df_ns  <- fd[fd$Status == "NS", ]
      if (nrow(df_ns) > 10000) {
        df_ns <- df_ns[sample(nrow(df_ns), 10000), ]
      }
      viz <- rbind(df_sig, df_ns)

      p <- ggplot2::ggplot(
        viz, ggplot2::aes(
          x = meth.diff, y = logQ, color = Status,
          text = paste0("Chr: ", chr,
                        "<br>Pos: ", start,
                        "<br>Diff: ", round(meth.diff, 2), "%",
                        "<br>Q: ", signif(qvalue, 3)))
      ) +
        ggplot2::geom_point(alpha = 0.6, size = 1.2) +
        ggplot2::scale_color_manual(
          values = c("Hyper" = "#c0392b",
                     "Hypo" = "#2980b9",
                     "NS" = "grey70")) +
        ggplot2::geom_vline(
          xintercept = c(-input$diff_cutoff, input$diff_cutoff),
          linetype = "dashed", colour = "grey40") +
        ggplot2::geom_hline(
          yintercept = -log10(input$qvalue_cutoff),
          linetype = "dashed", colour = "grey40") +
        ggplot2::labs(
          title = sprintf("%s - %d Hyper / %d Hypo",
                          input$model_select,
                          sum(viz$Status == "Hyper"),
                          sum(viz$Status == "Hypo")),
          x = "Methylation Diff (%)",
          y = "-log10(q-value)") +
        ggplot2::theme_minimal()

      plotly::ggplotly(p, tooltip = "text")
    })

    # ── Chromosome distribution ──
    output$chr_plot <- plotly::renderPlotly({
      df <- filtered_data()
      if (nrow(df) == 0L) {
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "No DMCs to show."))
      }

      df$Direction <- ifelse(df$meth.diff > 0, "Hyper", "Hypo")
      df$chr_clean <- gsub("chr", "", df$chr, ignore.case = TRUE)
      df$chr_f <- factor(df$chr_clean,
                         levels = c(as.character(1:22), "X", "Y", "M"))

      chr_counts <- as.data.frame(table(df$chr_f, df$Direction),
                                  stringsAsFactors = FALSE)
      colnames(chr_counts) <- c("Chr", "Direction", "N")
      chr_counts <- chr_counts[chr_counts$N > 0, ]

      p <- ggplot2::ggplot(chr_counts,
                           ggplot2::aes(Chr, N, fill = Direction)) +
        ggplot2::geom_col(position = "dodge", alpha = 0.85) +
        ggplot2::scale_fill_manual(values = c(Hyper = "#c0392b",
                                              Hypo = "#2980b9")) +
        ggplot2::labs(title = "DMC Chromosome Distribution",
                      x = "Chromosome", y = "Count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       legend.position = "bottom")

      plotly::ggplotly(p)
    })

    # ── meth.diff histogram ──
    output$diff_hist <- plotly::renderPlotly({
      df <- filtered_data()
      if (nrow(df) == 0L) {
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "No DMCs to show."))
      }

      df$Direction <- ifelse(df$meth.diff > 0, "Hyper", "Hypo")

      p <- ggplot2::ggplot(df, ggplot2::aes(meth.diff, fill = Direction)) +
        ggplot2::geom_histogram(bins = 50, alpha = 0.7) +
        ggplot2::scale_fill_manual(values = c(Hyper = "#c0392b",
                                              Hypo = "#2980b9")) +
        ggplot2::geom_vline(xintercept = c(-input$diff_cutoff, input$diff_cutoff),
                            linetype = "dashed") +
        ggplot2::labs(title = "Methylation Difference Distribution",
                      x = "Methylation Diff (%)", y = "Count") +
        ggplot2::theme_minimal()

      plotly::ggplotly(p)
    })

    # ── Table status ──
    output$table_status <- renderUI({
      df <- filtered_data()
      if (nrow(df) == 0L) {
        div(class = "alert alert-warning",
            style = "margin:10px;",
            icon("exclamation-triangle"),
            strong(" No DMCs pass current filters."),
            br(),
            "Try relaxing the q-value or difference threshold.")
      }
    })

    # ── Table ──
    output$table_dmc <- DT::renderDT({
      df <- filtered_data()
      if (nrow(df) == 0L) {
        return(DT::datatable(
          data.frame(Message = "No DMCs found with current filters."),
          options = list(dom = "t")))
      }

      # Add locus ID column
      df$Locus <- paste0(df$chr, ":", df$start, "-", df$end)

      # Round numeric columns
      df$meth.diff <- round(df$meth.diff, 2)
      df$qvalue    <- signif(df$qvalue, 3)
      df$pvalue    <- signif(df$pvalue, 3)

      DT::datatable(
        df, extensions = "Buttons",
        options = list(
          pageLength = 15, scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv"),
          order = list(list(which(colnames(df) == "qvalue") - 1, "asc"))))
    })

    # ── Save filtered results ──
    observeEvent(input$btn_save_filtered, {
      df <- filtered_data()

      if (nrow(df) == 0L) {
        showNotification("No DMCs to save. Adjust your filters.",
                         type = "warning")
        return()
      }

      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }

      safe <- gsub("[^A-Za-z0-9_.-]", "_", input$model_select)
      base <- paste0("Filtered_", safe,
                     "_q", input$qvalue_cutoff,
                     "_diff", input$diff_cutoff,
                     "_", input$type_select)

      # Excel
      xlsx_path <- file.path(output_dir, paste0(base, ".xlsx"))
      tryCatch({
        writexl::write_xlsx(df, xlsx_path)
      }, error = function(e) {
        showNotification(paste("Excel error:", e$message), type = "error")
      })

      # BED
      bed_path <- file.path(output_dir, paste0(base, ".bed"))
      bed_df <- data.frame(
        chrom      = df$chr,
        chromStart = df$start - 1L,
        chromEnd   = df$end,
        name       = paste0("DMC_", seq_len(nrow(df))),
        score      = pmin(round(abs(df$meth.diff) * 10), 1000L),
        strand     = df$strand
      )

      tryCatch({
        write.table(bed_df, bed_path, quote = FALSE,
                    sep = "\t", row.names = FALSE,
                    col.names = FALSE)
        showNotification(
          paste0("Saved ", nrow(df), " DMCs: ",
                 basename(xlsx_path), " & ", basename(bed_path)),
          type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("BED error:", e$message), type = "error")
      })
    })

    # ── Save as signature ──
    observeEvent(input$btn_save_sig, {
      df <- filtered_data()

      if (nrow(df) == 0L) {
        showNotification("No DMCs to create a signature. Adjust filters.",
                         type = "warning")
        return()
      }

      sig_name <- trimws(input$sig_name)
      if (!nzchar(sig_name)) {
        sig_name <- paste0("User_",
                           gsub("[^A-Za-z0-9_]", "_", input$model_select),
                           "_q", input$qvalue_cutoff,
                           "_d", input$diff_cutoff)
      }

      # Score and rank
      df$score <- -log10(pmax(df$qvalue, 1e-300)) * abs(df$meth.diff)
      df <- df[order(df$score, decreasing = TRUE), ]

      max_cpgs <- input$sig_max_cpgs
      if (nrow(df) > max_cpgs) {
        df <- df[seq_len(max_cpgs), ]
      }

      # Write to signatures directory (parent of output_dir)
      sig_dir <- file.path(dirname(output_dir), "Signatures")
      if (!dir.exists(sig_dir)) dir.create(sig_dir, recursive = TRUE)

      sig_ids <- paste0(df$chr, ":", df$start, "-", df$end)
      safe_name <- gsub("[^A-Za-z0-9_.-]", "_", sig_name)
      sig_path <- file.path(sig_dir, paste0(safe_name, ".txt"))
      writeLines(sig_ids, sig_path)

      showNotification(
        paste0("Signature saved: ", basename(sig_path),
               " (", length(sig_ids), " loci)"),
        type = "message", duration = 8)
    })

    # ── Close button ──
    observeEvent(input$btn_close, {
      message("[mETHYLotest] DiffMeth Explorer closed by user.")
      shiny::stopApp()
    })

    # ── Auto-close on browser disconnect ──
    session$onSessionEnded(function() {
      message("[mETHYLotest] Browser disconnected. Stopping DiffMeth Explorer.")
      shiny::stopApp()
    })
  }

  message("[mETHYLotest] Starting DiffMeth Explorer...")
  message("[mETHYLotest] Click 'Close & Continue' or close the browser tab to proceed.")

  invisible(shiny::runApp(
    shiny::shinyApp(ui, server),
    launch.browser = TRUE,
    quiet = TRUE
  ))
}
