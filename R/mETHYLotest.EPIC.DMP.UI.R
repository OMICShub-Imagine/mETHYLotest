#' mETHYLotest EPIC DMP Explorer
#'
#' @description
#' Shiny Dashboard for interactive DMP exploration.
#' Style matches \code{mETHYLotest.EPIC.ProjectUI}: skin = "purple".
#'
#' Features: Volcano, Density, ChromoMap, Context distributions,
#' Genomic Window, Interactive Heatmap, Signature export.
#'
#' @param dmp_results        Named list of DMP data frames (ChAMP output).
#' @param beta               Numeric matrix (probes x samples).
#' @param pheno              Phenotype vector (optional).
#' @param ages               Age vector (optional).
#' @param arraytype          \code{"EPICv1"}, \code{"EPICv2"}, or \code{"450K"}.
#' @param path_to_episignature Output directory for signatures.
#'
#' @return Path to episignature directory (invisibly).
#' @export
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import plotly
#' @import ggplot2
#' @importFrom DT DTOutput renderDT datatable dataTableProxy selectRows
#' @importFrom heatmaply heatmaply
#' @importFrom chromoMap chromoMap chromoMapOutput renderChromoMap
#' @importFrom tidyr pivot_longer
#' @importFrom scales percent_format comma
#' @importFrom utils combn
mETHYLotest.EPIC.DMP.UI <- function(dmp_results,
                                   beta,
                                   pheno                = NULL,
                                   ages                 = NULL,
                                   arraytype            = "EPICv1",
                                   path_to_episignature = file.path(
                                     getwd(), "EPIC", "Episignatures")) {

  for (pkg in c("shiny", "shinydashboard", "shinyjs", "plotly", "DT",
                "ggplot2", "heatmaply", "chromoMap", "tidyr", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", pkg))
  }

  if (missing(beta))
    stop("[mETHYLotest DMP UI] 'beta' matrix is required.")
  if (!is.null(pheno) && ncol(beta) != length(pheno))
    warning("[mETHYLotest DMP UI] Dimension mismatch pheno/beta.")
  if (!is.null(ages) && ncol(beta) != length(ages))
    warning("[mETHYLotest DMP UI] Dimension mismatch ages/beta.")

  # ── Annotation loading ────────────────────────────────────────────────────
  # Use the same Locations dataset approach proven in sex prediction
  # (ChAMPdata AnnoEPICv2$Annotation has no chr column)
  probe_features <- NULL

  anno_pkg <- if (grepl("EPICv2", arraytype, ignore.case = TRUE))
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
  else if (grepl("EPIC", arraytype, ignore.case = TRUE))
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  else if (grepl("mouse", arraytype, ignore.case = TRUE))
    "IlluminaMouseMethylationmanifest"
  else
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"

  if (requireNamespace(anno_pkg, quietly = TRUE)) {
    local_env <- new.env(parent = emptyenv())
    tryCatch({
      utils::data("Locations", package = anno_pkg, envir = local_env)
      locs <- local_env[["Locations"]]
      if (!is.null(locs) && "chr" %in% colnames(locs)) {
        probe_features <- data.frame(
          CpG       = rownames(locs),
          CHR       = as.character(locs[["chr"]]),
          MAPINFO   = as.numeric(locs[["pos"]]),
          CHR_Clean = gsub("chr", "", as.character(locs[["chr"]]),
                           ignore.case = TRUE),
          stringsAsFactors = FALSE
        )
        message("[mETHYLotest DMP UI] Annotation: ", anno_pkg,
                " (", nrow(probe_features), " probes)")
      }
    }, error = function(e)
      message("[mETHYLotest DMP UI] Annotation load failed: ", e$message))
  }

  # Fallback: extract chr/pos from DMP results themselves
  if (is.null(probe_features)) {
    message("[mETHYLotest DMP UI] Using annotation from DMP results (fallback).")
    combined <- do.call(rbind, lapply(dmp_results, function(d) {
      cols <- intersect(c("CHR", "MAPINFO"), colnames(d))
      if (length(cols) == 2L) d[, cols, drop = FALSE] else NULL
    }))
    if (!is.null(combined)) {
      combined <- combined[!duplicated(rownames(combined)), ]
      probe_features <- data.frame(
        CpG       = rownames(combined),
        CHR       = as.character(combined$CHR),
        MAPINFO   = as.numeric(combined$MAPINFO),
        CHR_Clean = gsub("chr", "", as.character(combined$CHR),
                         ignore.case = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── ChromoMap chromosome reference (hg38) ─────────────────────────────────
  chr_ref <- data.frame(
    chrom = c(paste0("chr", 1:22), "chrX", "chrY"),
    start = 1L,
    end   = c(248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
              170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
              135086622L, 133275309L, 114364328L, 107043718L, 101991189L,
              90338345L,  83257441L,  80373285L,  58617616L,  64444167L,
              46709983L,  50818468L,  156040895L, 57227415L),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------------------------
  # UI — shinydashboard skin purple
  # --------------------------------------------------------------------------
  ui <- shinydashboard::dashboardPage(
    skin = "purple",

    shinydashboard::dashboardHeader(title = "mETHYLotest DMP Explorer"),

    shinydashboard::dashboardSidebar(
      width = 280,

      div(style = "padding:15px;",

          h5("Analysis Settings", style = "color:#b8c7ce;"),

          selectInput("selected_comparison", "Comparison:",
                      choices = names(dmp_results)),

          hr(style = "border-color:#444;"),

          h5("Thresholds", style = "color:#b8c7ce;"),
          sliderInput("pval_cutoff", "Adj. P-value <",
                      min = 0, max = 0.2, value = 0.05, step = 0.001),
          sliderInput("raw_pval_cutoff", "Raw P-value <",
                      min = 0, max = 0.2, value = 0.05, step = 0.001),
          sliderInput("deltabeta_cutoff", "|Delta Beta| >",
                      min = 0, max = 1, value = 0.1, step = 0.01),

          hr(style = "border-color:#444;"),

          h5("Summary", style = "color:#b8c7ce;"),
          uiOutput("summary_callout")
      ),

      shinydashboard::sidebarMenu(
        id = "main_tabs",
        shinydashboard::menuItem("Volcano & Density", tabName = "tab_volcano",
                                 icon = icon("chart-bar")),
        shinydashboard::menuItem("Chromosomes",       tabName = "tab_chromo",
                                 icon = icon("dna")),
        shinydashboard::menuItem("Distributions",     tabName = "tab_dist",
                                 icon = icon("th")),
        shinydashboard::menuItem("Genomic Window",    tabName = "tab_window",
                                 icon = icon("search-location")),
        shinydashboard::menuItem("Table & Signatures", tabName = "tab_table",
                                 icon = icon("table"),
                                 badgeLabel = "Export",
                                 badgeColor = "green"),
        shinydashboard::menuItem("Heatmap",           tabName = "tab_heatmap",
                                 icon = icon("border-all"))
      ),

      div(style = "padding:15px; margin-top:10px;",
          actionButton("close_btn", " Close & Exit",
                       icon  = icon("power-off"),
                       class = "btn-danger btn-block"))
    ),

    shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),

      shinydashboard::tabItems(

        # ── Volcano & Density ──────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tab_volcano",
                                fluidRow(
                                  shinydashboard::box(
                                    title = "Volcano Plot", status = "primary",
                                    solidHeader = TRUE, width = 7,
                                    plotly::plotlyOutput("volcanoPlot", height = "500px")
                                  ),
                                  shinydashboard::box(
                                    title = "Delta Beta Density", status = "info",
                                    solidHeader = TRUE, width = 5,
                                    plotly::plotlyOutput("plot_density", height = "500px")
                                  )
                                )
        ),

        # ── Chromosomes ────────────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tab_chromo",
                                fluidRow(
                                  shinydashboard::box(
                                    title = "Interactive Chromosomal Map (Filtered DMPs)",
                                    status = "primary", solidHeader = TRUE, width = 12,
                                    chromoMap::chromoMapOutput("karyoPlot", height = "600px")
                                  )
                                )
        ),

        # ── Distributions ──────────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tab_dist",
                                fluidRow(
                                  shinydashboard::box(
                                    title = "CGI Context", status = "warning",
                                    solidHeader = TRUE, width = 4,
                                    plotly::plotlyOutput("plot_cgi", height = "350px")
                                  ),
                                  shinydashboard::box(
                                    title = "Genomic Feature", status = "warning",
                                    solidHeader = TRUE, width = 4,
                                    plotly::plotlyOutput("plot_feature", height = "350px")
                                  ),
                                  shinydashboard::box(
                                    title = "Chromosomal Distribution", status = "warning",
                                    solidHeader = TRUE, width = 4,
                                    plotly::plotlyOutput("plot_chr", height = "350px")
                                  )
                                ),
                                fluidRow(
                                  shinydashboard::box(
                                    title = "Feature x CGI (Significant DMPs)",
                                    status = "info", solidHeader = TRUE, width = 12,
                                    plotly::plotlyOutput("plot_feature_cgi", height = "450px")
                                  )
                                )
        ),

        # ── Genomic Window ─────────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tab_window",
                                fluidRow(
                                  shinydashboard::box(
                                    title = "Region Viewer", status = "primary",
                                    solidHeader = TRUE, width = 12,
                                    fluidRow(
                                      column(4, selectizeInput("target_cpg_selector", "Target CpG:",
                                                               choices = NULL, multiple = FALSE,
                                                               options = list(
                                                                 placeholder = "Type CpG ID..."))),
                                      column(4, sliderInput("window_bp",
                                                            "Window (+/- bp):",
                                                            min = 1000, max = 50000,
                                                            value = 5000, step = 1000)),
                                      column(4, selectizeInput("sample_filter",
                                                               "Filter Samples:",
                                                               choices = NULL, multiple = TRUE,
                                                               options = list(
                                                                 placeholder = "All")))
                                    ),
                                    if (!is.null(ages))
                                      checkboxInput("hide_missing_age",
                                                    "Exclude samples with missing Age",
                                                    value = FALSE),
                                    verbatimTextOutput("selected_cpg_info"),
                                    plotly::plotlyOutput("windowPlot", height = "600px")
                                  )
                                )
        ),

        # ── Table & Signatures ─────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tab_table",
                                fluidRow(
                                  shinydashboard::box(
                                    title = "Signature Manager", status = "success",
                                    solidHeader = TRUE, width = 12,
                                    fluidRow(
                                      column(4,
                                             textInput("sig_name", "Signature Name:",
                                                       placeholder = "My_Signature")
                                      ),
                                      column(3,
                                             br(),
                                             actionButton("save_sig_btn", " Save Selection (.txt)",
                                                          icon  = icon("save"),
                                                          class = "btn-success btn-block")
                                      ),
                                      column(2, br(),
                                             actionButton("select_all_btn", "Select All",
                                                          icon  = icon("check-square"),
                                                          class = "btn-info btn-block")),
                                      column(2, br(),
                                             actionButton("deselect_all_btn", "Deselect All",
                                                          icon  = icon("square"),
                                                          class = "btn-warning btn-block"))
                                    )
                                  )
                                ),
                                fluidRow(
                                  shinydashboard::box(
                                    title = "DMP Results Table", status = "primary",
                                    solidHeader = TRUE, width = 12,
                                    DT::DTOutput("dmpTable")
                                  )
                                )
        ),

        # ── Heatmap ────────────────────────────────────────────────────────
        shinydashboard::tabItem(tabName = "tab_heatmap",
                                fluidRow(
                                  shinydashboard::box(
                                    title = "Clustering Heatmap (Selected rows from Table)",
                                    status = "primary", solidHeader = TRUE, width = 12,
                                    div(class = "alert alert-info",
                                        style = "font-size:12px; padding:8px;",
                                        icon("info-circle"),
                                        " Select >= 2 rows in the 'Table & Signatures' tab,",
                                        " then switch to this tab."),
                                    plotly::plotlyOutput("heatmapPlot", height = "800px")
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

    observe({
      updateSelectizeInput(session, "sample_filter",
                           choices = colnames(beta), server = TRUE)
    })

    # ── Filtered data ────────────────────────────────────────────────────────
    data_selected <- reactive({
      req(input$selected_comparison)
      df <- dmp_results[[input$selected_comparison]]

      # Normalise column names
      if ("cgi" %in% colnames(df))
        colnames(df)[colnames(df) == "cgi"] <- "CGI"
      if (!"CGI"     %in% colnames(df)) df$CGI     <- "N/A"
      if (!"feature" %in% colnames(df)) df$feature <- "N/A"

      if (!"P.Value" %in% colnames(df))
        df$P.Value <- if ("P.Val" %in% colnames(df)) df$P.Val
      else df$adj.P.Val

      if (!"deltaBeta" %in% colnames(df)) {
        if ("logFC" %in% colnames(df)) df$deltaBeta <- df$logFC
        else stop("[mETHYLotest DMP UI] Missing 'deltaBeta' or 'logFC'.")
      }
      df$logFC <- df$deltaBeta

      # Status (all three thresholds)
      df$Status <- ifelse(
        df$adj.P.Val < input$pval_cutoff &
          df$P.Value  < input$raw_pval_cutoff &
          df$deltaBeta > input$deltabeta_cutoff,
        "Hyper",
        ifelse(
          df$adj.P.Val < input$pval_cutoff &
            df$P.Value  < input$raw_pval_cutoff &
            df$deltaBeta < -input$deltabeta_cutoff,
          "Hypo",
          "Not Sig"
        )
      )

      chr_clean <- gsub("chr", "", as.character(df$CHR), ignore.case = TRUE)
      chr_lvls  <- c(as.character(1:22), "X", "Y", "M", "MT")
      df$CHR_Factor <- factor(chr_clean, levels = chr_lvls)
      df
    })

    data_sig <- reactive({
      df <- data_selected()
      df[df$Status != "Not Sig", ]
    })

    data_for_plots <- reactive({
      df <- data_sig()
      if (nrow(df) == 0L) return(NULL)
      df_all        <- df
      df_all$Status <- "All"
      out           <- rbind(df, df_all)
      out$Status    <- factor(out$Status, levels = c("Hyper", "Hypo", "All"))
      out
    })

    data_table_view <- reactive({
      df   <- data_sig()
      cols <- c("deltaBeta", "adj.P.Val", "P.Value", "CHR", "MAPINFO",
                "gene", "feature", "CGI", "Status")
      df[, intersect(c(cols, grep("_AVG$", colnames(df), value = TRUE)),
                     colnames(df)), drop = FALSE]
    })

    observe({
      df <- data_table_view()
      updateSelectizeInput(session, "target_cpg_selector",
                           choices = rownames(df), server = TRUE)
    })

    # ── Sidebar summary callout ────────────────────────────────────────────
    output$summary_callout <- renderUI({
      df    <- data_selected()
      n_sig <- sum(df$Status != "Not Sig")
      n_h   <- sum(df$Status == "Hyper")
      n_ho  <- sum(df$Status == "Hypo")
      div(class = "callout callout-info",
          style = paste0("padding:8px; border-left-color:#00c0ef;",
                         " background-color:#222d32 !important;",
                         " border:1px solid #444;"),
          p(style = "color:#b8c7ce;", strong("Total DMPs: "),  nrow(df)),
          p(style = "color:#b8c7ce;", strong("Significant: "), n_sig),
          p(style = "color:#e74c3c;", strong("Hyper: "),        n_h),
          p(style = "color:#3498db;", strong("Hypo: "),         n_ho)
      )
    })

    # ── Volcano ───────────────────────────────────────────────────────────
    output$volcanoPlot <- plotly::renderPlotly({
      df     <- data_selected()
      df_sig <- df[df$Status != "Not Sig", ]
      # Downsample non-significant
      df_ns  <- df[df$Status == "Not Sig", ]
      if (nrow(df_ns) > 5000L)
        df_ns <- df_ns[sample(nrow(df_ns), 5000L), ]
      df_viz <- rbind(df_sig, df_ns)

      plotly::plot_ly(
        data         = df_viz,
        x            = ~deltaBeta,
        y            = ~ -log10(adj.P.Val),
        color        = ~Status,
        colors       = c("Hyper" = "#d9534f", "Hypo" = "#0275d8",
                         "Not Sig" = "grey70"),
        text         = ~paste0("CpG: ", rownames(df_viz),
                               "<br>Gene: ",
                               if ("gene" %in% colnames(df_viz))
                                 df_viz$gene else "N/A"),
        type         = "scatter",
        mode         = "markers",
        marker       = list(size = 5, opacity = 0.7),
        hovertemplate = "%{text}<br>dB: %{x:.3f}<br>-log10p: %{y:.2f}<extra></extra>"
      ) %>%
        plotly::layout(
          title  = list(text = "Volcano Plot", x = 0.5),
          xaxis  = list(title = "Delta Beta",
                        zeroline = TRUE,
                        zerolinecolor = "#aaa"),
          yaxis  = list(title = "-log10 Adj. P-value")
        )
    })

    # ── Density ───────────────────────────────────────────────────────────
    output$plot_density <- plotly::renderPlotly({
      df <- data_sig()
      if (nrow(df) == 0L) return(plotly::plotly_empty())

      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = deltaBeta, fill = Status)
      ) +
        ggplot2::geom_density(alpha = 0.6) +
        ggplot2::scale_fill_manual(
          values = c("Hyper" = "#d9534f", "Hypo" = "#0275d8")
        ) +
        ggplot2::labs(x = "Delta Beta", y = "Density") +
        ggplot2::theme_minimal()
      plotly::ggplotly(p)
    })

    # ── ChromoMap ─────────────────────────────────────────────────────────
    output$karyoPlot <- chromoMap::renderChromoMap({
      df <- data_sig()
      if (nrow(df) == 0L || is.null(probe_features)) return(NULL)

      df$CHR_Clean <- paste0("chr",
                             gsub("chr", "", df$CHR, ignore.case = TRUE))
      df_map <- df[df$CHR_Clean %in% chr_ref$chrom, ]
      if (nrow(df_map) == 0L) return(NULL)

      annot <- data.frame(
        Element = rownames(df_map),
        Chrom   = df_map$CHR_Clean,
        Start   = as.integer(df_map$MAPINFO),
        End     = as.integer(df_map$MAPINFO + 1L),
        Data    = as.numeric(df_map$deltaBeta),
        stringsAsFactors = FALSE
      )
      annot <- stats::na.omit(annot)

      chromoMap::chromoMap(
        list(chr_ref),
        list(annot),
        data_based_color_map = TRUE,
        data_type   = "numeric",
        data_colors = list(c("blue", "white", "red")),
        chr_color   = c("#d3d3d3"),
        legend      = TRUE
      )
    })

    # ── Context barplots ──────────────────────────────────────────────────
    make_bar <- function(aes_x) {
      df <- data_for_plots()
      if (is.null(df)) return(plotly::plotly_empty())
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[aes_x]], fill = Status)
      ) +
        ggplot2::geom_bar(position = "dodge") +
        ggplot2::scale_fill_manual(
          values = c("Hyper" = "#d9534f", "Hypo" = "#0275d8",
                     "All" = "#555555")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
      plotly::ggplotly(p)
    }

    output$plot_cgi     <- plotly::renderPlotly({ make_bar("CGI")     })
    output$plot_feature <- plotly::renderPlotly({ make_bar("feature") })

    output$plot_chr <- plotly::renderPlotly({
      df <- data_for_plots()
      if (is.null(df)) return(plotly::plotly_empty())
      df <- df[!is.na(df$CHR_Factor), ]
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = CHR_Factor, fill = Status)
      ) +
        ggplot2::geom_bar(position = "dodge") +
        ggplot2::scale_fill_manual(
          values = c("Hyper" = "#d9534f", "Hypo" = "#0275d8",
                     "All" = "#555555")
        ) +
        ggplot2::labs(x = "Chromosome", y = "Count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
      plotly::ggplotly(p)
    })

    output$plot_feature_cgi <- plotly::renderPlotly({
      df <- data_sig()
      if (nrow(df) == 0L) return(plotly::plotly_empty())
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = feature, fill = CGI)
      ) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::facet_wrap(~Status) +
        ggplot2::scale_fill_viridis_d() +
        ggplot2::labs(y = "Proportion") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
      plotly::ggplotly(p)
    })

    # ── DMP Table ─────────────────────────────────────────────────────────
    output$dmpTable <- DT::renderDT({
      dat       <- data_table_view()
      cols_pval <- intersect(c("adj.P.Val", "P.Value"), colnames(dat))
      cols_beta <- intersect(c("deltaBeta",
                               grep("_AVG$", colnames(dat), value = TRUE)),
                             colnames(dat))

      DT::datatable(
        dat,
        rownames   = TRUE,
        filter     = "top",
        selection  = "multiple",
        extensions = "Buttons",
        options    = list(
          pageLength = 15,
          dom        = "Bfrtip",
          scrollX    = TRUE,
          buttons    = list(
            list(extend = "csv",   text = "CSV"),
            list(extend = "excel", text = "Excel"),
            list(extend = "copy",  text = "Copy")
          )
        ),
        class = "stripe hover compact"
      ) %>%
        DT::formatRound(columns = cols_beta, digits = 4L) %>%
        DT::formatSignif(columns = cols_pval, digits = 3L) %>%
        DT::formatStyle(
          "deltaBeta",
          color      = DT::styleInterval(0, c("#0275d8", "#d9534f")),
          fontWeight = "bold"
        ) %>%
        DT::formatStyle(
          "Status",
          backgroundColor = DT::styleEqual(
            c("Hyper",   "Hypo"),
            c("#fde8e8", "#e8f0fe")
          )
        )
    })

    # ── Heatmap ───────────────────────────────────────────────────────────
    output$heatmapPlot <- plotly::renderPlotly({
      idx <- input$dmpTable_rows_selected
      if (is.null(idx) || length(idx) < 2L)
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "Select >= 2 rows in the Table tab."))

      sel   <- rownames(data_table_view())[idx]
      valid <- intersect(sel, rownames(beta))

      if (length(valid) < 2L)
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "Not enough valid probes."))

      heatmaply::heatmaply(
        beta[valid, , drop = FALSE],
        limits    = c(0, 1),
        colors    = grDevices::colorRampPalette(
          c("navy", "white", "firebrick"))(100L),
        main      = paste0("Heatmap: ", length(valid), " DMPs"),
        xlab      = "Samples",
        ylab      = "CpG Probes"
      )
    })

    # ── Genomic Window ────────────────────────────────────────────────────
    selected_target_info <- reactive({
      tid <- input$target_cpg_selector
      if (is.null(tid) || tid == "" || is.null(probe_features)) return(NULL)
      row <- probe_features[probe_features$CpG == tid, ]
      if (nrow(row) == 0L) return(NULL)
      list(id  = tid,
           chr = as.character(row$CHR_Clean[1L]),
           pos = as.numeric(row$MAPINFO[1L]))
    })

    output$selected_cpg_info <- renderText({
      s <- selected_target_info()
      if (is.null(s)) "Select a CpG from the dropdown."
      else paste0("Target: ", s$id,
                  " | Chr: ", s$chr,
                  " | Pos: ", format(s$pos, big.mark = ","))
    })

    output$windowPlot <- plotly::renderPlotly({
      sel <- selected_target_info()
      if (is.null(sel) || is.null(probe_features))
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "Waiting for selection..."))

      w_start    <- sel$pos - input$window_bp
      w_end      <- sel$pos + input$window_bp
      neighbors  <- probe_features[
        probe_features$CHR_Clean == sel$chr &
          !is.na(probe_features$MAPINFO) &
          probe_features$MAPINFO >= w_start &
          probe_features$MAPINFO <= w_end, ]
      valid_ids  <- intersect(neighbors$CpG, rownames(beta))

      if (length(valid_ids) == 0L)
        return(plotly::plotly_empty() %>%
                 plotly::layout(title = "No probes in this region."))

      mat       <- beta[valid_ids, , drop = FALSE]
      pos_map   <- setNames(
        probe_features$MAPINFO[match(valid_ids, probe_features$CpG)],
        valid_ids
      )

      # Long format without tidyr pivot issues with older R
      long_data <- data.frame(
        CpG    = rep(valid_ids, ncol(mat)),
        Pos    = rep(pos_map,   ncol(mat)),
        Sample = rep(colnames(mat), each = length(valid_ids)),
        Beta   = as.vector(mat),
        stringsAsFactors = FALSE
      )

      # Sample filter
      if (!is.null(input$sample_filter) && length(input$sample_filter) > 0L)
        long_data <- long_data[long_data$Sample %in% input$sample_filter, ]

      # Age mapping (fixed: vectorised, not inline string concat)
      if (!is.null(ages) && length(ages) == ncol(beta)) {
        age_map        <- data.frame(Sample = colnames(beta),
                                     Age    = ages,
                                     stringsAsFactors = FALSE)
        long_data      <- merge(long_data, age_map, by = "Sample",
                                all.x = TRUE)
        if (isTRUE(input$hide_missing_age))
          long_data <- long_data[!is.na(long_data$Age) &
                                   long_data$Age != "", ]
      } else {
        long_data$Age <- NA_character_
      }

      # Pheno mapping
      if (!is.null(pheno) && length(pheno) == ncol(beta)) {
        pheno_map  <- data.frame(Sample = colnames(beta),
                                 Group  = as.character(pheno),
                                 stringsAsFactors = FALSE)
        long_data  <- merge(long_data, pheno_map, by = "Sample", all.x = TRUE)
      } else {
        long_data$Group <- "All"
      }

      # Hover text (vectorised)
      long_data$hover_txt <- paste0(
        "CpG: ",    long_data$CpG,
        "<br>Sample: ", long_data$Sample,
        "<br>Group: ",  long_data$Group,
        if (!all(is.na(long_data$Age)))
          paste0("<br>Age: ", long_data$Age) else "",
        "<br>Pos: ",   format(long_data$Pos, big.mark = ","),
        "<br>Beta: ",  round(long_data$Beta, 3L)
      )

      p <- ggplot2::ggplot(
        long_data,
        ggplot2::aes(x = Pos, y = Beta, color = Group, fill = Group)
      ) +
        ggplot2::geom_point(
          ggplot2::aes(text = hover_txt),
          alpha = 0.6, size = 1.5
        ) +
        ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                             alpha = 0.2, linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = sel$pos,
                            linetype = "dashed", colour = "black") +
        ggplot2::scale_y_continuous(
          labels = scales::percent_format(accuracy = 1),
          limits = c(0, 1.05)
        ) +
        ggplot2::scale_x_continuous(labels = scales::comma) +
        ggplot2::labs(
          title = paste0("Region: chr", sel$chr, "  ",
                         format(w_start, big.mark = ","), " — ",
                         format(w_end,   big.mark = ",")),
          x = "Genomic Position (bp)",
          y = "Methylation Beta"
        ) +
        ggplot2::scale_colour_brewer(palette = "Set1") +
        ggplot2::scale_fill_brewer(palette   = "Set1") +
        ggplot2::theme_minimal(base_size = 12)

      plotly::ggplotly(p, tooltip = "text")
    })

    # ── Table actions ─────────────────────────────────────────────────────
    proxy <- DT::dataTableProxy("dmpTable")
    observeEvent(input$select_all_btn, {
      DT::selectRows(proxy, input$dmpTable_rows_all)
    })
    observeEvent(input$deselect_all_btn, {
      DT::selectRows(proxy, NULL)
    })

    # ── Save signature ────────────────────────────────────────────────────
    observeEvent(input$save_sig_btn, {
      idx  <- input$dmpTable_rows_selected
      name <- trimws(input$sig_name)

      if (is.null(idx) || length(idx) == 0L || nchar(name) == 0L) {
        showNotification("Provide a name and select at least one row.",
                         type = "error")
        return()
      }

      cpgs  <- rownames(data_table_view())[idx]
      fname <- gsub("[^a-zA-Z0-9_\\-]", "_", name)

      if (!dir.exists(path_to_episignature))
        dir.create(path_to_episignature, recursive = TRUE)

      fpath <- file.path(path_to_episignature, paste0(fname, ".txt"))
      writeLines(cpgs, fpath)
      message("[mETHYLotest DMP UI] Signature saved: ", fpath)

      showModal(modalDialog(
        title    = "Signature Saved",
        tagList(
          p(strong("Name: "), name),
          p(strong("CpGs: "), length(cpgs)),
          p(strong("Path: "), tags$code(
            style = "word-break:break-all;", fpath))
        ),
        easyClose = TRUE,
        footer    = modalButton("Close")
      ))
    })

    # ── Close ─────────────────────────────────────────────────────────────
    observeEvent(input$close_btn, {
      shiny::stopApp(returnValue = path_to_episignature)
    })
    session$onSessionEnded(function() {
      shiny::stopApp(returnValue = path_to_episignature)
    })
  }

  message("[mETHYLotest] Starting DMP Explorer UI...")
  invisible(shiny::runApp(shiny::shinyApp(ui, server),
                          host = "0.0.0.0",
                          launch.browser = TRUE))
}
