#' Launch mETHYLotest Project Setup UI (WGBS & Long Read)
#'
#' @description
#' Launches a Shiny Dashboard interface to initialize or reload a mETHYLotest
#' project for WGBS and Long Read methylation analysis.
#'
#' @return Absolute path of the project directory (invisibly).
#' @export
#' @import shinydashboard shinyjs
#' @importFrom shiny shinyApp runApp fluidRow column div br hr h3 h4 h5 p div
#'   strong span tags icon helpText textInput numericInput selectInput
#'   sliderInput checkboxInput checkboxGroupInput conditionalPanel radioButtons
#'   splitLayout actionButton wellPanel uiOutput renderUI tagList
#'   showNotification observeEvent observe reactive req reactiveVal isolate
#'   updateTextInput updateNumericInput updateSelectInput updateCheckboxInput
#'   updateRadioButtons updateSliderInput stopApp
#' @importFrom shinyFiles shinyDirButton shinyFilesButton shinyDirChoose
#'   shinyFileChoose getVolumes parseDirPath parseFilePaths
#' @importFrom fs path_home
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom parallel detectCores
mETHYLotest.NGS.ProjectUI <- function(prefill_pheno = NULL) {

  for (pkg in c("shinyFiles", "writexl", "parallel",
                "shinydashboard", "shinyjs"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", pkg))

  total_cores       <- parallel::detectCores()
  max_cores_allowed <- max(1L, total_cores - 2L)
  default_cores     <- max(1L, floor(max_cores_allowed / 2L))

  TECH_PRESETS <- list(
    modkit = list(
      label = "Nanopore - Modkit (bedMethyl)",
      description = paste(
        "bedMethyl format from Oxford Nanopore modkit.",
        "Col 10 = N_valid_cov, col 11 = percent_modified."),
      col_labels = c("chrom", "start", "end", "name", "score", "strand",
                     "thickStart", "thickEnd", "color", "N_valid_cov",
                     "percent_modified", "N_mod", "N_canonical",
                     "N_other_mod", "N_delete", "N_fail", "N_diff",
                     "N_nocall"),
      chr = 1L, start = 2L, end = 3L, strand = 6L,
      coverage = 10L, freqC = 11L, fraction = FALSE, offset = 1L),
    f5c = list(
      label = "Nanopore - f5c / Nanopolish (methylation_frequency)",
      description = paste(
        "methylation_frequency.tsv from f5c or Nanopolish.",
        "Col 8 = methylated_frequency (fraction 0-1)."),
      col_labels = c("chromosome", "start", "end", "strand",
                     "num_cpgs_in_group", "called_sites",
                     "called_sites_methylated",
                     "methylated_frequency", "group_sequence"),
      chr = 1L, start = 2L, end = 3L, strand = 4L,
      coverage = 6L, freqC = 8L, fraction = TRUE, offset = 1L),
    pbcpg = list(
      label = "PacBio - pb-CpG-tools (combined.bed)",
      description = paste(
        "combined.bed / hap1.bed / hap2.bed from pb-CpG-tools.",
        "Col 10 = coverage, col 11 = percent methylation."),
      col_labels = c("chrom", "chromStart", "chromEnd", "modification",
                     "score", "strand", "thickStart", "thickEnd",
                     "itemRgb", "coverage", "percentMeth"),
      chr = 1L, start = 2L, end = 3L, strand = 6L,
      coverage = 10L, freqC = 11L, fraction = FALSE, offset = 1L)
  )

  fmt_bool <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"
  fmt_vec  <- function(v) {
    if (is.null(v) || length(v) == 0L) return("NULL")
    paste0('c("', paste(v, collapse = '", "'), '")')
  }
  fmt_sep <- function(key) {
    switch(key, tab = '"\\t"', space = '" "',
           comma = '","', semicolon = '";"', '"\\t"')
  }
  sep_to_key <- function(s) {
    switch(s, "\t" = "tab", " " = "space",
           "," = "comma", ";" = "semicolon", "tab")
  }
  null_or <- function(a, b) {
    if (!is.null(a) && length(a) > 0L) a else b
  }

  # =========================================================================
  # UI
  # =========================================================================
  ui <- shinydashboard::dashboardPage(
    skin = "blue",
    shinydashboard::dashboardHeader(
      title = "mETHYLotest NGS Setup"),

    shinydashboard::dashboardSidebar(
      width = 300,
      div(style = "padding:15px;",
          textInput("project_name", "Project Name",
                    placeholder = "My_NGS_Project"),
          strong("Base Directory:"),
          splitLayout(
            cellWidths = c("75%", "25%"),
            textInput("base_dir", label = NULL, value = getwd()),
            shinyFiles::shinyDirButton("btn_base_dir", "Browse",
                                       "Select",
                                       icon = icon("folder-open"),
                                       class = "btn-primary btn-sm")
          ),
          uiOutput("project_status_alert_sidebar")
      ),
      shinydashboard::sidebarMenu(
        id = "tabs",
        shinydashboard::menuItem("1. Data Source",
                                 tabName = "data",
                                 icon = icon("database")),
        shinydashboard::menuItem("2. Import Options",
                                 tabName = "import",
                                 icon = icon("file-import")),
        shinydashboard::menuItem("3. QC & Filtering",
                                 tabName = "qc",
                                 icon = icon("filter")),
        shinydashboard::menuItem("4. Analysis Options",
                                 tabName = "analysis",
                                 icon = icon("cogs")),
        shinydashboard::menuItem("5. Review & Launch",
                                 tabName = "launch",
                                 icon = icon("rocket"))
      )
    ),

    shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),
      tags$head(tags$style(HTML("
        .param-section {
          margin-bottom: 15px; padding: 10px;
          background: #f9f9f9; border-radius: 5px;
          border-left: 3px solid #3c8dbc;
        }
      "))),

      shinydashboard::tabItems(

        # ================================================================
        # Tab 1 : Data Source
        # ================================================================
        shinydashboard::tabItem(
          tabName = "data",
          fluidRow(
            shinydashboard::box(
              title = "Input Data", status = "primary",
              solidHeader = TRUE, width = 12,
              uiOutput("locked_warning_msg"),

              h4(icon("file-excel"), " Phenotype / Metadata File"),
              p("Select the Excel file describing your samples.",
                " Required columns: ",
                tags$code("Sample_ID"), ", ",
                tags$code("Path"), ", ",
                tags$code("Treatment"), "."),

              div(
                class = "alert alert-info",
                style = "font-size:12px; padding:10px;",
                icon("info-circle"),
                " Any additional columns (e.g. ",
                tags$code("Batch"), ", ", tags$code("Sex"),
                ", ", tags$code("Age"),
                ") will be available for batch correction",
                " and covariate adjustment."
              ),

              splitLayout(
                cellWidths = c("85%", "15%"),
                textInput("pheno_path", label = NULL,
                          placeholder = "/path/to/samples.xlsx"),
                shinyFiles::shinyFilesButton(
                  "btn_pheno_path", "Browse", "Select",
                  multiple = FALSE, icon = icon("file-excel"),
                  class = "btn-primary")
              ),
              br(),
              actionButton("load_check_btn",
                           "Load & Validate Metadata",
                           icon = icon("check-circle"),
                           class = "btn-success btn-lg"),
              hr(),
              uiOutput("load_status_msg")
            )
          )
        ),

        # ================================================================
        # Tab 2 : Import Options
        # ================================================================
        shinydashboard::tabItem(
          tabName = "import",
          fluidRow(

            shinydashboard::box(
              title = "methylKit Import Parameters",
              status = "warning", solidHeader = TRUE, width = 6,

              div(class = "param-section",
                  textInput("mk_assembly", "Genome Assembly",
                            value = "hg38"),
                  selectInput("mk_context", "Methylation Context",
                              choices = c("CpG", "CHG", "CHH"),
                              selected = "CpG"),
                  checkboxInput("save_raw_rds",
                                "Save raw imported object (.rds)?",
                                value = TRUE)
              ),

              hr(),
              strong("File Format / Pipeline"),
              br(), br(),

              radioButtons(
                "pipeline_type", label = NULL,
                choices = c(
                  "Standard (Bismark, AMP...)"         = "standard",
                  "Custom column mapping (Nanopore, PacBio...)" = "custom"
                ),
                selected = "standard"
              ),

              conditionalPanel(
                "input.pipeline_type === 'standard'",
                selectInput("mk_pipeline", "Pipeline:",
                            choices = c("bismarkCytosineReport",
                                        "bismarkCoverage",
                                        "bismark", "amp"),
                            selected = "bismarkCytosineReport")
              ),

              conditionalPanel(
                "input.pipeline_type === 'custom'",
                wellPanel(
                  style = "background:#fdfbe6; padding:12px;",

                  selectInput(
                    "tech_preset",
                    label = tagList(icon("magic"),
                                    " Quick preset (optional):"),
                    choices = c(
                      "-- Manual --"                    = "manual",
                      "Nanopore - Modkit (bedMethyl)"   = "modkit",
                      "Nanopore - f5c / Nanopolish"     = "f5c",
                      "PacBio - pb-CpG-tools"           = "pbcpg"
                    ),
                    selected = "manual"
                  ),

                  uiOutput("preset_info_box"),
                  hr(style = "margin:8px 0;"),

                  strong("Column mapping:"),
                  helpText("1-based column numbers. 0 = omit."),

                  fluidRow(
                    column(4, numericInput("col_chr", "Chromosome",
                                           value = 1, min = 1)),
                    column(4, numericInput("col_start", "Start",
                                           value = 2, min = 1)),
                    column(4, numericInput("col_end", "End (0=omit)",
                                           value = 3, min = 0))
                  ),
                  fluidRow(
                    column(4, numericInput("col_coverage", "Coverage",
                                           value = 7, min = 1)),
                    column(4, numericInput("col_freqC", "Meth. freq.",
                                           value = 5, min = 1)),
                    column(4, numericInput("col_strand", "Strand (0=omit)",
                                           value = 6, min = 0))
                  ),

                  checkboxInput(
                    "col_fraction",
                    span(icon("percent"),
                         " Values are fractions (0-1) not percentages"),
                    value = FALSE),

                  uiOutput("pipeline_list_preview"),

                  hr(style = "margin:10px 0;"),
                  strong(icon("arrows-alt-h"), " Coordinate Offset"),
                  br(), br(),

                  div(
                    class = "alert alert-warning",
                    style = "font-size:12px; padding:10px;",
                    icon("exclamation-triangle"),
                    strong(" Coordinate system:"),
                    tags$table(
                      class = "table table-condensed table-bordered",
                      style = paste0("font-size:11px; margin:5px 0;",
                                     " background:#fff; color:#333;"),
                      tags$thead(tags$tr(
                        tags$th("Tool"), tags$th("System"),
                        tags$th("Offset"))),
                      tags$tbody(
                        tags$tr(tags$td("modkit"),
                                tags$td("0-based BED"),
                                tags$td(strong("+1",
                                               style = "color:#c0392b;"))),
                        tags$tr(tags$td("f5c / Nanopolish"),
                                tags$td("0-based BED"),
                                tags$td(strong("+1",
                                               style = "color:#c0392b;"))),
                        tags$tr(tags$td("pb-CpG-tools"),
                                tags$td("0-based BED"),
                                tags$td(strong("+1",
                                               style = "color:#c0392b;"))),
                        tags$tr(tags$td("Bismark"),
                                tags$td("1-based"),
                                tags$td(strong("0",
                                               style = "color:#27ae60;")))
                      )
                    )
                  ),

                  numericInput("coord_offset",
                               "Offset value (post-import):",
                               value = 0L, min = -100L,
                               max = 100L, step = 1L)
                )
              )
            ),

            shinydashboard::box(
              title = "Advanced Import (methRead)",
              status = "danger", solidHeader = TRUE, width = 6,
              collapsible = TRUE, collapsed = TRUE,

              div(class = "alert alert-warning",
                  style = "margin-bottom:12px;",
                  icon("exclamation-triangle"),
                  strong(" Low-level parameters."),
                  " Defaults cover most cases."),

              fluidRow(
                column(4,
                       selectInput("mk_dbtype", "DB backend (dbtype)",
                                   choices = c("None (RAM)" = "none",
                                               "Tabix"      = "tabix"),
                                   selected = "none")
                ),
                column(4,
                       checkboxInput("mk_header", "File has header",
                                     value = FALSE),
                       numericInput("mk_skip", "Lines to skip",
                                    value = 0, min = 0)
                ),
                column(4,
                       selectInput("mk_sep", "Separator",
                                   choices = c("Tab" = "tab",
                                               "Space" = "space",
                                               "Comma" = "comma",
                                               "Semicolon" = "semicolon"),
                                   selected = "tab"),
                       selectInput("mk_resolution", "Resolution",
                                   choices = c("Base" = "base",
                                               "Region" = "region"),
                                   selected = "base")
                )
              )
            )
          ),

          fluidRow(
            shinydashboard::box(
              title = "System Resources",
              status = "info", solidHeader = TRUE, width = 12,
              fluidRow(
                column(6,
                       sliderInput(
                         "num_cores",
                         paste0("Processing Cores (Max: ",
                                max_cores_allowed, ")"),
                         min = 1, max = max_cores_allowed,
                         value = default_cores, step = 1),
                       helpText(paste("Detected:", total_cores))
                )
              )
            )
          )
        ),

        # ================================================================
        # Tab 3 : QC & Filtering
        # ================================================================
        shinydashboard::tabItem(
          tabName = "qc",
          fluidRow(
            shinydashboard::box(
              title = "Initial QC Parameters",
              status = "warning", solidHeader = TRUE, width = 12,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px;",
                  icon("info-circle"),
                  "These are the ", strong("starting values"),
                  " for the iterative QC loop.",
                  " You can refine them interactively during the",
                  " pipeline execution."),

              fluidRow(
                column(4,
                       numericInput(
                         "qc_hi_perc",
                         "High coverage percentile (hi.perc)",
                         value = 99.9, min = 90, max = 100,
                         step = 0.1),
                       helpText("Bases with coverage above this",
                                " percentile are removed",
                                " (PCR duplicates).")
                ),
                column(4,
                       numericInput(
                         "qc_lo_count",
                         "Minimum coverage (lo.count)",
                         value = 10, min = 0, max = 100,
                         step = 1),
                       helpText("Minimum read count per base.",
                                " Used for both import and QC filtering.")
                ),
                column(4,
                       checkboxInput(
                         "qc_use_lo_perc",
                         "Use low coverage percentile instead",
                         value = FALSE),
                       conditionalPanel(
                         "input.qc_use_lo_perc == true",
                         numericInput(
                           "qc_lo_perc",
                           "Low coverage percentile (lo.perc)",
                           value = NULL, min = 0, max = 50,
                           step = 0.1)
                       )
                )
              )
            )
          )
        ),

        # ================================================================
        # Tab 4 : Analysis Options
        # ================================================================
        shinydashboard::tabItem(
          tabName = "analysis",

          # Unite
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("object-group"),
                              " Sample Unification (unite)"),
              status = "primary", solidHeader = TRUE, width = 12,

              div(class = "param-section",
                  checkboxInput("unite_destrand",
                                "Merge both strands (destrand)",
                                value = FALSE),
                  helpText("If TRUE, reads from both strands of a",
                           " CpG are merged.",
                           " Usually TRUE for CpG context, FALSE for",
                           " CHG/CHH.")
              )
            )
          ),

          # Coverage Normalization
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("balance-scale"),
                              " Coverage Normalization"),
              status = "primary", solidHeader = TRUE, width = 12,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px;",
                  icon("info-circle"),
                  " Normalizes coverage between samples using ",
                  tags$code("methylKit::normalizeCoverage()"), ".",
                  " Recommended when samples have very different",
                  " sequencing depths."),

              checkboxInput("do_normalize_coverage",
                            strong(" Enable Coverage Normalization"),
                            value = FALSE),

              conditionalPanel(
                "input.do_normalize_coverage == true",
                div(class = "param-section",
                    selectInput("normalize_cov_method",
                                "Method",
                                choices = c("median", "mean"),
                                selected = "median"),
                    helpText(
                      tags$b("median:"),
                      " Scales to median coverage (robust).", br(),
                      tags$b("mean:"),
                      " Scales to mean coverage.")
                )
              )
            )
          ),

          # Clustering
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("sitemap"),
                              " Clustering & Correlation"),
              status = "info", solidHeader = TRUE, width = 12,

              div(class = "param-section",
                  fluidRow(
                    column(6,
                           selectInput("cluster_dist",
                                       "Distance metric (dist)",
                                       choices = c("correlation",
                                                   "euclidean",
                                                   "maximum",
                                                   "manhattan"),
                                       selected = "correlation")
                    ),
                    column(6,
                           selectInput("cluster_method",
                                       "Clustering method",
                                       choices = c("ward",
                                                   "complete",
                                                   "average",
                                                   "single"),
                                       selected = "ward")
                    )
                  )
              )
            )
          ),

          # Tiling Windows
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("th"), " Tiling Windows (DMR detection)"),
              status = "primary", solidHeader = TRUE, width = 12,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px;",
                  icon("info-circle"),
                  " Aggregates per-base methylation into fixed-size",
                  " genomic windows for region-level differential",
                  " analysis. Recommended for WGBS data."),

              checkboxInput("do_tiling", strong(" Enable Tiling Windows"),
                            value = FALSE),

              conditionalPanel(
                "input.do_tiling == true",
                div(class = "param-section",
                    fluidRow(
                      column(4,
                             numericInput("tiling_win_size",
                                          "Window size (bp)",
                                          value = 1000, min = 100,
                                          max = 10000, step = 100),
                             helpText("Default 1000bp. Larger = smoother.")
                      ),
                      column(4,
                             numericInput("tiling_step_size",
                                          "Step size (bp)",
                                          value = 1000, min = 100,
                                          max = 10000, step = 100),
                             helpText("Non-overlapping if step = window.")
                      ),
                      column(4,
                             numericInput("tiling_min_cov",
                                          "Min CpGs per window",
                                          value = 3, min = 1,
                                          max = 20, step = 1),
                             helpText("Windows with fewer CpGs are dropped.")
                      )
                    )
                )
              )
            )
          ),

          # Segmentation
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("puzzle-piece"),
                              " Methylation Segmentation (methSeg)"),
              status = "warning", solidHeader = TRUE, width = 12,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px;",
                  icon("info-circle"),
                  " Segments the genome into regions of similar",
                  " methylation levels using a HMM approach.",
                  " Useful for identifying UMRs, LMRs, PMDs, etc."),

              checkboxInput("do_segmentation",
                            strong(" Enable Segmentation"),
                            value = FALSE),

              conditionalPanel(
                "input.do_segmentation == true",
                div(class = "param-section",
                    fluidRow(
                      column(6,
                             numericInput("seg_min_seg",
                                          "Min segment length (bp)",
                                          value = 500, min = 100,
                                          max = 10000, step = 100)
                      ),
                      column(6,
                             numericInput("seg_G",
                                          "Number of HMM states (G)",
                                          value = 4, min = 2,
                                          max = 8, step = 1),
                             helpText("2=hypo/hyper, 3=add intermediate,",
                                      " 4=UMR/LMR/PMD/HMR")
                      )
                    )
                )
              )
            )
          ),

          # Differential Methylation
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("not-equal"),
                              " Differential Methylation (calculateDiffMeth)"),
              status = "success", solidHeader = TRUE, width = 12,

              div(class = "param-section",
                  fluidRow(
                    column(4,
                           selectInput("diff_overdispersion",
                                       "Overdispersion model",
                                       choices = c("MN", "shrinkMN",
                                                   "none"),
                                       selected = "MN"),
                           helpText(
                             tags$b("MN:"),
                             " Beta-binomial (recommended).", br(),
                             tags$b("shrinkMN:"),
                             " Shrunken estimate.", br(),
                             tags$b("none:"),
                             " Binomial (no overdispersion).")
                    ),
                    column(4,
                           selectInput("diff_test",
                                       "Statistical test",
                                       choices = c("Chisq", "F",
                                                   "midPval", "fast.fisher"),
                                       selected = "Chisq"),
                           helpText("Chi-squared is the default for",
                                    " logistic regression.")
                    ),
                    column(4,
                           numericInput("diff_cutoff",
                                        "Min methylation difference (%)",
                                        value = 10, min = 0,
                                        max = 100, step = 1),
                           helpText("Minimum absolute methylation",
                                    " difference to call a DMC.")
                    )
                  ),
                  fluidRow(
                    column(4,
                           numericInput("diff_qvalue",
                                        "Q-value cutoff",
                                        value = 0.05, min = 0,
                                        max = 1, step = 0.01),
                           helpText("FDR-adjusted p-value threshold.")
                    )
                  )
              )
            )
          ),

          # Annotation
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("tag"),
                              " Genomic Annotation (annotatr)"),
              status = "warning", solidHeader = TRUE, width = 12,

              div(class = "param-section",
                  fluidRow(
                    column(6,
                           numericInput("annot_qval_cutoff",
                                        "Annotation Q-value cutoff",
                                        value = 0.05, min = 0,
                                        max = 1, step = 0.01)
                    )
                  )
              )
            )
          )
        ),

        # ================================================================
        # Tab 5 : Review & Launch
        # ================================================================
        shinydashboard::tabItem(
          tabName = "launch",
          uiOutput("launch_ui_dynamic")
        )
      )
    )
  )

  # =========================================================================
  # Server
  # =========================================================================
  server <- function(input, output, session) {

    loaded_data         <- reactiveVal(NULL)
    pheno_columns       <- reactiveVal(NULL)
    final_path          <- reactiveVal(NULL)
    is_existing_project <- reactiveVal(FALSE)

    os_volumes <- shinyFiles::getVolumes()()
    my_roots   <- c(Home = fs::path_home(), os_volumes)

    LOCK_IDS <- c(
      "project_name", "base_dir",
      "pheno_path", "btn_pheno_path", "load_check_btn",
      "mk_assembly", "mk_context", "save_raw_rds",
      "pipeline_type", "mk_pipeline", "tech_preset",
      "do_normalize_coverage", "normalize_cov_method",
      "col_chr", "col_start", "col_end",
      "col_coverage", "col_freqC", "col_strand",
      "col_fraction", "coord_offset",
      "mk_dbtype", "mk_header", "mk_skip", "mk_sep", "mk_resolution",
      "num_cores",
      "qc_hi_perc", "qc_lo_count", "qc_use_lo_perc", "qc_lo_perc",
      "unite_destrand",
      "cluster_dist", "cluster_method",
      "diff_overdispersion", "diff_test", "diff_qvalue",
      "annot_diff_cutoff", "annot_qval_cutoff"
    )

    # ══════════════════════════════════════════════════════════════
    # Auto-fill pheno path (demo mode) — user still clicks Load
    # ══════════════════════════════════════════════════════════════
    observe({
      if (!is.null(prefill_pheno) && file.exists(prefill_pheno)) {
        updateTextInput(session, "pheno_path", value = prefill_pheno)
        showNotification("Demo metadata path set. Click 'Load & Validate Metadata' to continue.",
                         type = "message", duration = 8)
      }
    })

    # -- Preset auto-fill --
    observeEvent(input$tech_preset, {
      req(input$tech_preset != "manual")
      p <- TECH_PRESETS[[input$tech_preset]]
      if (is.null(p)) return()
      updateNumericInput(session, "col_chr",      value = p$chr)
      updateNumericInput(session, "col_start",    value = p$start)
      updateNumericInput(session, "col_end",      value = p$end)
      updateNumericInput(session, "col_coverage", value = p$coverage)
      updateNumericInput(session, "col_freqC",    value = p$freqC)
      updateNumericInput(session, "col_strand",   value = p$strand)
      updateCheckboxInput(session, "col_fraction", value = p$fraction)
      updateNumericInput(session, "coord_offset", value = p$offset)
    })

    # -- Header auto-toggle for custom --
    observeEvent(input$pipeline_type, {
      updateCheckboxInput(session, "mk_header",
                          value = (input$pipeline_type == "custom"))
    })

    # -- Preset info box --
    output$preset_info_box <- renderUI({
      req(input$tech_preset)
      if (input$tech_preset == "manual") return(NULL)
      p <- TECH_PRESETS[[input$tech_preset]]
      if (is.null(p)) return(NULL)

      get_role <- function(i) {
        if (i == p$chr) return("Chromosome")
        if (i == p$start) return("Start")
        if (p$end > 0L && i == p$end) return("End")
        if (p$strand > 0L && i == p$strand) return("Strand")
        if (i == p$coverage) return("Coverage")
        if (i == p$freqC) return("Meth. Freq.")
        ""
      }

      rows <- lapply(seq_along(p$col_labels), function(i) {
        role <- get_role(i)
        st <- if (nchar(role) > 0)
          "background:#fff3cd; color:#000;" else "color:#000;"
        tags$tr(style = st,
                tags$td(i), tags$td(p$col_labels[i]),
                tags$td(strong(role)))
      })

      div(class = "alert alert-info",
          style = "font-size:12px; padding:8px; margin-top:6px;",
          icon("info-circle"), strong(" ", p$label), br(),
          p$description,
          hr(style = "margin:6px 0;"),
          tags$table(
            class = "table table-condensed table-bordered",
            style = "font-size:11px; background:#fff; color:#000;",
            tags$thead(tags$tr(tags$th("#"), tags$th("Column"),
                               tags$th("Used as"))),
            tags$tbody(rows)
          )
      )
    })

    # -- Pipeline list preview --
    output$pipeline_list_preview <- renderUI({
      req(input$pipeline_type == "custom")
      req(input$col_chr, input$col_start,
          input$col_coverage, input$col_freqC)

      end_p <- if (isTRUE(input$col_end > 0L))
        sprintf(",\n  end.col      = %dL", as.integer(input$col_end))
      else ""
      str_p <- if (isTRUE(input$col_strand > 0L))
        sprintf(",\n  strand.col   = %dL", as.integer(input$col_strand))
      else ""

      code <- sprintf(
        paste0("list(\n  fraction     = %s,\n",
               "  chr.col      = %dL,\n",
               "  start.col    = %dL%s%s,\n",
               "  coverage.col = %dL,\n",
               "  freqC.col    = %dL\n)"),
        fmt_bool(input$col_fraction),
        as.integer(input$col_chr), as.integer(input$col_start),
        end_p, str_p,
        as.integer(input$col_coverage), as.integer(input$col_freqC))

      tagList(
        hr(style = "margin:8px 0;"),
        strong(icon("code"), " Generated pipeline list:"),
        tags$pre(style = paste0("background:#f0f0f0; padding:8px;",
                                " font-size:11px; border-radius:4px;"),
                 code))
    })

    # -- Project detection --
    output$project_status_alert_sidebar <- renderUI({
      if (is_existing_project())
        div(class = "callout callout-success",
            style = "padding:10px; margin-top:10px;",
            h5(icon("check-circle"), " Project Found"),
            tags$small("Locked."))
      else
        div(class = "callout callout-info",
            style = "padding:10px; margin-top:10px;",
            h5(icon("plus"), " New Project"))
    })

    output$locked_warning_msg <- renderUI({
      if (is_existing_project())
        div(class = "alert alert-warning",
            icon("lock"), strong(" READ-ONLY"),
            " — existing project detected.")
    })

    observe({
      req(input$base_dir)
      path <- trimws(input$base_dir)
      if (!nchar(path)) return()

      detected_path <- NULL
      if (file.exists(path) && !dir.exists(path) &&
          tools::file_ext(path) == "R" &&
          basename(path) == "project_config.R")
        detected_path <- dirname(dirname(path))
      else if (dir.exists(path) &&
               file.exists(file.path(path, "Results",
                                     "project_config.R")))
        detected_path <- path

      if (!is.null(detected_path)) {
        if (!is_existing_project()) {
          is_existing_project(TRUE)
          updateTextInput(session, "project_name",
                          value = basename(detected_path))
          updateTextInput(session, "base_dir",
                          value = dirname(detected_path))

          cfg_path <- file.path(detected_path, "Results",
                                "project_config.R")
          tryCatch({
            env <- new.env()
            source(cfg_path, local = env)
            cfg <- env$project_config

            .ut <- function(id, k = id)
              if (!is.null(cfg[[k]]))
                updateTextInput(session, id, value = cfg[[k]])
            .un <- function(id, k = id)
              if (!is.null(cfg[[k]]))
                updateNumericInput(session, id, value = cfg[[k]])
            .us <- function(id, k = id)
              if (!is.null(cfg[[k]]))
                updateSelectInput(session, id, selected = cfg[[k]])
            .uc <- function(id, k = id)
              if (!is.null(cfg[[k]]))
                updateCheckboxInput(session, id, value = cfg[[k]])

            .ut("mk_assembly", "assembly")
            .us("mk_context", "context")
            if (!is.null(cfg$min_coverage))
              updateNumericInput(session, "qc_lo_count",
                                 value = cfg$min_coverage)
            .ut("pheno_path", "pheno_file")
            .uc("save_raw_rds", "save_raw_obj")

            if (is.list(cfg$pipeline)) {
              updateRadioButtons(session, "pipeline_type",
                                 selected = "custom")
              p <- cfg$pipeline
              .un("col_chr"); .un("col_start")
              if (!is.null(p$end.col))
                updateNumericInput(session, "col_end",
                                   value = p$end.col)
              .un("col_coverage"); .un("col_freqC")
              if (!is.null(p$strand.col))
                updateNumericInput(session, "col_strand",
                                   value = p$strand.col)
              .uc("col_fraction", "fraction")
              .un("coord_offset")
            } else {
              updateRadioButtons(session, "pipeline_type",
                                 selected = "standard")
              .us("mk_pipeline", "pipeline")
            }

            .us("mk_dbtype")
            .uc("mk_header", "header")
            .un("mk_skip", "skip")
            if (!is.null(cfg$sep))
              updateSelectInput(session, "mk_sep",
                                selected = sep_to_key(cfg$sep))
            .us("mk_resolution", "resolution")
            .un("num_cores")

            # QC params
            .un("qc_hi_perc"); .un("qc_lo_count")
            if (!is.null(cfg$qc_lo_perc)) {
              updateCheckboxInput(session, "qc_use_lo_perc",
                                  value = TRUE)
              updateNumericInput(session, "qc_lo_perc",
                                 value = cfg$qc_lo_perc)
            }

            # Analysis params
            .uc("unite_destrand")
            .uc("do_normalize_coverage")
            .us("normalize_cov_method")
            .us("cluster_dist"); .us("cluster_method")
            .us("diff_overdispersion"); .us("diff_test")
            .un("diff_cutoff"); .un("diff_qvalue")
            .un("annot_qval_cutoff")

            for (id in LOCK_IDS) shinyjs::disable(id)

            if (file.exists(cfg$pheno_file))
              try(loaded_data(readxl::read_excel(cfg$pheno_file)),
                  silent = TRUE)

            showNotification("Existing project loaded.",
                             type = "message", duration = 5)
            shinydashboard::updateTabItems(session, "tabs",
                                           "launch")
          }, error = function(e)
            showNotification(paste("Error:", e$message),
                             type = "error"))
        }
      } else {
        cur <- file.path(input$base_dir, input$project_name,
                         "Results", "project_config.R")
        if (is_existing_project() && !file.exists(cur)) {
          is_existing_project(FALSE)
          for (id in LOCK_IDS) shinyjs::enable(id)
          loaded_data(NULL); pheno_columns(NULL)
        }
      }
    })

    # -- File choosers --
    shinyFiles::shinyDirChoose(input, "btn_base_dir",
                               roots = my_roots, session = session)
    observeEvent(input$btn_base_dir, {
      p <- shinyFiles::parseDirPath(my_roots, input$btn_base_dir)
      if (length(p) > 0)
        updateTextInput(session, "base_dir", value = as.character(p))
    })

    shinyFiles::shinyFileChoose(input, "btn_pheno_path",
                                roots = my_roots, session = session,
                                filetypes = c("xlsx", "xls", "csv"))
    observeEvent(input$btn_pheno_path, {
      f <- shinyFiles::parseFilePaths(my_roots, input$btn_pheno_path)
      if (nrow(f) > 0)
        updateTextInput(session, "pheno_path",
                        value = as.character(f$datapath))
    })

    session$onSessionEnded(function() {
      stopApp(returnValue = isolate(final_path()))
    })

    # -- Data loading --
    observeEvent(input$load_check_btn, {
      req(input$pheno_path)
      if (!file.exists(input$pheno_path)) {
        showNotification("File not found!", type = "error")
        return()
      }
      tryCatch({
        df <- readxl::read_excel(input$pheno_path)
        missing <- setdiff(c("Sample_ID", "Path", "Treatment"),
                           colnames(df))
        if (length(missing) > 0) {
          showNotification(paste("Missing:", paste(missing,
                                                   collapse = ", ")),
                           type = "error")
          loaded_data(NULL)
        } else {
          loaded_data(df)
          pheno_columns(setdiff(colnames(df),
                                c("Sample_ID", "Path", "Treatment")))
          output$load_status_msg <- renderUI(
            div(class = "alert alert-success", icon("check"),
                paste(" Loaded:", nrow(df), "samples.")))
          shinydashboard::updateTabItems(session, "tabs", "import")
        }
      }, error = function(e)
        showNotification(paste("Error:", e$message), type = "error"))
    })

    # -- Launch tab --
    output$launch_ui_dynamic <- renderUI({

      if (is_existing_project())
        return(tagList(
          div(class = "alert alert-info",
              icon("info-circle"), " Project locked."),
          fluidRow(shinydashboard::box(
            title = "Project Status", status = "success",
            solidHeader = TRUE, width = 12,
            div(style = "text-align:center; padding:10px;",
                h3(input$project_name),
                hr(),
                actionButton("relaunch_btn", "Re-launch",
                             icon = icon("rocket"),
                             class = "btn-primary btn-lg",
                             style = "width:50%;"))
          )),
          if (!is.null(loaded_data()))
            fluidRow(shinydashboard::box(
              title = "Samples", status = "primary",
              solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              DT::DTOutput("selection_table")))
        ))

      df <- loaded_data()
      if (is.null(df))
        return(fluidRow(shinydashboard::box(
          width = 12, status = "danger",
          h3(icon("arrow-left"), " No Data"),
          p("Load metadata first."))))

      cols_avail <- pheno_columns()

      tagList(
        fluidRow(shinydashboard::box(
          title = paste("Samples (", nrow(df), ")"),
          status = "primary", solidHeader = TRUE, width = 12,
          DT::DTOutput("selection_table"))),

        fluidRow(shinydashboard::box(
          title = "Statistical Adjustments",
          status = "info", solidHeader = TRUE, width = 12,
          if (length(cols_avail) > 0) tagList(
            fluidRow(
              column(6,
                     checkboxInput("do_batch_correction",
                                   span(icon("layer-group"),
                                        strong(" Batch Correction")),
                                   value = FALSE),
                     conditionalPanel(
                       "input.do_batch_correction == true",
                       checkboxGroupInput("batch_cols",
                                          "Batch variable(s):",
                                          choices = cols_avail,
                                          inline = TRUE))
              ),
              column(6,
                     checkboxInput("do_covariates",
                                   span(icon("project-diagram"),
                                        strong(" Covariates")),
                                   value = FALSE),
                     conditionalPanel(
                       "input.do_covariates == true",
                       checkboxGroupInput("covariate_cols",
                                          "Covariate(s):",
                                          choices = cols_avail,
                                          inline = TRUE))
              )
            )
          ) else
            div(class = "text-warning",
                icon("exclamation-circle"),
                " No extra columns.")
        )),

        fluidRow(shinydashboard::box(
          title = "Summary", status = "info",
          solidHeader = TRUE, width = 12, collapsible = TRUE,
          fluidRow(
            column(3, h5(strong("Import")), tags$ul(
              tags$li(paste("Assembly:", input$mk_assembly)),
              tags$li(paste("Context:", input$mk_context)),
              tags$li(paste("Min cov:", input$qc_lo_count, "(from QC)"))
            )),
            column(3, h5(strong("QC")), tags$ul(
              tags$li(paste("Hi perc:", input$qc_hi_perc)),
              tags$li(paste("Lo count:", input$qc_lo_count))
            )),
            column(3, h5(strong("DiffMeth")), tags$ul(
              tags$li(paste("Test:", input$diff_test)),
              tags$li(paste("Diff:", input$diff_cutoff, "%")),
              tags$li(paste("Q-val:", input$diff_qvalue))
            )),
            column(3, h5(strong("System")), tags$ul(
              tags$li(paste("Cores:", input$num_cores))
            ))
          )
        )),

        fluidRow(column(12, align = "center", hr(),
                        actionButton("save_generate_btn",
                                     "Generate Project",
                                     icon = icon("save"),
                                     class = "btn-success btn-lg",
                                     style = "width:50%;")))
      )
    })

    output$selection_table <- DT::renderDT({
      req(loaded_data())
      DT::datatable(
        loaded_data(),
        selection = list(
          mode = if (is_existing_project()) "none" else "multiple",
          selected = seq_len(nrow(loaded_data())),
          target = "row"),
        options = list(dom = "tp", pageLength = 10, scrollX = TRUE))
    })

    observeEvent(input$relaunch_btn, {
      d <- file.path(input$base_dir, input$project_name)
      final_path(d); stopApp(d)
    })

    # -- Generate config --
    observeEvent(input$save_generate_btn, {
      req(input$project_name, input$base_dir, loaded_data())

      sel <- input$selection_table_rows_selected
      if (length(sel) == 0L) {
        showNotification("Select at least one sample.",
                         type = "error"); return()
      }

      final_df    <- loaded_data()[sel, , drop = FALSE]
      project_dir <- file.path(input$base_dir, input$project_name)
      res_dir     <- file.path(project_dir, "Results")
      data_dir    <- file.path(project_dir, "data")
      rds_dir     <- file.path(project_dir, "data", "interim")
      config_file <- file.path(res_dir, "project_config.R")

      if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
      if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)

      pheno_out <- file.path(data_dir, "selected_samples.xlsx")
      if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
      writexl::write_xlsx(final_df, pheno_out)

      pipeline_str <- if (input$pipeline_type == "standard") {
        sprintf('"%s"', input$mk_pipeline)
      } else {
        end_p <- if (isTRUE(input$col_end > 0L))
          sprintf(",\n  end.col      = %dL",
                  as.integer(input$col_end)) else ""
        str_p <- if (isTRUE(input$col_strand > 0L))
          sprintf(",\n  strand.col   = %dL",
                  as.integer(input$col_strand)) else ""
        sprintf(
          paste0("list(\n  fraction     = %s,\n",
                 "  chr.col      = %dL,\n",
                 "  start.col    = %dL%s%s,\n",
                 "  coverage.col = %dL,\n",
                 "  freqC.col    = %dL\n)"),
          fmt_bool(input$col_fraction),
          as.integer(input$col_chr),
          as.integer(input$col_start), end_p, str_p,
          as.integer(input$col_coverage),
          as.integer(input$col_freqC))
      }

      offset_val <- if (input$pipeline_type == "custom")
        as.integer(input$coord_offset) else 0L
      db_str <- if (input$mk_dbtype == "none") "NA"
      else sprintf('"%s"', input$mk_dbtype)
      batch_vec <- if (isTRUE(input$do_batch_correction))
        input$batch_cols else NULL
      cov_vec <- if (isTRUE(input$do_covariates))
        input$covariate_cols else NULL
      lo_perc_str <- if (isTRUE(input$qc_use_lo_perc) &&
                         !is.null(input$qc_lo_perc))
        as.character(input$qc_lo_perc) else "NULL"

      lines <- c(
        "# mETHYLotest NGS Project Configuration",
        sprintf("# Generated: %s", Sys.time()),
        sprintf("# Samples:   %d", nrow(final_df)),
        "",
        "project_config <- list()",
        "",
        "# Paths",
        sprintf('project_config$project_name <- "%s"',
                input$project_name),
        sprintf('project_config$project_dir  <- "%s"', project_dir),
        sprintf('project_config$res_dir      <- "%s"', res_dir),
        sprintf('project_config$rds_dir      <- "%s"', rds_dir),
        sprintf('project_config$pheno_file   <- "%s"', pheno_out),
        "",
        "# Column mapping",
        'project_config$col_sampleID  <- "Sample_ID"',
        'project_config$col_file_path <- "Path"',
        'project_config$col_treatment <- "Treatment"',
        "",
        "# Import (methylKit::methRead)",
        sprintf('project_config$assembly     <- "%s"',
                input$mk_assembly),
        sprintf('project_config$context      <- "%s"',
                input$mk_context),
        sprintf("project_config$min_coverage <- %s", input$qc_lo_count),
        sprintf("project_config$pipeline     <- %s", pipeline_str),
        sprintf("project_config$coord_offset <- %dL", offset_val),
        sprintf("project_config$dbtype       <- %s", db_str),
        sprintf("project_config$header       <- %s",
                fmt_bool(input$mk_header)),
        sprintf("project_config$skip         <- %dL",
                as.integer(input$mk_skip)),
        sprintf("project_config$sep          <- %s",
                fmt_sep(input$mk_sep)),
        sprintf('project_config$resolution   <- "%s"',
                input$mk_resolution),
        "",
        "# QC initial parameters",
        sprintf("project_config$qc_hi_perc  <- %s",
                input$qc_hi_perc),
        sprintf("project_config$qc_lo_count <- %s",
                input$qc_lo_count),
        sprintf("project_config$qc_lo_perc  <- %s", lo_perc_str),
        "",
        "# Unite",
        sprintf("project_config$unite_destrand <- %s",
                fmt_bool(input$unite_destrand)),
        "",
        "# Coverage Normalization",
        sprintf("project_config$do_normalize_coverage <- %s",
                fmt_bool(input$do_normalize_coverage)),
        sprintf('project_config$normalize_cov_method  <- "%s"',
                null_or(input$normalize_cov_method, "median")),
        "",
        "# Clustering",
        sprintf('project_config$cluster_dist   <- "%s"',
                input$cluster_dist),
        sprintf('project_config$cluster_method <- "%s"',
                input$cluster_method),
        "",
        "# Tiling Windows",
        sprintf("project_config$do_tiling        <- %s",
                fmt_bool(input$do_tiling)),
        sprintf("project_config$tiling_win_size  <- %s",
                null_or(input$tiling_win_size, 1000)),
        sprintf("project_config$tiling_step_size <- %s",
                null_or(input$tiling_step_size, 1000)),
        sprintf("project_config$tiling_min_cov   <- %s",
                null_or(input$tiling_min_cov, 3)),
        "",
        "# Segmentation",
        sprintf("project_config$do_segmentation <- %s",
                fmt_bool(input$do_segmentation)),
        sprintf("project_config$seg_min_seg     <- %s",
                null_or(input$seg_min_seg, 500)),
        sprintf("project_config$seg_G           <- %s",
                null_or(input$seg_G, 4)),
        "# Differential Methylation (calculateDiffMeth)",
        sprintf('project_config$diff_overdispersion <- "%s"',
                input$diff_overdispersion),
        sprintf('project_config$diff_test           <- "%s"',
                input$diff_test),
        sprintf("project_config$diff_cutoff         <- %s",
                input$diff_cutoff),
        sprintf("project_config$diff_qvalue         <- %s",
                input$diff_qvalue),
        "",
        "# Annotation",
        sprintf("project_config$annot_diff_cutoff <- %s",
                input$diff_cutoff),
        sprintf("project_config$annot_qval_cutoff <- %s",
                input$annot_qval_cutoff),
        "",
        "# Statistical adjustments",
        sprintf("project_config$perform_batch_correction <- %s",
                fmt_bool(input$do_batch_correction)),
        sprintf("project_config$batch_cols               <- %s",
                fmt_vec(batch_vec)),
        sprintf("project_config$perform_covariate_adj    <- %s",
                fmt_bool(input$do_covariates)),
        sprintf("project_config$covariates               <- %s",
                fmt_vec(cov_vec)),
        "",
        "# System",
        sprintf("project_config$num_cores    <- %s", input$num_cores),
        sprintf("project_config$save_raw_obj <- %s",
                fmt_bool(input$save_raw_rds))
      )

      tryCatch({
        writeLines(lines, config_file)
        final_path(project_dir)
        showNotification("Project created.", type = "message")
      }, error = function(e)
        showNotification(paste("Failed:", e$message),
                         type = "error"))
    })
  }

  runApp(shinyApp(ui, server), host = "0.0.0.0", launch.browser = TRUE)
}
