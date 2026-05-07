#' Launch mETHYLotest Project Setup UI (EPIC Array)
#'
#' @description
#' Launches a Shiny Dashboard interface to initialize or reload a mETHYLotest
#' EPIC project for array-based methylation analysis (ChAMP pipeline).
#'
#' @details
#' Pipeline order:
#' \enumerate{
#'   \item \code{champ.load} — Import IDAT files
#'   \item \code{champ.filter} — Quality control filtering
#'   \item \code{champ.impute} — Imputation of missing values
#'   \item \code{champ.norm} — Normalisation
#'   \item \code{champ.refbase} — Cell type correction (blood only)
#'   \item \code{champ.runCombat} — Batch correction
#'   \item \code{champ.DMP} — Differentially methylated positions
#'   \item \code{champ.DMR} — Differentially methylated regions
#'   \item \code{champ.Block} — Block methylation analysis
#'   \item \code{champ.GSEA} — Gene set enrichment analysis
#'   \item \code{champ.CNA} — Copy number aberration analysis
#' }
#'
#' Array type (technology) is detected automatically from \code{Sample_Plate}
#' and stored in \code{cfg$technology} by the pipeline — no manual selection needed.
#'
#' @return Absolute path of the project directory (invisibly).
#' @export
#' @import shiny shinydashboard shinyjs
#' @importFrom shinyFiles shinyDirButton shinyFilesButton shinyDirChoose
#'   shinyFileChoose getVolumes parseDirPath parseFilePaths
#' @importFrom fs path_home
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom parallel detectCores
mETHYLotest.EPIC.ProjectUI <- function(prefill_pheno = NULL,
                                       prefill_idat_dir = NULL) {

  for (pkg in c("shinyFiles", "writexl", "parallel", "shinydashboard", "shinyjs")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", pkg))
  }

  total_cores       <- parallel::detectCores()
  max_cores_allowed <- max(1L, total_cores - 2L)
  default_cores     <- max(1L, floor(max_cores_allowed / 2L))

  ILLUMINA_COLS <- c("Sample_Name", "Sample_Plate", "Sample_Group",
                     "Pool_ID", "Project", "Sample_Well",
                     "Sentrix_ID", "Sentrix_Position")
  REQUIRED_COLS <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position",
                     "Sample_Plate", "Sample_Group")

  fmt_bool <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"
  fmt_vec  <- function(v) {
    if (is.null(v) || length(v) == 0L) return("NULL")
    paste0('c("', paste(v, collapse = '", "'), '")')
  }
  null_or <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b

  # =========================================================================
  # UI
  # =========================================================================
  ui <- shinydashboard::dashboardPage(
    skin = "purple",
    shinydashboard::dashboardHeader(title = "mETHYLotest EPIC Setup"),

    shinydashboard::dashboardSidebar(
      width = 300,
      div(style = "padding:15px;",
          textInput("project_name", "Project Name",
                    placeholder = "My_EPIC_Project"),
          strong("Base Directory:"),
          splitLayout(
            cellWidths = c("75%", "25%"),
            textInput("base_dir", label = NULL, value = getwd()),
            shinyFiles::shinyDirButton("btn_base_dir", "Browse", "Select",
                                       icon  = icon("folder-open"),
                                       class = "btn-primary btn-sm")
          ),
          uiOutput("project_status_alert_sidebar")
      ),
      shinydashboard::sidebarMenu(
        id = "tabs",
        shinydashboard::menuItem("1. Data Source",
                                 tabName = "data", icon = icon("database")),
        shinydashboard::menuItem("2. Import & QC",
                                 tabName = "import_qc", icon = icon("filter")),
        shinydashboard::menuItem("3. Imputation & Normalisation",
                                 tabName = "impute_norm", icon = icon("sliders-h")),
        shinydashboard::menuItem("4. Corrections",
                                 tabName = "corrections", icon = icon("balance-scale")),
        shinydashboard::menuItem("5. Downstream Analyses",
                                 tabName = "downstream", icon = icon("chart-bar")),
        shinydashboard::menuItem("6. Review & Launch",
                                 tabName = "launch", icon = icon("rocket"))
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
        .ss-table { width:100%; border-collapse:collapse;
                    color:#333; background:#fff; }
        .ss-table th { padding:8px 10px; border-bottom:2px solid #999;
                       font-weight:600; }
        .ss-table td { padding:6px 10px;
                       border-bottom:1px solid #e0e0e0; }
        .ss-table tr:nth-child(even) td { background:#f9f9f9; }
        .tech-badge {
          display:inline-block; padding:3px 8px; border-radius:4px;
          font-size:11px; font-weight:600; background:#8e44ad;
          color:#fff; margin-left:6px;
        }
      "))),

      shinydashboard::tabItems(

        # ==================================================================
        # Tab 1 : Data Source
        # ==================================================================
        shinydashboard::tabItem(
          tabName = "data",
          fluidRow(
            shinydashboard::box(
              title = "Input Data", status = "primary",
              solidHeader = TRUE, width = 12,
              uiOutput("locked_warning_msg"),

              h4(icon("file-excel"), " Sample Sheet"),
              p("Select the Excel file (.xlsx / .xls) describing your samples.",
                "This file should follow the standard Illumina sample sheet",
                " format."),

              div(
                style = paste0("border:1px solid #ddd; border-radius:5px; ",
                               "padding:15px; margin-bottom:15px;",
                               " background:#fff;"),
                strong(icon("table"),
                       " Standard Illumina Sample Sheet Columns:"),
                tags$table(
                  class = "ss-table",
                  style = "margin-top:10px; font-size:12px;",
                  tags$thead(tags$tr(
                    tags$th("Column"), tags$th("Description"),
                    tags$th(style = "text-align:center;", "Required"))),
                  tags$tbody(
                    tags$tr(
                      tags$td(tags$code("Sample_Name")),
                      tags$td("Unique sample identifier"),
                      tags$td(style = "text-align:center;",
                              icon("check-circle",
                                   style = "color:#27ae60;"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Sample_Plate")),
                      tags$td(HTML(paste0(
                        "Array platform: <b>450K</b>, <b>EPICv1</b>, <b>EPICv2</b>",
                        " or <b>Mouse</b>.",
                        " Automatically inferred by the pipeline as",
                        " <code>cfg$technology</code>."))),
                      tags$td(style = "text-align:center;",
                              icon("check-circle",
                                   style = "color:#27ae60;"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Sample_Group")),
                      tags$td(HTML(paste0(
                        "Biological group (e.g. CTL, CASE).",
                        " Default comparison variable for DMP/DMR."))),
                      tags$td(style = "text-align:center;",
                              icon("check-circle",
                                   style = "color:#27ae60;"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Sentrix_ID")),
                      tags$td("BeadChip barcode (e.g. 207127940107)"),
                      tags$td(style = "text-align:center;",
                              icon("check-circle",
                                   style = "color:#27ae60;"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Sentrix_Position")),
                      tags$td("Position on chip (e.g. R06C01)"),
                      tags$td(style = "text-align:center;",
                              icon("check-circle",
                                   style = "color:#27ae60;"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Pool_ID")),
                      tags$td("Pool identifier (if applicable)"),
                      tags$td(style = "text-align:center; color:#999;",
                              tags$em("Optional"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Project")),
                      tags$td("Project name"),
                      tags$td(style = "text-align:center; color:#999;",
                              tags$em("Optional"))
                    ),
                    tags$tr(
                      tags$td(tags$code("Sample_Well")),
                      tags$td("Well position on plate"),
                      tags$td(style = "text-align:center; color:#999;",
                              tags$em("Optional"))
                    )
                  )
                )
              ),

              div(
                class = "alert alert-warning",
                style = "font-size:12px; padding:10px; margin-bottom:12px;",
                icon("exclamation-triangle"),
                strong(" Sample_Plate: "),
                "must contain ", tags$code("450K"), ", ",
                tags$code("EPICv1"), ", ", tags$code("EPICv2"),
                " or ", tags$code("Mouse"), ".",
                br(),
                icon("info-circle"),
                " The array type is ", strong("automatically detected"),
                " by the pipeline — no manual selection is required."
              ),

              div(
                style = paste0("border:1px solid #d5e8d4; border-radius:5px;",
                               " padding:12px; margin-bottom:15px;",
                               " background:#f0faf0;"),
                icon("lightbulb", style = "color:#27ae60;"),
                strong(" Additional columns "),
                "(e.g. ", tags$code("Batch"), ", ", tags$code("Sex"),
                ", ", tags$code("Age"), ", ...)",
                " will be detected and made available for:",
                tags$ul(style = "margin:6px 0 0 0;",
                        tags$li("Batch correction (ComBat)"),
                        tags$li("Selection as analysis variable for DMP/DMR"))
              ),

              splitLayout(
                cellWidths = c("85%", "15%"),
                textInput("pheno_path", label = NULL,
                          placeholder = "/path/to/sample_sheet.xlsx"),
                shinyFiles::shinyFilesButton("btn_pheno_path", "Browse",
                                             "Select", multiple = FALSE,
                                             icon  = icon("file-excel"),
                                             class = "btn-primary")
              ),
              br(),
              h4(icon("folder-open"), " IDAT Files Directory"),
              p("Directory containing all .idat files."),
              splitLayout(
                cellWidths = c("85%", "15%"),
                textInput("idat_dir", label = NULL,
                          placeholder = "/path/to/idat/"),
                shinyFiles::shinyDirButton("btn_idat_dir", "Browse",
                                           "Select",
                                           icon  = icon("folder-open"),
                                           class = "btn-primary")
              ),
              br(),
              actionButton("load_check_btn",
                           "Load & Validate Sample Sheet",
                           icon  = icon("check-circle"),
                           class = "btn-success btn-lg"),
              hr(),
              uiOutput("load_status_msg")
            )
          )
        ),

        # ==================================================================
        # Tab 2 : Import & QC
        # ==================================================================
        shinydashboard::tabItem(
          tabName = "import_qc",
          fluidRow(
            shinydashboard::box(
              title = "Quality Control Filters (champ.filter)",
              status = "warning", solidHeader = TRUE, width = 6,

              div(class = "param-section",
                  strong("Detection P-value"), br(), br(),
                  checkboxInput("filter_det_p",
                                "Filter failed positions (filterDetP)",
                                value = TRUE),
                  conditionalPanel(
                    "input.filter_det_p == true",
                    numericInput("det_p_cut",
                                 "Probe cutoff (detPcut)",
                                 value = 0.01, min = 0,
                                 max = 1, step = 0.001),
                    numericInput("sample_det_p_cut",
                                 "Sample failure rate (SampleCutoff)",
                                 value = 0.1, min = 0,
                                 max = 1, step = 0.01),
                    helpText("Samples exceeding this fraction of failed",
                             " probes are removed.")
                  )
              ),

              div(class = "param-section",
                  strong("Probe Filters"), br(), br(),
                  checkboxInput("filter_beads",
                                "Low bead-count (filterBeads)",
                                value = TRUE),
                  conditionalPanel(
                    "input.filter_beads == true",
                    numericInput("bead_cutoff",
                                 "Bead fraction cutoff",
                                 value = 0.05, min = 0,
                                 max = 1, step = 0.01)
                  ),
                  checkboxInput("filter_no_cg",
                                "Non-CpG probes (filterNoCG)",
                                value = TRUE),
                  checkboxInput("filter_snps",
                                "SNP-associated probes (filterSNPs)",
                                value = TRUE),
                  checkboxInput("filter_multi_hit",
                                "Multi-hit probes (filterMultiHit)",
                                value = TRUE),
                  checkboxInput("filter_xy",
                                "Sex-chromosome probes (filterXY)",
                                value = FALSE)
              )
            ),

            shinydashboard::box(
              title = "Import Options (champ.load)",
              status = "info", solidHeader = TRUE, width = 6,

              div(class = "param-section",
                  selectInput("load_method", "Loading Method",
                              choices = c(
                                "ChAMP (fast, no detP matrix)" = "ChAMP",
                                "minfi (slower, provides detP)"= "minfi"),
                              selected = "ChAMP"),
                  helpText(icon("info-circle"),
                           "minfi required for detection p-value QC plots",
                           " and SWAN normalisation.")
              ),

              div(class = "param-section",
                  checkboxInput("load_force",
                                "Force loading (force = TRUE)",
                                value = FALSE),
                  helpText("Allow loading if some sample sheet checks",
                           " fail.")
              ),

              div(class = "param-section",
                  checkboxInput("load_autoimpute",
                                "Auto-impute missing values at import",
                                value = TRUE),
                  conditionalPanel(
                    "input.load_autoimpute == true",
                    selectInput("load_imputation_method",
                                "Imputation method",
                                choices  = c("KNN"  = "knn",
                                             "Mean" = "mean"),
                                selected = "knn")
                  )
              ),

              div(class = "param-section",
                  conditionalPanel(
                    "input.filter_snps == true",
                    selectInput("load_population",
                                "SNP population filter",
                                choices = c(
                                  "None (all SNPs)" = "",
                                  "African (AFR)"   = "AFR",
                                  "American (AMR)"  = "AMR",
                                  "Asian (ASN)"     = "ASN",
                                  "European (EUR)"  = "EUR",
                                  "Mixed (MXL)"     = "MXL"),
                                selected = ""),
                    helpText("Population-specific SNP filter.",
                             " Empty = all known SNPs.")
                  )
              )
            )
          )
        ),

        # ==================================================================
        # Tab 3 : Imputation & Normalisation
        # ==================================================================
        shinydashboard::tabItem(
          tabName = "impute_norm",
          fluidRow(

            shinydashboard::box(
              title = "1. Imputation (champ.impute)",
              status = "warning", solidHeader = TRUE, width = 6,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px; margin-bottom:10px;",
                  icon("info-circle"),
                  tags$code("champ.impute()"),
                  " handles remaining missing values in the beta matrix ",
                  strong("before"), " normalisation."),

              checkboxInput("do_impute",
                            span(icon("wrench"),
                                 strong(" Run champ.impute")),
                            value = FALSE),

              conditionalPanel(
                "input.do_impute == true",
                div(class = "param-section",
                    selectInput("impute_method", "Method",
                                choices  = c(
                                  "K-Nearest Neighbours (KNN)" = "KNN",
                                  "Combined (Combine)"         = "Combine"),
                                selected = "Combine"),
                    helpText(
                      tags$b("KNN:"),
                      " Imputes using K nearest probes/samples.", br(),
                      tags$b("Combine:"),
                      " Combined strategy for robust imputation."),
                    numericInput("impute_k",
                                 "Number of neighbours (k)",
                                 value = 5, min = 1,
                                 max = 50, step = 1),
                    numericInput("impute_probe_cutoff",
                                 "Max NA fraction per probe (ProbeCutoff)",
                                 value = 0.2, min = 0,
                                 max = 1, step = 0.05),
                    helpText("Probes exceeding this fraction are removed."),
                    numericInput("impute_sample_cutoff",
                                 "Max NA fraction per sample (SampleCutoff)",
                                 value = 0.1, min = 0,
                                 max = 1, step = 0.05),
                    helpText("Samples exceeding this fraction are removed.")
                )
              )
            ),

            shinydashboard::box(
              title = "2. Normalisation (champ.norm)",
              status = "primary", solidHeader = TRUE, width = 6,

              selectInput("norm_method", "Normalisation Method",
                          choices  = c("BMIQ", "PBC", "SWAN",
                                       "illumina", "none"),
                          selected = "BMIQ"),
              uiOutput("norm_method_info"),
              hr(),
              checkboxInput("plot_norm",
                            "Save normalisation density plots",
                            value = FALSE),
              hr(),
              numericInput("norm_cores", "Cores for normalisation",
                           value = default_cores, min = 1,
                           max = max_cores_allowed, step = 1),
              helpText("BMIQ can be parallelised across samples.")
            )
          )
        ),

        # ==================================================================
        # Tab 4 : Corrections
        # ==================================================================
        shinydashboard::tabItem(
          tabName = "corrections",

          fluidRow(
            shinydashboard::box(
              title = tagList(icon("tint"),
                              " Cell Type Correction (champ.refbase)"),
              status = "warning", solidHeader = TRUE, width = 12,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px; margin-bottom:10px;",
                  icon("info-circle"),
                  strong(" Blood samples only."), br(),
                  tags$code("champ.refbase()"),
                  " estimates cell type proportions (CD4T, CD8T, NK,",
                  " Bcell, Mono, Gran) using the Houseman algorithm with",
                  " a blood reference dataset.",
                  br(), br(),
                  tags$small(
                    icon("exclamation-triangle"),
                    " Not applicable to non-blood tissues.", br(),
                    strong("Required: "),
                    tags$code("FlowSorted.Blood.450k"),
                    " reference package.")
              ),

              checkboxInput("do_refbase",
                            span(icon("tint"),
                                 strong(" Enable cell type correction",
                                        " (blood only)")),
                            value = FALSE),
              conditionalPanel(
                "input.do_refbase == true",
                div(class = "param-section",
                    checkboxInput("refbase_save_plots",
                                  "Save cell type proportion plots",
                                  value = TRUE)
                )
              )
            )
          ),

          fluidRow(
            shinydashboard::box(
              title = tagList(icon("layer-group"),
                              " Batch Correction (champ.runCombat)"),
              status = "info", solidHeader = TRUE, width = 12,
              uiOutput("batch_correction_ui")
            )
          )
        ),

        # ==================================================================
        # Tab 5 : Downstream Analyses
        # ==================================================================
        shinydashboard::tabItem(
          tabName = "downstream",

          div(class = "alert alert-info",
              style = paste0("margin:10px; font-size:12px; padding:10px;",
                             " border-radius:5px;"),
              icon("microchip"),
              strong(" Array type (technology)"),
              " is detected automatically from ",
              tags$code("Sample_Plate"),
              " and stored in ",
              tags$code("cfg$technology"),
              " by the pipeline. No manual selection is needed here."),

          # Analysis variable
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("flask"), " Analysis Variable"),
              status = "primary", solidHeader = TRUE, width = 12,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px; margin-bottom:10px;",
                  icon("info-circle"),
                  "This variable is used as ", tags$code("pheno"),
                  " in champ.DMP, champ.DMR and champ.Block.", br(),
                  "It can be categorical (e.g. Sample_Group) or",
                  " numeric (e.g. Age)."),

              fluidRow(
                column(4,
                       selectInput("analysis_pheno",
                                   "Analysis variable (pheno)",
                                   choices  = c("Sample_Group"),
                                   selected = "Sample_Group")
                ),
                column(8, uiOutput("compare_group_ui"))
              )
            )
          ),

          # DMP
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("map-pin"), " DMP (champ.DMP)"),
              status = "primary", solidHeader = TRUE, width = 12,
              collapsible = TRUE,

              checkboxInput("do_dmp",
                            span(strong(" Differentially Methylated",
                                        " Positions")),
                            value = TRUE),

              conditionalPanel(
                "input.do_dmp == true",
                div(class = "param-section",
                    fluidRow(
                      column(6,
                             numericInput("dmp_adj_p_val",
                                          "Adjusted P-value (adjPVal)",
                                          value = 0.05, min = 0,
                                          max = 1, step = 0.01)
                      ),
                      column(6,
                             selectInput("dmp_adjust_method",
                                         "Adjustment method (adjust.method)",
                                         choices = c("BH", "bonferroni",
                                                     "holm", "hochberg",
                                                     "hommel", "BY",
                                                     "fdr", "none"),
                                         selected = "BH")
                      )
                    )
                )
              )
            )
          ),

          # DMR
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("ruler-combined"),
                              " DMR (champ.DMR)"),
              status = "success", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,

              checkboxInput("do_dmr",
                            span(strong(" Differentially Methylated",
                                        " Regions")),
                            value = FALSE),

              conditionalPanel(
                "input.do_dmr == true",
                div(class = "param-section",
                    fluidRow(
                      column(3,
                             selectInput("dmr_method", "Method",
                                         choices = c("Bumphunter",
                                                     "DMRcate",
                                                     "ProbeLasso"),
                                         selected = "Bumphunter")
                      ),
                      column(3,
                             numericInput("dmr_min_probes",
                                          "Min probes (minProbes)",
                                          value = 7, min = 2,
                                          max = 50, step = 1)
                      ),
                      column(3,
                             numericInput("dmr_adj_p_val",
                                          "Adj. P-value (adjPval)",
                                          value = 0.05, min = 0,
                                          max = 1, step = 0.01)
                      ),
                      column(3,
                             numericInput("dmr_cores", "Cores",
                                          value = default_cores, min = 1,
                                          max = max_cores_allowed, step = 1)
                      )
                    ),

                    conditionalPanel(
                      "input.dmr_method == 'Bumphunter'",
                      hr(style = "margin:8px 0;"),
                      strong("Bumphunter parameters"),
                      fluidRow(
                        column(4,
                               numericInput("dmr_bh_cutoff", "Cutoff",
                                            value = 0.1, min = 0,
                                            max = 1, step = 0.01),
                               helpText("Effect-size cutoff for",
                                        " candidate regions.")
                        ),
                        column(4,
                               numericInput("dmr_bh_max_gap",
                                            "Max gap (maxGap, bp)",
                                            value = 300, min = 50,
                                            max = 5000, step = 50)
                        ),
                        column(4,
                               numericInput("dmr_bh_B",
                                            "Permutations (B)",
                                            value = 250, min = 0,
                                            max = 10000, step = 50)
                        )
                      ),
                      checkboxInput("dmr_bh_smooth",
                                    "Smoothing (smooth)", value = TRUE),
                      checkboxInput("dmr_bh_pick_cutoff",
                                    "Auto-pick cutoff (pickCutoff)",
                                    value = TRUE)
                    ),

                    conditionalPanel(
                      "input.dmr_method == 'DMRcate'",
                      hr(style = "margin:8px 0;"),
                      strong("DMRcate parameters"),
                      fluidRow(
                        column(4,
                               numericInput("dmr_cate_lambda",
                                            "Lambda (bandwidth)",
                                            value = 1000, min = 100,
                                            max = 10000, step = 100)
                        ),
                        column(4,
                               numericInput("dmr_cate_c",
                                            "Scaling factor (C)",
                                            value = 2, min = 1,
                                            max = 10, step = 1)
                        ),
                        column(4,
                               numericInput("dmr_cate_dist",
                                            "Merging distance (dist, bp)",
                                            value = 1000, min = 100,
                                            max = 10000, step = 100)
                        )
                      )
                    ),

                    conditionalPanel(
                      "input.dmr_method == 'ProbeLasso'",
                      hr(style = "margin:8px 0;"),
                      strong("ProbeLasso parameters"),
                      fluidRow(
                        column(4,
                               numericInput("dmr_pl_min_sep",
                                            "Min DMR separation (bp)",
                                            value = 1000, min = 100,
                                            max = 50000, step = 100)
                        ),
                        column(4,
                               numericInput("dmr_pl_min_size",
                                            "Min DMR size (bp)",
                                            value = 50, min = 10,
                                            max = 10000, step = 10)
                        ),
                        column(4,
                               numericInput("dmr_pl_adj_p_probe",
                                            "Probe adj. P-value",
                                            value = 0.05, min = 0,
                                            max = 1, step = 0.01)
                        )
                      )
                    )
                )
              )
            )
          ),

          # Block
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("th-large"),
                              " Block (champ.Block)"),
              status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px; margin-bottom:10px;",
                  icon("info-circle"),
                  "Identifies large-scale (>= 5 kb) blocks of",
                  " hypo/hyper-methylation using open-sea probes.",
                  br(),
                  tags$small(
                    "Parameters: ",
                    tags$code("maxClusterGap"),
                    " (max gap between probes to belong to the same block),",
                    " ", tags$code("bpSpan"),
                    " (smoothing span), ",
                    tags$code("minNum"),
                    " (min probes per block),",
                    " ", tags$code("B"),
                    " (permutations)."
                  )),

              checkboxInput("do_block",
                            span(strong(" Block Methylation Analysis")),
                            value = FALSE),

              conditionalPanel(
                "input.do_block == true",
                div(class = "param-section",
                    fluidRow(
                      column(3,
                             numericInput(
                               "block_max_cluster_gap",
                               "Max cluster gap / bpSpan (bp)",
                               value = 250000, min = 10000,
                               max = 2000000, step = 10000),
                             helpText("Used as both ",
                                      tags$code("maxClusterGap"),
                                      " and ",
                                      tags$code("bpSpan"), ".")
                      ),
                      column(3,
                             numericInput("block_min_num",
                                          "Min probes (minNum)",
                                          value = 5, min = 2,
                                          max = 100, step = 1)
                      ),
                      column(3,
                             numericInput("block_B",
                                          "Permutations (B)",
                                          value = 500, min = 0,
                                          max = 10000, step = 50)
                      ),
                      column(3,
                             numericInput("block_cores", "Cores",
                                          value = default_cores, min = 1,
                                          max = max_cores_allowed, step = 1)
                      )
                    )
                )
              )
            )
          ),

          # GSEA
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("project-diagram"),
                              " GSEA (champ.GSEA)"),
              status = "warning", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px; margin-bottom:10px;",
                  icon("info-circle"),
                  "Gene Set Enrichment Analysis on DMP results."),

              checkboxInput("do_gsea",
                            span(strong(" Gene Set Enrichment Analysis")),
                            value = FALSE),

              conditionalPanel(
                "input.do_gsea == true",
                uiOutput("gsea_dmp_warning"),
                div(class = "param-section",
                    fluidRow(
                      column(6,
                             numericInput("gsea_adj_p_val",
                                          "Adj. P-value (adjPval)",
                                          value = 0.05, min = 0,
                                          max = 1, step = 0.01)
                      ),
                      column(6,
                             selectInput("gsea_method", "Method",
                                         choices = c("fisher", "gometh"),
                                         selected = "fisher"),
                             helpText(
                               tags$b("fisher:"),
                               " Fisher's exact test.", br(),
                               tags$b("gometh:"),
                               " Adjusts for probe count bias.")
                      )
                    )
                )
              )
            )
          ),

          # CNA
          fluidRow(
            shinydashboard::box(
              title = tagList(icon("dna"), " CNA (champ.CNA)"),
              status = "danger", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,

              div(class = "alert alert-info",
                  style = "font-size:12px; padding:10px; margin-bottom:10px;",
                  icon("info-circle"),
                  "Copy Number Aberration detection from methylation",
                  " intensity data.",
                  br(),
                  tags$small(icon("exclamation-triangle"),
                             " Requires ",
                             tags$code("load_method = 'minfi'"),
                             " to obtain intensity data.")
              ),

              checkboxInput("do_cna",
                            span(strong(" Copy Number Analysis")),
                            value = FALSE),

              conditionalPanel(
                "input.do_cna == true",
                div(class = "param-section",
                    fluidRow(
                      column(4,
                             textInput("cna_control_group",
                                       "Control group (controlGroup)",
                                       value = "CTL"),
                             helpText("Sample_Group label for",
                                      " reference samples.")
                      ),
                      column(4,
                             numericInput("cna_freq_threshold",
                                          "Freq. threshold (freqThreshold)",
                                          value = 0.3, min = 0,
                                          max = 1, step = 0.05)
                      ),
                      column(4,
                             selectInput("cna_genome_build",
                                         "Genome build",
                                         choices = c("hg19", "hg38"),
                                         selected = "hg19")
                      )
                    ),
                    fluidRow(
                      column(4,
                             checkboxInput("cna_sample_cna",
                                           "Per-sample CNA plots",
                                           value = TRUE)
                      ),
                      column(4,
                             checkboxInput("cna_group_freq_plots",
                                           "Group frequency plots",
                                           value = TRUE)
                      ),
                      column(4,
                             numericInput("cna_cores", "Cores",
                                          value = default_cores, min = 1,
                                          max = max_cores_allowed, step = 1)
                      )
                    )
                )
              )
            )
          )
        ),

        # ==================================================================
        # Tab 6 : Review & Launch
        # ==================================================================
        shinydashboard::tabItem(tabName = "launch",
                                uiOutput("launch_ui_dynamic"))
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

    # All lockable IDs — no *_array_type entries
    LOCK_IDS <- c(
      "project_name", "base_dir",
      "pheno_path", "idat_dir",
      "btn_pheno_path", "btn_idat_dir", "btn_base_dir",
      "load_check_btn",
      "filter_det_p", "det_p_cut", "sample_det_p_cut",
      "filter_beads", "bead_cutoff",
      "filter_no_cg", "filter_snps", "filter_multi_hit", "filter_xy",
      "load_method", "load_force", "load_autoimpute",
      "load_imputation_method", "load_population",
      "do_impute", "impute_method", "impute_k",
      "impute_probe_cutoff", "impute_sample_cutoff",
      "norm_method", "plot_norm", "norm_cores",
      "do_refbase", "refbase_save_plots",
      "do_batch_correction", "combat_logit_transform",
      "analysis_pheno", "do_specific_comparison",
      "compare_group_1", "compare_group_2",
      "do_dmp", "dmp_adj_p_val", "dmp_adjust_method",
      "do_dmr", "dmr_method", "dmr_min_probes", "dmr_adj_p_val",
      "dmr_cores", "dmr_bh_cutoff", "dmr_bh_max_gap", "dmr_bh_B",
      "dmr_bh_smooth", "dmr_bh_pick_cutoff",
      "dmr_cate_lambda", "dmr_cate_c", "dmr_cate_dist",
      "dmr_pl_min_sep", "dmr_pl_min_size", "dmr_pl_adj_p_probe",
      "do_block", "block_max_cluster_gap", "block_min_num",
      "block_B", "block_cores",
      "do_gsea", "gsea_adj_p_val", "gsea_method",
      "do_cna", "cna_control_group", "cna_freq_threshold",
      "cna_genome_build",
      "cna_sample_cna", "cna_group_freq_plots", "cna_cores",
      "num_cores", "save_raw_obj"
    )

    # ══════════════════════════════════════════════════════════════
    # Auto-fill paths (demo mode) — user still clicks Load
    # ══════════════════════════════════════════════════════════════
    observe({
      filled <- FALSE
      if (!is.null(prefill_pheno) && file.exists(prefill_pheno)) {
        updateTextInput(session, "pheno_path", value = prefill_pheno)
        filled <- TRUE
      }
      if (!is.null(prefill_idat_dir) && dir.exists(prefill_idat_dir)) {
        updateTextInput(session, "idat_dir", value = prefill_idat_dir)
        filled <- TRUE
      }
      if (filled) {
        showNotification(
          "Demo paths set. Click 'Load & Validate Sample Sheet' to continue.",
          type = "message", duration = 8)
      }
    })

    output$norm_method_info <- renderUI({
      req(input$norm_method)
      if (input$norm_method == "SWAN" &&
          !is.null(input$load_method) &&
          input$load_method == "ChAMP") {
        div(class = "alert alert-warning",
            style = "font-size:11px; padding:6px; margin-top:4px;",
            icon("exclamation-triangle"),
            " SWAN requires minfi loading method.")
      } else if (input$norm_method == "none") {
        div(class = "alert alert-info",
            style = "font-size:11px; padding:6px; margin-top:4px;",
            icon("info-circle"),
            " No normalisation. Raw beta values will be used.")
      } else {
        helpText(icon("info-circle"),
                 switch(input$norm_method,
                        BMIQ     = paste("Beta-Mixture Quantile.",
                                         "Recommended for EPIC arrays."),
                        PBC      = "Peak-Based Correction.",
                        SWAN     = paste("Subset Within Array Normalisation",
                                         "(requires minfi loading)."),
                        illumina = "Illumina internal normalisation.",
                        ""))
      }
    })

    output$gsea_dmp_warning <- renderUI({
      if (!isTRUE(input$do_dmp))
        div(class = "alert alert-danger",
            style = "font-size:12px; padding:8px; margin-bottom:10px;",
            icon("exclamation-triangle"),
            strong(" DMP analysis is not enabled."),
            " GSEA requires DMP results as input.")
    })

    observe({
      cols <- pheno_columns()
      if (!is.null(cols))
        updateSelectInput(session, "analysis_pheno",
                          choices  = c("Sample_Group", cols),
                          selected = "Sample_Group")
    })

    output$compare_group_ui <- renderUI({
      df        <- loaded_data()
      pheno_col <- input$analysis_pheno
      if (is.null(df) || is.null(pheno_col) ||
          !pheno_col %in% colnames(df)) return(NULL)

      vals <- unique(as.character(df[[pheno_col]]))
      if (length(vals) > 20)
        return(helpText(icon("info-circle"),
                        "High-cardinality / numeric variable.",
                        " Pairwise comparison not applicable."))
      if (length(vals) < 2)
        return(helpText(icon("exclamation-triangle"),
                        "Only one level detected."))

      tagList(
        checkboxInput("do_specific_comparison",
                      "Specific pairwise comparison (compare.group)",
                      value = FALSE),
        conditionalPanel(
          "input.do_specific_comparison == true",
          fluidRow(
            column(6,
                   selectInput("compare_group_1",
                               "Group 1 (reference)",
                               choices = vals, selected = vals[1])
            ),
            column(6,
                   selectInput("compare_group_2",
                               "Group 2 (comparison)",
                               choices = vals,
                               selected = if (length(vals) >= 2)
                                 vals[2] else vals[1])
            )
          ),
          helpText("Leave unchecked for all pairwise comparisons.")
        )
      )
    })

    output$batch_correction_ui <- renderUI({
      cols <- pheno_columns()
      if (is.null(cols) || length(cols) == 0L)
        return(div(class = "alert alert-warning",
                   icon("exclamation-circle"),
                   " No extra columns available. Add columns (e.g. ",
                   tags$code("Batch"), ", ", tags$code("Slide"),
                   ") to enable batch correction."))

      tagList(
        checkboxInput("do_batch_correction",
                      span(icon("layer-group"),
                           strong(" Enable ComBat Batch Correction")),
                      value = FALSE),
        conditionalPanel(
          "input.do_batch_correction == true",
          div(class = "param-section",
              checkboxGroupInput("batch_cols",
                                 "Batch variable(s) (batchname):",
                                 choices = cols, inline = TRUE),
              hr(style = "margin:8px 0;"),
              selectInput("biological_variable",
                          tagList(icon("flask"),
                                  " Variable to protect (variablename):"),
                          choices  = c("Sample_Group", cols),
                          selected = "Sample_Group"),
              helpText("ComBat preserves variance from this variable."),
              checkboxInput("combat_logit_transform",
                            "Logit-transform betas (logitTrans)",
                            value = TRUE),
              helpText("Recommended for beta-value data.")
          )
        )
      )
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
               file.exists(file.path(path, "Results", "project_config.R")))
        detected_path <- path

      if (!is.null(detected_path)) {
        if (!is_existing_project()) {
          is_existing_project(TRUE)
          updateTextInput(session, "project_name",
                          value = basename(detected_path))
          updateTextInput(session, "base_dir",
                          value = dirname(detected_path))

          cfg_path <- file.path(detected_path, "Results", "project_config.R")
          tryCatch({
            env <- new.env()
            source(cfg_path, local = env)
            cfg <- env$project_config

            .ut <- function(id, key = id)
              if (!is.null(cfg[[key]]))
                updateTextInput(session, id, value = cfg[[key]])
            .uc <- function(id, key = id)
              if (!is.null(cfg[[key]]))
                updateCheckboxInput(session, id, value = cfg[[key]])
            .un <- function(id, key = id)
              if (!is.null(cfg[[key]]))
                updateNumericInput(session, id, value = cfg[[key]])
            .us <- function(id, key = id)
              if (!is.null(cfg[[key]]))
                updateSelectInput(session, id, selected = cfg[[key]])

            .ut("pheno_path", "pheno_file")
            .ut("idat_dir")

            for (id in c("filter_det_p", "filter_beads", "filter_no_cg",
                         "filter_snps", "filter_multi_hit", "filter_xy"))
              .uc(id)
            for (id in c("det_p_cut", "sample_det_p_cut", "bead_cutoff"))
              .un(id)

            .us("load_method"); .uc("load_force")
            .uc("load_autoimpute"); .us("load_imputation_method")
            .us("load_population")

            .uc("do_impute"); .us("impute_method")
            for (id in c("impute_k", "impute_probe_cutoff",
                         "impute_sample_cutoff"))
              .un(id)

            .us("norm_method"); .uc("plot_norm"); .un("norm_cores")

            .uc("do_refbase"); .uc("refbase_save_plots")
            if (!is.null(cfg$perform_batch_correction))
              updateCheckboxInput(session, "do_batch_correction",
                                  value = cfg$perform_batch_correction)
            .uc("combat_logit_transform")

            .us("analysis_pheno")
            if (!is.null(cfg$compare_group) &&
                length(cfg$compare_group) == 2L) {
              updateCheckboxInput(session, "do_specific_comparison",
                                  value = TRUE)
              updateSelectInput(session, "compare_group_1",
                                selected = cfg$compare_group[1L])
              updateSelectInput(session, "compare_group_2",
                                selected = cfg$compare_group[2L])
            }

            .uc("do_dmp"); .un("dmp_adj_p_val")
            .us("dmp_adjust_method")

            .uc("do_dmr"); .us("dmr_method")
            for (id in c("dmr_min_probes", "dmr_adj_p_val", "dmr_cores",
                         "dmr_bh_cutoff", "dmr_bh_max_gap", "dmr_bh_B",
                         "dmr_cate_lambda", "dmr_cate_c", "dmr_cate_dist",
                         "dmr_pl_min_sep", "dmr_pl_min_size",
                         "dmr_pl_adj_p_probe"))
              .un(id)
            .uc("dmr_bh_smooth"); .uc("dmr_bh_pick_cutoff")

            .uc("do_block")
            for (id in c("block_max_cluster_gap", "block_min_num",
                         "block_B", "block_cores"))
              .un(id)

            .uc("do_gsea"); .un("gsea_adj_p_val"); .us("gsea_method")

            .uc("do_cna"); .ut("cna_control_group")
            .un("cna_freq_threshold"); .us("cna_genome_build")
            .uc("cna_sample_cna"); .uc("cna_group_freq_plots")
            .un("cna_cores")

            if (!is.null(cfg$num_cores))
              updateSliderInput(session, "num_cores",
                                value = cfg$num_cores)
            .uc("save_raw_obj")

            for (id in LOCK_IDS) shinyjs::disable(id)

            pheno_out <- file.path(detected_path, "data",
                                   "selected_samples.xlsx")
            if (file.exists(pheno_out))
              try({
                d <- readxl::read_excel(pheno_out)
                loaded_data(d)
                pheno_columns(setdiff(
                  colnames(d),
                  c("Sample_Name", "Sample_Group")))
              }, silent = TRUE)

            showNotification("Existing project loaded. Ready to re-launch.",
                             type = "message", duration = 5)
            shinydashboard::updateTabItems(session, "tabs", "launch")

          }, error = function(e)
            showNotification(paste("Error reading config:", e$message),
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

    output$project_status_alert_sidebar <- renderUI({
      if (is_existing_project())
        div(class = "callout callout-success",
            style = "padding:10px; margin-top:10px; border-left-color:#00a65a;",
            h5(icon("check-circle"), " Project Found"),
            tags$small("Configuration loaded & locked."))
      else
        div(class = "callout callout-info",
            style = "padding:10px; margin-top:10px; border-left-color:#00c0ef;",
            h5(icon("plus"), " New Project"),
            tags$small("Configure and generate."))
    })

    output$locked_warning_msg <- renderUI({
      if (is_existing_project())
        div(class = "alert alert-warning",
            icon("lock"), strong(" READ-ONLY:"),
            " Existing project. Parameters locked.")
    })

    shinyFiles::shinyDirChoose(input, "btn_base_dir",
                               roots = my_roots, session = session)
    observeEvent(input$btn_base_dir, {
      p <- shinyFiles::parseDirPath(my_roots, input$btn_base_dir)
      if (length(p) > 0)
        updateTextInput(session, "base_dir", value = as.character(p))
    })

    shinyFiles::shinyFileChoose(input, "btn_pheno_path",
                                roots = my_roots, session = session,
                                filetypes = c("xlsx", "xls"))
    observeEvent(input$btn_pheno_path, {
      f <- shinyFiles::parseFilePaths(my_roots, input$btn_pheno_path)
      if (nrow(f) > 0)
        updateTextInput(session, "pheno_path",
                        value = as.character(f$datapath))
    })

    shinyFiles::shinyDirChoose(input, "btn_idat_dir",
                               roots = my_roots, session = session)
    observeEvent(input$btn_idat_dir, {
      p <- shinyFiles::parseDirPath(my_roots, input$btn_idat_dir)
      if (length(p) > 0)
        updateTextInput(session, "idat_dir", value = as.character(p))
    })

    session$onSessionEnded(function() {
      stopApp(returnValue = isolate(final_path()))
    })

    observeEvent(input$load_check_btn, {
      req(input$pheno_path, input$idat_dir)

      if (!file.exists(input$pheno_path)) {
        showNotification("Sample sheet not found!", type = "error")
        return()
      }
      if (!dir.exists(input$idat_dir)) {
        showNotification("IDAT directory not found!", type = "error")
        return()
      }

      tryCatch({
        df      <- readxl::read_excel(input$pheno_path)
        missing <- setdiff(REQUIRED_COLS, colnames(df))

        if (length(missing) > 0) {
          showNotification(
            paste("Missing columns:", paste(missing, collapse = ", ")),
            type = "error")
          output$load_status_msg <- renderUI(
            div(class = "alert alert-danger",
                icon("exclamation-triangle"),
                " Missing: ", paste(missing, collapse = ", ")))
          loaded_data(NULL)
          return()
        }

        valid_plates <- c("450K", "EPICv1", "EPICv2", "Mouse")
        bad_plates   <- unique(
          df$Sample_Plate[!df$Sample_Plate %in% valid_plates])
        bad_pos <- df$Sentrix_Position[
          !grepl("^R[0-9]{2}C[0-9]{2}$", df$Sentrix_Position)]

        adj_cols <- setdiff(colnames(df),
                            c("Sample_Name", "Sample_Group"))
        loaded_data(df)
        pheno_columns(adj_cols)

        output$load_status_msg <- renderUI({
          tagList(
            div(class = "alert alert-success", icon("check"),
                sprintf(" Sample sheet loaded: %d samples, %d array type(s).",
                        nrow(df), length(unique(df$Sample_Plate)))),

            if (length(bad_plates) > 0)
              div(class = "alert alert-danger",
                  style = "font-size:12px; padding:8px;",
                  icon("times-circle"),
                  strong(" Invalid Sample_Plate: "),
                  paste(bad_plates, collapse = ", "),
                  ". Expected: ", paste(valid_plates, collapse = ", ")),

            if (length(bad_pos) > 0)
              div(class = "alert alert-warning",
                  style = "font-size:12px; padding:8px;",
                  icon("exclamation-triangle"),
                  strong(" Sentrix_Position format: "),
                  length(bad_pos), " entries don't match ",
                  tags$code("RxxCxx")),

            {
              extra <- setdiff(colnames(df), ILLUMINA_COLS)
              if (length(extra) > 0)
                div(class = "alert alert-info",
                    style = "font-size:12px; padding:8px;",
                    icon("info-circle"),
                    " Additional columns: ",
                    tags$strong(paste(extra, collapse = ", ")),
                    br(),
                    tags$small("Available for batch correction and",
                               " analysis variable selection."))
              else
                div(class = "alert alert-warning",
                    style = "font-size:12px; padding:8px;",
                    icon("exclamation-circle"),
                    " No extra columns. Add ",
                    tags$code("Batch"), ", ", tags$code("Sex"), ", ",
                    tags$code("Age"), " etc. to enable adjustments.")
            }
          )
        })

        if (length(bad_plates) == 0)
          shinydashboard::updateTabItems(session, "tabs", "import_qc")

      }, error = function(e)
        showNotification(paste("Error:", e$message), type = "error"))
    })

    output$launch_ui_dynamic <- renderUI({

      if (is_existing_project()) {
        return(tagList(
          div(class = "alert alert-info", style = "margin-bottom:15px;",
              icon("info-circle"),
              " Project locked. Click Re-launch to proceed."),
          fluidRow(
            shinydashboard::box(
              title = "Project Status", status = "success",
              solidHeader = TRUE, width = 12,
              div(style = "text-align:center; padding:10px;",
                  h3(input$project_name),
                  p(style = "color:gray;",
                    file.path(input$base_dir, input$project_name)),
                  hr(),
                  actionButton("relaunch_btn",
                               "Re-launch EPIC Project",
                               icon  = icon("rocket"),
                               class = "btn-primary btn-lg",
                               style = "width:50%;margin-bottom:20px;"))
            )
          ),
          if (!is.null(loaded_data()))
            fluidRow(
              shinydashboard::box(
                title = "Samples (Read-Only)", status = "primary",
                solidHeader = TRUE, width = 12,
                collapsible = TRUE, collapsed = TRUE,
                DT::DTOutput("selection_table"))
            )
        ))
      }

      df <- loaded_data()
      if (is.null(df))
        return(fluidRow(shinydashboard::box(
          width = 12, status = "danger",
          h3(icon("arrow-left"), " No Data Loaded"),
          p("Go to 'Data Source' and load your sample sheet."))))

      tagList(
        fluidRow(
          shinydashboard::box(
            title = paste("Select Samples (", nrow(df), "total)"),
            status = "primary", solidHeader = TRUE, width = 12,
            p("Click rows to select/deselect."),
            DT::DTOutput("selection_table"))
        ),

        fluidRow(
          shinydashboard::box(
            title = "Configuration Summary", status = "info",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            fluidRow(
              column(3,
                     h5(strong("Import & QC")),
                     tags$ul(
                       tags$li(paste("Method:",
                                     null_or(input$load_method, "ChAMP"))),
                       tags$li(paste("DetP:",   fmt_bool(input$filter_det_p))),
                       tags$li(paste("Beads:",  fmt_bool(input$filter_beads))),
                       tags$li(paste("SNPs:",   fmt_bool(input$filter_snps))),
                       tags$li(paste("XY:",     fmt_bool(input$filter_xy)))
                     )
              ),
              column(3,
                     h5(strong("Processing")),
                     tags$ul(
                       tags$li(paste("Impute:", fmt_bool(input$do_impute))),
                       tags$li(paste("Norm:",
                                     null_or(input$norm_method, "BMIQ"))),
                       tags$li(paste("RefBase:",
                                     fmt_bool(input$do_refbase))),
                       tags$li(paste("ComBat:",
                                     fmt_bool(input$do_batch_correction)))
                     )
              ),
              column(3,
                     h5(strong("Analyses")),
                     tags$ul(
                       tags$li(paste("Pheno:",
                                     null_or(input$analysis_pheno,
                                             "Sample_Group"))),
                       tags$li(paste("DMP:",   fmt_bool(input$do_dmp))),
                       tags$li(paste("DMR:",   fmt_bool(input$do_dmr))),
                       tags$li(paste("Block:", fmt_bool(input$do_block))),
                       tags$li(paste("GSEA:",  fmt_bool(input$do_gsea))),
                       tags$li(paste("CNA:",   fmt_bool(input$do_cna)))
                     )
              ),
              column(3,
                     h5(strong("System")),
                     tags$ul(
                       tags$li(paste("Cores:", input$num_cores)),
                       tags$li(paste("Save raw:",
                                     fmt_bool(input$save_raw_obj))),
                       tags$li(span("Technology: ",
                                    span(class = "tech-badge",
                                         "Auto-detected")))
                     )
              )
            )
          )
        ),

        fluidRow(
          shinydashboard::box(
            title = "System Resources", status = "info",
            solidHeader = TRUE, width = 12,
            fluidRow(
              column(4,
                     sliderInput("num_cores",
                                 label = paste0("Cores (Max: ",
                                                max_cores_allowed, ")"),
                                 min = 1, max = max_cores_allowed,
                                 value = default_cores, step = 1),
                     helpText(paste("Detected:", total_cores))
              ),
              column(4,
                     checkboxInput("save_raw_obj",
                                   "Save raw object (.rds)?",
                                   value = TRUE)
              )
            )
          )
        ),

        fluidRow(
          column(12, align = "center", hr(),
                 actionButton("save_generate_btn",
                              "Generate EPIC Project",
                              icon  = icon("save"),
                              class = "btn-success btn-lg",
                              style = "width:50%;margin-bottom:30px;"))
        )
      )
    })

    output$selection_table <- DT::renderDT({
      req(loaded_data())
      DT::datatable(
        loaded_data(),
        selection = list(
          mode     = if (is_existing_project()) "none" else "multiple",
          selected = seq_len(nrow(loaded_data())),
          target   = "row"),
        options = list(dom = "ftp", pageLength = 10, scrollX = TRUE))
    })

    observeEvent(input$relaunch_btn, {
      project_dir <- file.path(input$base_dir, input$project_name)
      final_path(project_dir)
      stopApp(project_dir)
    })

    observeEvent(input$save_generate_btn, {
      req(input$project_name, input$base_dir, loaded_data())

      sel <- input$selection_table_rows_selected
      if (length(sel) == 0L) {
        showNotification("Select at least one sample.", type = "error")
        return()
      }

      final_df    <- loaded_data()[sel, , drop = FALSE]
      project_dir <- file.path(input$base_dir, input$project_name)
      res_dir     <- file.path(project_dir, "Results")
      data_dir    <- file.path(project_dir, "data")
      rds_dir     <- file.path(data_dir, "interim")
      config_file <- file.path(res_dir, "project_config.R")

      if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
      if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)

      pheno_out <- file.path(data_dir, "selected_samples.xlsx")
      writexl::write_xlsx(final_df, pheno_out)

      cg_str <- if (isTRUE(input$do_specific_comparison) &&
                    !is.null(input$compare_group_1) &&
                    !is.null(input$compare_group_2) &&
                    input$compare_group_1 != input$compare_group_2)
        sprintf('c("%s", "%s")',
                input$compare_group_1, input$compare_group_2)
      else "NULL"

      batch_vec <- if (isTRUE(input$do_batch_correction))
        input$batch_cols else NULL
      bio_var <- null_or(
        if (isTRUE(input$do_batch_correction)) input$biological_variable,
        "Sample_Group")

      snp_val <- trimws(null_or(input$load_population, ""))
      snp_str <- if (!nchar(snp_val)) "NULL"
      else sprintf('"%s"', snp_val)

      lines <- c(
        "# mETHYLotest EPIC Project Configuration",
        sprintf("# Generated: %s", Sys.time()),
        sprintf("# Samples:   %d", nrow(final_df)),
        "",
        "project_config <- list()",
        "",
        "# Paths",
        sprintf('project_config$project_name <- "%s"', input$project_name),
        sprintf('project_config$project_dir  <- "%s"', project_dir),
        sprintf('project_config$res_dir      <- "%s"', res_dir),
        sprintf('project_config$rds_dir      <- "%s"', rds_dir),
        sprintf('project_config$pheno_file   <- "%s"', pheno_out),
        sprintf('project_config$idat_dir     <- "%s"', input$idat_dir),
        "",
        "# Technology — set automatically by pipeline (cfg$technology)",
        "project_config$technology <- NA",
        "",
        "# Column mapping",
        'project_config$col_sample_name  <- "Sample_Name"',
        'project_config$col_sample_plate <- "Sample_Plate"',
        'project_config$col_sample_group <- "Sample_Group"',
        'project_config$col_sentrix_id   <- "Sentrix_ID"',
        'project_config$col_sentrix_pos  <- "Sentrix_Position"',
        "",
        "# QC (champ.filter)",
        sprintf("project_config$filter_det_p     <- %s",
                fmt_bool(input$filter_det_p)),
        sprintf("project_config$det_p_cut        <- %s", input$det_p_cut),
        sprintf("project_config$sample_det_p_cut <- %s",
                null_or(input$sample_det_p_cut, 0.1)),
        sprintf("project_config$filter_beads     <- %s",
                fmt_bool(input$filter_beads)),
        sprintf("project_config$bead_cutoff      <- %s", input$bead_cutoff),
        sprintf("project_config$filter_no_cg     <- %s",
                fmt_bool(input$filter_no_cg)),
        sprintf("project_config$filter_snps      <- %s",
                fmt_bool(input$filter_snps)),
        sprintf("project_config$filter_multi_hit <- %s",
                fmt_bool(input$filter_multi_hit)),
        sprintf("project_config$filter_xy        <- %s",
                fmt_bool(input$filter_xy)),
        "",
        "# Import (champ.load)",
        sprintf('project_config$load_method            <- "%s"',
                null_or(input$load_method, "ChAMP")),
        sprintf("project_config$load_force             <- %s",
                fmt_bool(input$load_force)),
        sprintf("project_config$load_autoimpute        <- %s",
                fmt_bool(null_or(input$load_autoimpute, TRUE))),
        sprintf('project_config$load_imputation_method <- "%s"',
                null_or(input$load_imputation_method, "knn")),
        sprintf("project_config$load_population        <- %s", snp_str),
        "",
        "# Imputation (champ.impute) — before normalisation",
        sprintf("project_config$do_impute            <- %s",
                fmt_bool(input$do_impute)),
        sprintf('project_config$impute_method        <- "%s"',
                null_or(input$impute_method, "Combine")),
        sprintf("project_config$impute_k             <- %s",
                null_or(input$impute_k, 5)),
        sprintf("project_config$impute_probe_cutoff  <- %s",
                null_or(input$impute_probe_cutoff, 0.2)),
        sprintf("project_config$impute_sample_cutoff <- %s",
                null_or(input$impute_sample_cutoff, 0.1)),
        "",
        "# Normalisation (champ.norm)",
        sprintf('project_config$norm_method    <- "%s"', input$norm_method),
        sprintf("project_config$normalize_data <- %s",
                fmt_bool(input$norm_method != "none")),
        sprintf("project_config$plot_norm      <- %s",
                fmt_bool(input$plot_norm)),
        sprintf("project_config$norm_cores     <- %s",
                null_or(input$norm_cores, default_cores)),
        "",
        "# Cell type correction (champ.refbase) — blood only",
        sprintf("project_config$do_refbase         <- %s",
                fmt_bool(input$do_refbase)),
        sprintf("project_config$refbase_save_plots <- %s",
                fmt_bool(null_or(input$refbase_save_plots, TRUE))),
        "",
        "# Batch correction (champ.runCombat)",
        sprintf("project_config$perform_batch_correction <- %s",
                fmt_bool(input$do_batch_correction)),
        sprintf("project_config$batch_cols               <- %s",
                fmt_vec(batch_vec)),
        sprintf('project_config$biological_variable      <- "%s"', bio_var),
        sprintf("project_config$combat_logit_transform   <- %s",
                fmt_bool(input$combat_logit_transform)),
        "",
        "# Analysis variable (shared by DMP / DMR / Block)",
        sprintf('project_config$analysis_pheno <- "%s"',
                null_or(input$analysis_pheno, "Sample_Group")),
        sprintf("project_config$compare_group  <- %s", cg_str),
        "",
        "# DMP (champ.DMP)",
        sprintf("project_config$do_dmp            <- %s",
                fmt_bool(input$do_dmp)),
        sprintf("project_config$dmp_adj_p_val     <- %s",
                null_or(input$dmp_adj_p_val, 0.05)),
        sprintf('project_config$dmp_adjust_method <- "%s"',
                null_or(input$dmp_adjust_method, "BH")),
        "",
        "# DMR (champ.DMR)",
        sprintf("project_config$do_dmr             <- %s",
                fmt_bool(input$do_dmr)),
        sprintf('project_config$dmr_method         <- "%s"',
                null_or(input$dmr_method, "Bumphunter")),
        sprintf("project_config$dmr_min_probes     <- %s",
                null_or(input$dmr_min_probes, 7)),
        sprintf("project_config$dmr_adj_p_val      <- %s",
                null_or(input$dmr_adj_p_val, 0.05)),
        sprintf("project_config$dmr_cores          <- %s",
                null_or(input$dmr_cores, default_cores)),
        sprintf("project_config$dmr_bh_cutoff      <- %s",
                null_or(input$dmr_bh_cutoff, 0.1)),
        sprintf("project_config$dmr_bh_max_gap     <- %s",
                null_or(input$dmr_bh_max_gap, 300)),
        sprintf("project_config$dmr_bh_B           <- %s",
                null_or(input$dmr_bh_B, 250)),
        sprintf("project_config$dmr_bh_smooth      <- %s",
                fmt_bool(null_or(input$dmr_bh_smooth, TRUE))),
        sprintf("project_config$dmr_bh_pick_cutoff <- %s",
                fmt_bool(null_or(input$dmr_bh_pick_cutoff, TRUE))),
        sprintf("project_config$dmr_cate_lambda    <- %s",
                null_or(input$dmr_cate_lambda, 1000)),
        sprintf("project_config$dmr_cate_c         <- %s",
                null_or(input$dmr_cate_c, 2)),
        sprintf("project_config$dmr_cate_dist      <- %s",
                null_or(input$dmr_cate_dist, 1000)),
        sprintf("project_config$dmr_pl_min_sep     <- %s",
                null_or(input$dmr_pl_min_sep, 1000)),
        sprintf("project_config$dmr_pl_min_size    <- %s",
                null_or(input$dmr_pl_min_size, 50)),
        sprintf("project_config$dmr_pl_adj_p_probe <- %s",
                null_or(input$dmr_pl_adj_p_probe, 0.05)),
        "",
        "# Block (champ.Block)",
        "# block_max_cluster_gap is used for both maxClusterGap and bpSpan",
        sprintf("project_config$do_block             <- %s",
                fmt_bool(input$do_block)),
        sprintf("project_config$block_max_cluster_gap <- %s",
                null_or(input$block_max_cluster_gap, 250000)),
        sprintf("project_config$block_min_num        <- %s",
                null_or(input$block_min_num, 5)),
        sprintf("project_config$block_B              <- %s",
                null_or(input$block_B, 500)),
        sprintf("project_config$block_cores          <- %s",
                null_or(input$block_cores, default_cores)),
        "",
        "# GSEA (champ.GSEA)",
        sprintf("project_config$do_gsea        <- %s",
                fmt_bool(input$do_gsea)),
        sprintf("project_config$gsea_adj_p_val <- %s",
                null_or(input$gsea_adj_p_val, 0.05)),
        sprintf('project_config$gsea_method    <- "%s"',
                null_or(input$gsea_method, "fisher")),
        "",
        "# CNA (champ.CNA)",
        sprintf("project_config$do_cna               <- %s",
                fmt_bool(input$do_cna)),
        sprintf('project_config$cna_control_group    <- "%s"',
                null_or(input$cna_control_group, "CTL")),
        sprintf("project_config$cna_freq_threshold   <- %s",
                null_or(input$cna_freq_threshold, 0.3)),
        sprintf('project_config$cna_genome_build     <- "%s"',
                null_or(input$cna_genome_build, "hg19")),
        sprintf("project_config$cna_sample_cna       <- %s",
                fmt_bool(null_or(input$cna_sample_cna, TRUE))),
        sprintf("project_config$cna_group_freq_plots <- %s",
                fmt_bool(null_or(input$cna_group_freq_plots, TRUE))),
        sprintf("project_config$cna_cores            <- %s",
                null_or(input$cna_cores, default_cores)),
        "",
        "# System",
        sprintf("project_config$num_cores    <- %s", input$num_cores),
        sprintf("project_config$save_raw_obj <- %s",
                fmt_bool(input$save_raw_obj))
      )

      tryCatch({
        writeLines(lines, config_file)
        final_path(project_dir)
        showNotification("EPIC project created.", type = "message")
      }, error = function(e)
        showNotification(paste("Save failed:", e$message), type = "error"))
    })
  }

  runApp(shinyApp(ui, server), host = "0.0.0.0", launch.browser = TRUE)
}
