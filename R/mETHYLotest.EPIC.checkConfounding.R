#' Check for confounding before running ComBat
#'
#' @description
#' Checks each batch variable for confounding with the biological variable:
#' \enumerate{
#'   \item Single level: batch has < 2 levels (nothing to correct).
#'   \item Perfect aliasing: a biological group exists in only one batch level.
#'   \item Total collinearity: ALL batch levels contain only one bio group.
#'   \item Singleton batches: a batch level has only 1 sample.
#' }
#'
#' Partial collinearity (some batch levels contain only one group but
#' not all) is reported as a WARNING but does NOT block ComBat.
#'
#' @param pd Data frame of phenotype data.
#' @param variablename Character. The biological variable to protect.
#' @param batchname Character vector. The batch variable(s) to check.
#'
#' @return A named list with:
#'   \describe{
#'     \item{confounded}{Character vector of confounded batch variables.}
#'     \item{clean}{Character vector of batch variables safe for ComBat.}
#'     \item{warnings}{Character vector of variables with partial issues.}
#'     \item{details}{Named list of diagnostic info per variable.}
#'   }
#'
#' @export
mETHYLotest.EPIC.checkConfounding <- function(pd, variablename, batchname) {

  if (is.null(variablename) || is.null(batchname) ||
      length(batchname) == 0L)
    return(list(confounded = character(0),
                clean      = character(0),
                warnings   = character(0),
                details    = list()))

  variablename <- as.character(variablename)[1L]
  batchname    <- as.character(batchname)

  if (!variablename %in% colnames(pd)) {
    warning("[mETHYLotest Confounding] Biological variable '",
            variablename, "' not found in pd.")
    return(list(confounded = character(0),
                clean      = character(0),
                warnings   = character(0),
                details    = list()))
  }

  confounded_list <- character(0)
  clean_list      <- character(0)
  warn_list       <- character(0)
  details_list    <- list()

  pd[[variablename]] <- droplevels(as.factor(pd[[variablename]]))

  message("[mETHYLotest Confounding] Checking ", length(batchname),
          " batch variable(s) against '", variablename, "'")
  message("[mETHYLotest Confounding] Biological groups: ",
          paste(levels(pd[[variablename]]), collapse = ", "),
          " (", nlevels(pd[[variablename]]), " levels)")
  message(strrep("-", 60))

  for (b_name in batchname) {

    diag <- list(variable = b_name, issues = character(0),
                 warnings = character(0), confounded = FALSE)

    # -- Column existence --
    if (!b_name %in% colnames(pd)) {
      diag$issues <- c(diag$issues,
                       paste0("Column '", b_name, "' not found in pd."))
      diag$confounded <- TRUE
      details_list[[b_name]] <- diag
      confounded_list <- c(confounded_list, b_name)
      message("[mETHYLotest Confounding] [SKIP] '", b_name,
              "': column not found.")
      next
    }

    pd[[b_name]] <- droplevels(as.factor(pd[[b_name]]))
    pd_clean <- pd[!is.na(pd[[variablename]]) & !is.na(pd[[b_name]]), ]

    if (nrow(pd_clean) == 0L) {
      diag$issues <- c(diag$issues, "No non-NA observations.")
      diag$confounded <- TRUE
      details_list[[b_name]] <- diag
      confounded_list <- c(confounded_list, b_name)
      message("[mETHYLotest Confounding] [SKIP] '", b_name,
              "': no valid observations.")
      next
    }

    is_confounded <- FALSE
    has_warnings  <- FALSE
    tbl <- table(Bio = pd_clean[[variablename]],
                 Batch = pd_clean[[b_name]])
    diag$contingency <- tbl
    diag$n_levels    <- ncol(tbl)

    # -- CHECK 0: Single level --
    if (ncol(tbl) < 2L) {
      msg <- paste0("'", b_name, "' has only ", ncol(tbl),
                    " level(s). Nothing to correct.")
      diag$issues <- c(diag$issues, msg)
      is_confounded <- TRUE
      message("[mETHYLotest Confounding] [FAIL] ", msg)
    }

    # -- CHECK 1: Perfect aliasing --
    # A bio group exists in only 1 batch level
    # This means ComBat CANNOT separate batch from biology for that group
    if (!is_confounded) {
      batches_per_bio <- rowSums(tbl > 0)
      nested_groups   <- names(batches_per_bio[batches_per_bio == 1L])

      if (length(nested_groups) > 0L) {
        msg <- paste0("'", b_name, "': group(s) ",
                      paste(nested_groups, collapse = ", "),
                      " found in only 1 batch level.",
                      " ComBat cannot separate batch from biology.")
        diag$issues <- c(diag$issues, msg)
        is_confounded <- TRUE
        message("[mETHYLotest Confounding] [FAIL] ", msg)
      }
    }

    # -- CHECK 2: Collinearity --
    # Some batch levels contain only 1 bio group
    # CRITICAL: only a problem if ALL levels are pure (total confounding)
    # If some levels are mixed, ComBat can still estimate the effect
    if (!is_confounded) {
      bio_per_batch    <- colSums(tbl > 0)
      pure_batches     <- names(bio_per_batch[bio_per_batch == 1L])
      mixed_batches    <- names(bio_per_batch[bio_per_batch > 1L])
      frac_pure        <- length(pure_batches) / ncol(tbl)

      if (length(pure_batches) > 0L) {
        if (length(mixed_batches) == 0L) {
          # ALL levels are pure = total confounding
          msg <- paste0("'", b_name,
                        "': ALL batch levels contain only 1 bio group.",
                        " Total confounding. ComBat impossible.")
          diag$issues <- c(diag$issues, msg)
          is_confounded <- TRUE
          message("[mETHYLotest Confounding] [FAIL] ", msg)
        } else {
          # Partial: some pure, some mixed = WARNING only
          msg <- paste0("'", b_name, "': ", length(pure_batches),
                        "/", ncol(tbl), " batch level(s) (",
                        paste(pure_batches, collapse = ", "),
                        ") contain only 1 bio group.",
                        " Mixed levels available: ",
                        paste(mixed_batches, collapse = ", "),
                        ". ComBat can proceed but results should be",
                        " interpreted with caution.")
          diag$warnings <- c(diag$warnings, msg)
          has_warnings <- TRUE
          message("[mETHYLotest Confounding] [WARN] ", msg)
        }
      }
    }

    # -- CHECK 3: Singleton batches --
    if (!is_confounded) {
      samples_per_batch <- colSums(tbl)
      singletons <- names(samples_per_batch[samples_per_batch < 2L])

      if (length(singletons) > 0L) {
        n_singletons <- length(singletons)
        n_total      <- ncol(tbl)

        if (n_singletons == n_total) {
          # ALL batches are singletons
          msg <- paste0("'", b_name,
                        "': ALL batch levels have only 1 sample.",
                        " ComBat impossible.")
          diag$issues <- c(diag$issues, msg)
          is_confounded <- TRUE
          message("[mETHYLotest Confounding] [FAIL] ", msg)
        } else {
          # Some singletons
          msg <- paste0("'", b_name, "': ", n_singletons,
                        "/", n_total, " batch level(s) (",
                        paste(singletons, collapse = ", "),
                        ") have only 1 sample.",
                        " These samples may have unstable corrections.")
          diag$warnings <- c(diag$warnings, msg)
          has_warnings <- TRUE
          message("[mETHYLotest Confounding] [WARN] ", msg)
        }
      }
    }

    diag$confounded <- is_confounded

    if (is_confounded) {
      confounded_list <- c(confounded_list, b_name)
      message("[mETHYLotest Confounding]   Contingency table for '",
              b_name, "':")
      tbl_str <- utils::capture.output(print(tbl))
      for (line in tbl_str)
        message("[mETHYLotest Confounding]     ", line)
    } else if (has_warnings) {
      clean_list <- c(clean_list, b_name)
      warn_list  <- c(warn_list, b_name)
      message("[mETHYLotest Confounding] [PASS with WARNINGS] '",
              b_name, "' (", ncol(tbl), " levels)")
      message("[mETHYLotest Confounding]   Contingency table:")
      tbl_str <- utils::capture.output(print(tbl))
      for (line in tbl_str)
        message("[mETHYLotest Confounding]     ", line)
    } else {
      clean_list <- c(clean_list, b_name)
      message("[mETHYLotest Confounding] [PASS] '", b_name,
              "' (", ncol(tbl), " levels, no issues)")
    }

    details_list[[b_name]] <- diag
  }

  message(strrep("-", 60))

  if (length(confounded_list) > 0L && length(clean_list) > 0L) {
    message("[mETHYLotest Confounding] PARTIAL CONFOUNDING:")
    message("[mETHYLotest Confounding]   Excluded: ",
            paste(confounded_list, collapse = ", "))
    message("[mETHYLotest Confounding]   Proceeding: ",
            paste(clean_list, collapse = ", "))
  } else if (length(confounded_list) > 0L) {
    message("[mETHYLotest Confounding] ALL CONFOUNDED.",
            " No variables available.")
  } else {
    message("[mETHYLotest Confounding] All variables passed.")
  }

  if (length(warn_list) > 0L)
    message("[mETHYLotest Confounding] Variables with warnings: ",
            paste(warn_list, collapse = ", "),
            " (proceed with caution)")

  message(strrep("=", 60))

  list(
    confounded = unique(confounded_list),
    clean      = unique(clean_list),
    warnings   = unique(warn_list),
    details    = details_list
  )
}
