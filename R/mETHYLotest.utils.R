#' Format bytes into a human-readable string
#' @param bytes A number (size in bytes)
#' @return A formatted string (e.g., "10.5 Mb")
#' @noRd
mf_format_bytes <- function(bytes) {
  if (is.na(bytes) || bytes == 0) return("0 bytes")
  units <- c("bytes", "Kb", "Mb", "Gb", "Tb", "Pb")
  i <- floor(log(bytes, 1024))
  size <- bytes / (1024 ^ i)
  paste(sprintf("%.1f", size), units[i + 1])
}
# ==============================================================================
# INTERNAL HELPERS (non-exported)
# ==============================================================================

#' @keywords internal
.epicv2_strip_suffixes <- function(load_v2, duplicate_strategy) {

  raw_names   <- rownames(load_v2$beta)
  clean_names <- sub("_.*", "", raw_names)
  n_stripped  <- sum(clean_names != raw_names)

  message("[mETHYLotest] EPICv2 probes with suffix stripped: ",
          n_stripped, " / ", length(raw_names))

  n_dups <- sum(duplicated(clean_names))

  if (n_dups > 0) {
    warning(
      "[mETHYLotest] ", n_dups, " duplicate EPICv2 probe ID(s) after suffix",
      " removal. Collapsing with strategy = '", duplicate_strategy, "'."
    )
    collapse <- .make_probe_collapser(duplicate_strategy, clean_names)
    for (slot in c("beta", "M", "intensity", "detP", "beadcount")) {
      if (!is.null(load_v2[[slot]]))
        load_v2[[slot]] <- collapse(load_v2[[slot]])
    }
  } else {
    for (slot in c("beta", "M", "intensity", "detP", "beadcount")) {
      if (!is.null(load_v2[[slot]]))
        rownames(load_v2[[slot]]) <- clean_names
    }
  }

  load_v2
}

#' @keywords internal
.make_probe_collapser <- function(strategy, clean_names) {

  if (strategy == "first") {
    keep        <- !duplicated(clean_names)
    dedup_names <- clean_names[keep]
    function(mat) {
      out <- mat[keep, , drop = FALSE]
      rownames(out) <- dedup_names
      out
    }

  } else {  # "mean"
    idx_by_probe <- split(seq_along(clean_names), clean_names)
    unique_names <- names(idx_by_probe)
    function(mat) {
      out <- matrix(
        NA_real_,
        nrow     = length(idx_by_probe),
        ncol     = ncol(mat),
        dimnames = list(unique_names, colnames(mat))
      )
      for (g in unique_names) {
        idx      <- idx_by_probe[[g]]
        out[g, ] <- if (length(idx) == 1L) mat[idx, ]
        else colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      }
      out
    }
  }
}


# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

#' Harmonize and Merge Illumina Methylation Array Datasets
#'
#' @description
#' Merges two or three ChAMP objects from different Illumina array generations
#' (450K, EPICv1, EPICv2) into a single harmonized object. Only probes present
#' in **all** provided arrays are retained (intersection).
#'
#' EPICv2 replicate probe suffixes (e.g. \code{_BC21}, \code{_TC21}) are
#' automatically stripped to recover EPICv1/450K-compatible probe IDs.
#' Resulting duplicate probes are collapsed according to
#' \code{duplicate_strategy}.
#'
#' @details
#' \strong{Why suffix stripping is safe for EPICv2:} Standard Illumina probe
#' identifiers (\code{cg}, \code{rs}, \code{ch} prefixes) never contain
#' underscores in their base ID. Only EPICv2 replicate suffixes introduce
#' underscores, so \code{sub("_.*", "")} unambiguously recovers the original
#' probe ID.
#'
#' \strong{Duplicate strategy (EPICv2 only):}
#' \itemize{
#'   \item \code{"mean"} (default): \code{colMeans} across replicate probes
#'     with \code{na.rm = TRUE}. Recommended — uses all available signal.
#'   \item \code{"first"}: keeps the first occurrence. Fast but arbitrary.
#' }
#'
#' \strong{Supported combinations:}
#' \itemize{
#'   \item EPICv1 + EPICv2
#'   \item EPICv1 + 450K
#'   \item EPICv2 + 450K
#'   \item EPICv1 + EPICv2 + 450K
#' }
#'
#' @param loads A \strong{named list} of ChAMP objects. Accepted names:
#'   \code{"EPICv1"}, \code{"EPICv2"}, \code{"450K"}.
#'   At least two elements are required.
#' @param duplicate_strategy Character. How to collapse EPICv2 replicate probes
#'   after suffix removal. One of \code{"mean"} (default) or \code{"first"}.
#'
#' @return A merged ChAMP-like list with \code{$beta}, \code{$pd}, and
#'   optionally \code{$M}, \code{$intensity}, \code{$detP}, \code{$beadcount}.
#'   Slots absent from any one array are silently dropped from the merged output.
#'
#' @seealso \code{\link{mETHYLotest.utils.HarmonizeEpicV2toV1}} (compatibility
#'   wrapper), \code{\link{mETHYLotest.EPIC.pipeline}}
#'
#' @export
mETHYLotest.utils.HarmonizeArrays <- function(
    loads,
    duplicate_strategy = c("mean", "first")
) {

  duplicate_strategy <- match.arg(duplicate_strategy)

  # ── 1. Input validation ──────────────────────────────────────────────────
  VALID_NAMES <- c("EPICv1", "EPICv2", "450K")

  if (!is.list(loads) || is.null(names(loads)))
    stop("'loads' must be a named list. ",
         "Valid names: ", paste(VALID_NAMES, collapse = ", "), ".")

  unknown <- setdiff(names(loads), VALID_NAMES)
  if (length(unknown) > 0)
    stop("Unknown array type(s): ", paste(unknown, collapse = ", "),
         ". Valid names: ", paste(VALID_NAMES, collapse = ", "), ".")

  if (length(loads) < 2)
    stop("At least 2 arrays are required. ",
         "For a single array no harmonization is needed.")

  for (nm in names(loads)) {
    if (is.null(loads[[nm]]$beta))
      stop("Object '", nm, "' does not contain a '$beta' matrix.")
  }

  # Duplicate sample names across arrays → corrupted cbind
  all_samples <- unlist(lapply(loads, function(x) colnames(x$beta)))
  dup_samples <- unique(all_samples[duplicated(all_samples)])
  if (length(dup_samples) > 0)
    stop(
      length(dup_samples), " duplicate sample name(s) across arrays: ",
      paste(dup_samples, collapse = ", "),
      ". Sample names must be unique across all arrays."
    )

  message("[mETHYLotest] Starting harmonization: ",
          paste(names(loads), collapse = " + "))

  # ── 2. EPICv2 suffix stripping ───────────────────────────────────────────
  if ("EPICv2" %in% names(loads)) {
    loads[["EPICv2"]] <- .epicv2_strip_suffixes(
      loads[["EPICv2"]], duplicate_strategy
    )
  }

  # ── 3. Common probes ─────────────────────────────────────────────────────
  probe_sets    <- lapply(loads, function(x) rownames(x$beta))
  probes_in_com <- Reduce(intersect, probe_sets)

  for (nm in names(loads))
    message("[mETHYLotest] ", nm, " probes: ", nrow(loads[[nm]]$beta))

  message("[mETHYLotest] Common probes retained: ", length(probes_in_com))

  for (nm in names(loads)) {
    dropped <- nrow(loads[[nm]]$beta) - length(probes_in_com)
    message("[mETHYLotest]   ", nm, "-only probes dropped: ", dropped)
  }

  if (length(probes_in_com) == 0)
    stop("No common probes found across: ",
         paste(names(loads), collapse = ", "), ".")

  # ── 4. Merge phenodata ───────────────────────────────────────────────────
  all_pd_cols    <- unique(unlist(lapply(loads, function(x) colnames(x$pd))))
  common_pd_cols <- Reduce(intersect, lapply(loads, function(x) colnames(x$pd)))
  dropped_cols   <- setdiff(all_pd_cols, common_pd_cols)

  if (length(dropped_cols) > 0)
    warning(
      "pd column(s) not present in all arrays (dropped): ",
      paste(dropped_cols, collapse = ", ")
    )

  merged_pd <- do.call(
    rbind,
    lapply(loads, function(x) x$pd[, common_pd_cols, drop = FALSE])
  )

  # ── 5. Merge matrices ────────────────────────────────────────────────────
  # Only merge a slot if ALL arrays have it; warn and skip otherwise.
  merge_slot <- function(slot_name) {
    present <- vapply(loads, function(x) !is.null(x[[slot_name]]), logical(1))
    if (all(present)) {
      do.call(cbind,
              lapply(loads, function(x) x[[slot_name]][probes_in_com, , drop = FALSE]))
    } else {
      if (any(present))
        message("[mETHYLotest] '$", slot_name,
                "' missing in: ",
                paste(names(loads)[!present], collapse = ", "),
                " — slot skipped.")
      NULL
    }
  }

  myLoad <- list(
    pd        = merged_pd,
    beta      = merge_slot("beta"),
    M         = merge_slot("M"),
    intensity = merge_slot("intensity"),
    detP      = merge_slot("detP"),
    beadcount = merge_slot("beadcount")
  )

  # Drop NULL optional slots
  myLoad <- myLoad[!vapply(myLoad, is.null, logical(1))]

  message("[mETHYLotest] Harmonization complete.")
  message("[mETHYLotest] Final dimensions: ",
          nrow(myLoad$beta), " probes x ", ncol(myLoad$beta), " samples.")

  invisible(myLoad)
}


# ==============================================================================
# BACKWARD-COMPATIBILITY WRAPPER
# ==============================================================================

#' Harmonize EPICv1 and EPICv2 (compatibility wrapper)
#'
#' @description
#' Wrapper around \code{\link{mETHYLotest.utils.HarmonizeArrays}} for the
#' EPICv1 + EPICv2 two-array case. Kept for backward compatibility.
#'
#' @param myLoad_v1 ChAMP object for EPICv1 samples.
#' @param myLoad_v2 ChAMP object for EPICv2 samples.
#' @param duplicate_strategy Passed to
#'   \code{\link{mETHYLotest.utils.HarmonizeArrays}}.
#'
#' @seealso \code{\link{mETHYLotest.utils.HarmonizeArrays}}
#' @export
mETHYLotest.utils.HarmonizeEpicV2toV1 <- function(
    myLoad_v1,
    myLoad_v2,
    duplicate_strategy = c("mean", "first")
) {
  .Deprecated("mETHYLotest.utils.HarmonizeArrays",
              msg = paste(
                "'mETHYLotest.utils.HarmonizeEpicV2toV1' is deprecated.",
                "Use 'mETHYLotest.utils.HarmonizeArrays' instead."
              ))
  mETHYLotest.utils.HarmonizeArrays(
    loads              = list(EPICv1 = myLoad_v1, EPICv2 = myLoad_v2),
    duplicate_strategy = match.arg(duplicate_strategy)
  )
}
