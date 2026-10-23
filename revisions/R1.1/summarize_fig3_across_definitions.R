#!/usr/bin/env Rscript

# Audit-only cross-definition summary for the six R1.1 rapid Fig. 3 reruns.
#
# This script deliberately does not fit a model, rerun a test, choose a new
# threshold, or make a figure.  It reads only finalized D1--D6 exports and
# writes descriptive counts/directions from the quantities already used in
# restored Fig. 3.  Output replacement is atomic: if any definition, table,
# schema, or developmental-evidence guard fails, no output is written.

options(stringsAsFactors = FALSE, warn = 1)

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

script_file <- function() {
  candidates <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(candidates) == 0L) {
    return(normalizePath(".", mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", candidates[[1L]]), mustWork = TRUE)
}

SCRIPT_PATH <- script_file()
SCRIPT_DIR <- if (dir.exists(SCRIPT_PATH)) SCRIPT_PATH else dirname(SCRIPT_PATH)

usage <- function() {
  paste(
    "Usage:",
    "  Rscript summarize_fig3_across_definitions.R",
    paste(
      "    [--input-root PATH] [--output PATH] [--consistency-output PATH]",
      "[--focus-output PATH]"
    ),
    "",
    "Defaults:",
    paste0("  --input-root ", file.path(SCRIPT_DIR, "full_original_layout")),
    paste0(
      "  --output ",
      file.path(SCRIPT_DIR, "full_original_layout", "Fig3_across_definitions_audit.csv")
    ),
    paste0(
      "  --consistency-output ",
      file.path(
        SCRIPT_DIR,
        "full_original_layout",
        "Fig3_across_definitions_direction_consistency.csv"
      )
    ),
    paste0(
      "  --focus-output ",
      file.path(
        SCRIPT_DIR,
        "full_original_layout",
        "Fig3_focus_genera_HDR_minus_LDR.csv"
      )
    ),
    sep = "\n"
  )
}

parse_args <- function(args) {
  parsed <- list(
    input_root = file.path(SCRIPT_DIR, "full_original_layout"),
    output = file.path(
      SCRIPT_DIR,
      "full_original_layout",
      "Fig3_across_definitions_audit.csv"
    ),
    consistency_output = file.path(
      SCRIPT_DIR,
      "full_original_layout",
      "Fig3_across_definitions_direction_consistency.csv"
    ),
    focus_output = file.path(
      SCRIPT_DIR,
      "full_original_layout",
      "Fig3_focus_genera_HDR_minus_LDR.csv"
    )
  )
  if (any(args %in% c("-h", "--help"))) {
    cat(usage(), "\n")
    quit(status = 0L)
  }
  valid <- c(
    "--input-root" = "input_root",
    "--output" = "output",
    "--consistency-output" = "consistency_output",
    "--focus-output" = "focus_output"
  )
  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (!token %in% names(valid)) {
      fail("Unknown argument: ", token, "\n", usage())
    }
    if (index == length(args)) {
      fail("Missing value after ", token, "\n", usage())
    }
    parsed[[valid[[token]]]] <- args[[index + 1L]]
    index <- index + 2L
  }
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_root <- normalizePath(args$input_root, mustWork = FALSE)
output_path <- normalizePath(args$output, mustWork = FALSE)
consistency_path <- normalizePath(args$consistency_output, mustWork = FALSE)
focus_path <- normalizePath(args$focus_output, mustWork = FALSE)

output_targets <- c(output_path, consistency_path, focus_path)
if (anyDuplicated(output_targets)) {
  fail("--output, --consistency-output, and --focus-output must be different files")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  fail("Required R package 'jsonlite' is unavailable; no output was written")
}

definition_slugs <- c(
  D1 = "D1_endpoint_abs_shannon",
  D2 = "D2_endpoint_bray",
  D3 = "D3_endpoint_jaccard",
  D4 = "D4_auc_abs_shannon",
  D5 = "D5_auc_bray",
  D6 = "D6_auc_jaccard"
)
EXPECTED_FIXED_TAXA <- 259L
FOCUS_GENERA <- c("Faecalibacterium", "Bacteroides", "Alistipes")

required_table_schemas <- list(
  "perturbation_effects_combined.csv" = c(
    "taxon", "LDR_median", "HDR_median", "HDR_minus_LDR_median"
  ),
  "growth_effects_combined.csv" = c(
    "taxon", "LDR_median", "HDR_median", "HDR_minus_LDR_median"
  ),
  "posterior_diagnostics.csv" = c(
    "group", "parameter_type", "parameter", "rhat_rank_folded", "n_chains",
    "draws_per_chain_used", "diagnostic_status", "evidence_tier"
  ),
  "LDR_module_assignments.csv" = c("taxon", "module"),
  "HDR_module_assignments.csv" = c("taxon", "module"),
  "Fig3C_top25_perturbation_dumbbell.csv" = c(
    "taxon", "LDR_median", "HDR_median", "contrast", "source_visible_slot"
  ),
  "Fig3D_growth_vs_perturbation.csv" = c(
    "taxon", "group", "perturbation_response", "growth_rate", "genus"
  ),
  "Fig3E_HDR_module_spearman_edges_drawn.csv" = c(
    "source", "target", "rho", "group", "observation_unit"
  ),
  "Fig3F_LDR_module_spearman_edges_drawn.csv" = c(
    "source", "target", "rho", "group", "observation_unit"
  )
)

as_scalar_number <- function(x, label) {
  value <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  if (length(value) != 1L || is.na(value)) {
    fail(label, " must contain exactly one numeric value")
  }
  value
}

as_scalar_character <- function(x, label) {
  value <- as.character(unlist(x, use.names = FALSE))
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    fail(label, " must contain exactly one non-empty character value")
  }
  value
}

read_json_object <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    fail(label, " is missing or empty: ", path)
  }
  result <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(error) fail("Cannot parse ", label, " at ", path, ": ", error$message)
  )
  if (!is.list(result)) {
    fail(label, " is not a JSON object: ", path)
  }
  result
}

read_required_csv <- function(path, columns, allow_empty = FALSE, sep = ",") {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    fail("Required table is missing or empty: ", path)
  }
  frame <- tryCatch(
    read.csv(
      path,
      sep = sep,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA")
    ),
    error = function(error) fail("Cannot read table ", path, ": ", error$message)
  )
  missing_columns <- setdiff(columns, names(frame))
  if (length(missing_columns) > 0L) {
    fail(
      "Table ", path, " lacks required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (!allow_empty && nrow(frame) == 0L) {
    fail("Required table has no data rows: ", path)
  }
  frame
}

require_unique_nonmissing <- function(frame, column, label) {
  values <- as.character(frame[[column]])
  if (anyNA(values) || any(!nzchar(values)) || anyDuplicated(values)) {
    fail(label, " must be non-missing and unique")
  }
  invisible(values)
}

require_finite <- function(frame, columns, label) {
  for (column in columns) {
    values <- suppressWarnings(as.numeric(frame[[column]]))
    if (length(values) != nrow(frame) || any(!is.finite(values))) {
      fail(label, "$", column, " must be numeric and finite in every row")
    }
  }
  invisible(TRUE)
}

# Exact R mirror of render_fig3.py::_pretty_taxon and ::_genus:
#   1. remove one leading "s__" prefix;
#   2. replace every underscore with a space;
#   3. remove literal square brackets;
#   4. take the first whitespace-delimited token.
# The finalized Fig3D export already contains a genus column.  Recomputing it
# here is an integrity check that prevents a later, broader taxonomic matching
# rule from silently entering the reviewer summary.
render_fig3_pretty_taxon <- function(taxa) {
  values <- sub("^s__", "", as.character(taxa))
  values <- gsub("_", " ", values, fixed = TRUE)
  values <- gsub("[", "", values, fixed = TRUE)
  gsub("]", "", values, fixed = TRUE)
}

render_fig3_genus <- function(taxa) {
  pretty <- trimws(render_fig3_pretty_taxon(taxa))
  vapply(
    strsplit(pretty, "[[:space:]]+"),
    function(parts) {
      if (length(parts) == 0L || !nzchar(parts[[1L]])) "" else parts[[1L]]
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

FOCUS_MATCHING_RULE <- paste(
  "Exact render_fig3.py _pretty_taxon/_genus rule:",
  "remove one leading s__, replace underscores with spaces, remove [ and ],",
  "then use the first whitespace-delimited token; exact genus equality only."
)

validate_definition <- function(definition_id, slug) {
  root <- file.path(input_root, slug)
  marker_path <- file.path(root, "mdsine2", "finalized.json")
  marker <- read_json_object(marker_path, paste0(definition_id, " finalized marker"))
  marker_id <- as_scalar_character(marker$definition_id, paste0(definition_id, " marker definition_id"))
  if (!identical(marker_id, definition_id)) {
    fail(definition_id, ": finalized marker identifies ", marker_id)
  }
  profile <- marker$profile
  if (!is.list(profile)) {
    fail(definition_id, ": finalized marker lacks profile object")
  }
  if (as_scalar_number(profile$chains_per_group, paste0(definition_id, " chains_per_group")) != 1) {
    fail(definition_id, ": expected exactly one rapid chain per group")
  }
  retained <- as_scalar_number(profile$retained_draws, paste0(definition_id, " retained_draws"))
  if (retained != 100) {
    fail(definition_id, ": expected 100 retained rapid-chain draws, observed ", retained)
  }
  tier <- as_scalar_character(profile$evidence_tier, paste0(definition_id, " evidence_tier"))
  if (!identical(tier, "pilot_single_or_short_chain")) {
    fail(definition_id, ": unexpected evidence tier ", tier)
  }
  boundary <- as_scalar_character(
    marker$interpretation_boundary,
    paste0(definition_id, " interpretation_boundary")
  )

  for (suffix in c("pdf", "tif")) {
    figure_path <- file.path(root, "figures", paste0("Fig3_restored.", suffix))
    if (!file.exists(figure_path) || file.info(figure_path)$size <= 0) {
      fail(definition_id, ": finalized figure is missing or empty: ", figure_path)
    }
  }

  qa_path <- file.path(root, "qa", "validation_checks.tsv")
  qa <- read_required_csv(
    qa_path,
    c("name", "status", "severity"),
    allow_empty = FALSE,
    sep = "\t"
  )
  error_checks <- qa[tolower(qa$severity) == "error", , drop = FALSE]
  if (nrow(error_checks) == 0L || any(toupper(error_checks$status) != "PASS")) {
    failed <- error_checks$name[toupper(error_checks$status) != "PASS"]
    fail(
      definition_id,
      ": strict Fig3 QA error-level checks did not all pass",
      if (length(failed)) paste0(": ", paste(failed, collapse = "; ")) else ""
    )
  }

  tables_dir <- file.path(root, "tables")
  tables <- lapply(names(required_table_schemas), function(filename) {
    read_required_csv(
      file.path(tables_dir, filename),
      required_table_schemas[[filename]],
      allow_empty = grepl("edges_drawn", filename, fixed = TRUE)
    )
  })
  names(tables) <- names(required_table_schemas)

  manifest_path <- file.path(tables_dir, "posterior_summary_manifest.json")
  manifest <- read_json_object(
    manifest_path,
    paste0(definition_id, " posterior summary manifest")
  )
  bf_threshold <- as_scalar_number(
    manifest$bayes_factor_threshold,
    paste0(definition_id, " BF threshold")
  )
  if (bf_threshold != 10) {
    fail(definition_id, ": expected the existing BF threshold 10, observed ", bf_threshold)
  }
  for (group in c("LDR", "HDR")) {
    group_manifest <- manifest$groups[[group]]
    if (!is.list(group_manifest)) {
      fail(definition_id, ": posterior manifest lacks group ", group)
    }
    group_tier <- as_scalar_character(
      group_manifest$evidence_tier,
      paste0(definition_id, " ", group, " manifest evidence_tier")
    )
    if (!identical(group_tier, tier)) {
      fail(definition_id, ": ", group, " manifest evidence tier disagrees with marker")
    }
    chains <- unlist(group_manifest$chains, use.names = FALSE)
    draws <- suppressWarnings(as.numeric(unlist(
      group_manifest$posterior_draws_by_chain,
      use.names = FALSE
    )))
    if (length(chains) != 1L || length(draws) != 1L || is.na(draws) || draws != retained) {
      fail(definition_id, ": ", group, " is not the expected one-chain/100-draw summary")
    }
  }

  diagnostics <- tables[["posterior_diagnostics.csv"]]
  observed_groups <- sort(unique(as.character(diagnostics$group)))
  if (!identical(observed_groups, c("HDR", "LDR"))) {
    fail(definition_id, ": diagnostics must contain exactly HDR and LDR")
  }
  if (any(suppressWarnings(as.numeric(diagnostics$n_chains)) != 1)) {
    fail(definition_id, ": diagnostics contain a non-single-chain row")
  }
  if (any(suppressWarnings(as.numeric(diagnostics$draws_per_chain_used)) != retained)) {
    fail(definition_id, ": diagnostics contain unexpected retained-draw counts")
  }
  if (any(as.character(diagnostics$evidence_tier) != tier)) {
    fail(definition_id, ": diagnostics contain an unexpected evidence tier")
  }
  if (any(as.character(diagnostics$diagnostic_status) != "rhat_unavailable_single_chain")) {
    fail(definition_id, ": diagnostics do not uniformly flag single-chain R-hat as unavailable")
  }
  rhat <- suppressWarnings(as.numeric(diagnostics$rhat_rank_folded))
  if (any(is.finite(rhat))) {
    fail(definition_id, ": finite R-hat values are not allowed for this one-chain audit")
  }

  perturbation <- tables[["perturbation_effects_combined.csv"]]
  growth <- tables[["growth_effects_combined.csv"]]
  require_unique_nonmissing(perturbation, "taxon", paste0(definition_id, " perturbation taxa"))
  require_unique_nonmissing(growth, "taxon", paste0(definition_id, " growth taxa"))
  require_finite(
    perturbation,
    c("LDR_median", "HDR_median", "HDR_minus_LDR_median"),
    paste0(definition_id, " perturbation")
  )
  require_finite(
    growth,
    c("LDR_median", "HDR_median", "HDR_minus_LDR_median"),
    paste0(definition_id, " growth")
  )
  if (!identical(sort(perturbation$taxon), sort(growth$taxon))) {
    fail(definition_id, ": growth and perturbation taxon universes differ")
  }
  if (nrow(perturbation) != EXPECTED_FIXED_TAXA) {
    fail(
      definition_id,
      ": expected the fixed ",
      EXPECTED_FIXED_TAXA,
      "-taxon Fig3 universe, observed ",
      nrow(perturbation)
    )
  }

  for (group in c("LDR", "HDR")) {
    assignments <- tables[[paste0(group, "_module_assignments.csv")]]
    require_unique_nonmissing(
      assignments,
      "taxon",
      paste0(definition_id, " ", group, " module taxa")
    )
    if (!identical(sort(assignments$taxon), sort(perturbation$taxon))) {
      fail(definition_id, ": ", group, " module taxon universe differs from effects")
    }
    if (anyNA(assignments$module) || any(!nzchar(as.character(assignments$module)))) {
      fail(definition_id, ": ", group, " contains missing module labels")
    }
  }

  top25 <- tables[["Fig3C_top25_perturbation_dumbbell.csv"]]
  require_unique_nonmissing(top25, "taxon", paste0(definition_id, " Fig3C taxa"))
  require_finite(
    top25,
    c("LDR_median", "HDR_median", "contrast", "source_visible_slot"),
    paste0(definition_id, " Fig3C")
  )
  slots <- sort(as.integer(top25$source_visible_slot))
  if (nrow(top25) != 25L || !identical(slots, 1:25)) {
    fail(definition_id, ": Fig3C must contain the exact 25 visible slots")
  }
  if (!all(top25$taxon %in% perturbation$taxon)) {
    fail(definition_id, ": Fig3C contains taxa outside the fixed effect universe")
  }

  fig3d <- tables[["Fig3D_growth_vs_perturbation.csv"]]
  require_finite(
    fig3d,
    c("perturbation_response", "growth_rate"),
    paste0(definition_id, " Fig3D")
  )
  if (!identical(sort(unique(as.character(fig3d$group))), c("HDR", "LDR"))) {
    fail(definition_id, ": Fig3D must contain exactly HDR and LDR")
  }
  expected_genus <- render_fig3_genus(fig3d$taxon)
  exported_genus <- as.character(fig3d$genus)
  if (anyNA(exported_genus) || !identical(exported_genus, expected_genus)) {
    mismatch <- which(is.na(exported_genus) | exported_genus != expected_genus)
    preview <- paste(head(as.character(fig3d$taxon[mismatch]), 5L), collapse = ", ")
    fail(
      definition_id,
      ": Fig3D genus labels do not exactly reproduce render_fig3.py for: ",
      preview
    )
  }
  for (focus_genus in FOCUS_GENERA) {
    focus <- fig3d[exported_genus == focus_genus, , drop = FALSE]
    counts <- table(factor(focus$group, levels = c("LDR", "HDR")))
    if (any(counts <= 0L) || counts[["LDR"]] != counts[["HDR"]]) {
      fail(
        definition_id,
        ": incomplete paired Fig3D points for focus genus ",
        focus_genus
      )
    }
    ldr_taxa <- sort(as.character(focus$taxon[focus$group == "LDR"]))
    hdr_taxa <- sort(as.character(focus$taxon[focus$group == "HDR"]))
    if (!identical(ldr_taxa, hdr_taxa)) {
      fail(
        definition_id,
        ": HDR/LDR taxon identities differ for focus genus ",
        focus_genus
      )
    }
  }
  expected_points <- 2L * nrow(perturbation)
  if (nrow(fig3d) != expected_points || anyDuplicated(fig3d[c("taxon", "group")])) {
    fail(
      definition_id,
      ": Fig3D must contain one point per taxon and group (expected ",
      expected_points,
      ")"
    )
  }

  network_files <- c(
    HDR = "Fig3E_HDR_module_spearman_edges_drawn.csv",
    LDR = "Fig3F_LDR_module_spearman_edges_drawn.csv"
  )
  for (group in names(network_files)) {
    network <- tables[[network_files[[group]]]]
    if (nrow(network) > 0L) {
      require_finite(network, "rho", paste0(definition_id, " ", group, " drawn network"))
      if (any(abs(as.numeric(network$rho)) > 1)) {
        fail(definition_id, ": ", group, " drawn network contains |rho| > 1")
      }
      if (any(as.character(network$group) != group)) {
        fail(definition_id, ": ", group, " drawn network has wrong group labels")
      }
      if (any(as.character(network$observation_unit) != "canonical_stage")) {
        fail(definition_id, ": ", group, " network does not use canonical stages")
      }
      edge_keys <- apply(
        cbind(
          pmin(as.character(network$source), as.character(network$target)),
          pmax(as.character(network$source), as.character(network$target))
        ),
        1L,
        paste,
        collapse = "--"
      )
      if (anyDuplicated(edge_keys)) {
        fail(definition_id, ": ", group, " drawn network has duplicate undirected edges")
      }
    }
  }

  list(
    definition_id = definition_id,
    slug = slug,
    root = root,
    tables_dir = tables_dir,
    marker = marker,
    evidence_tier = tier,
    retained_draws = retained,
    interpretation_boundary = boundary,
    bf_threshold = bf_threshold,
    tables = tables
  )
}

if (!dir.exists(input_root)) {
  fail("Input root does not exist: ", input_root, "; no output was written")
}

# Complete the gate for all six definitions before constructing or writing any
# result.  This ensures a missing late definition cannot leave a partial CSV.
validated <- Map(
  validate_definition,
  names(definition_slugs),
  unname(definition_slugs)
)
names(validated) <- names(definition_slugs)

reference_taxa <- sort(validated[["D1"]]$tables[["perturbation_effects_combined.csv"]]$taxon)
for (definition_id in names(validated)[-1L]) {
  observed_taxa <- sort(
    validated[[definition_id]]$tables[["perturbation_effects_combined.csv"]]$taxon
  )
  if (!identical(observed_taxa, reference_taxa)) {
    fail(
      definition_id,
      ": fixed Fig3 taxon universe differs from D1; no output was written"
    )
  }
}

evidence_boundary <- paste(
  "Pilot/developmental one-chain-per-group sensitivity output",
  "(140 total iterations, 40 burn-in, 100 retained draws);",
  "R-hat is unavailable and no formal multichain convergence is claimed.",
  "No measured qPCR: model parameters remain on a relative-abundance scale",
  "and must not be interpreted as absolute biomass rates."
)

direction <- function(value) {
  ifelse(value > 0, "positive", ifelse(value < 0, "negative", "zero"))
}

numeric_summary <- function(
  definition,
  panel,
  estimand,
  group_or_contrast,
  values,
  source_filename,
  source_column,
  selection_rule,
  direction_definition,
  focus_genus = ""
) {
  values <- as.numeric(values)
  if (any(!is.finite(values))) {
    fail(definition$definition_id, ": non-finite values reached audit summary for ", estimand)
  }
  n <- length(values)
  if (n == 0L) {
    statistics <- setNames(
      rep(NA_real_, 5L),
      c("median", "q25", "q75", "minimum", "maximum")
    )
    counts <- c(positive = 0L, negative = 0L, zero = 0L)
    median_direction <- "not_applicable_empty"
  } else {
    statistics <- c(
      median = median(values),
      q25 = unname(quantile(values, 0.25, names = FALSE, type = 7)),
      q75 = unname(quantile(values, 0.75, names = FALSE, type = 7)),
      minimum = min(values),
      maximum = max(values)
    )
    counts <- c(
      positive = sum(values > 0),
      negative = sum(values < 0),
      zero = sum(values == 0)
    )
    median_direction <- direction(statistics[["median"]])
  }
  source_path <- file.path(definition$tables_dir, source_filename)
  data.frame(
    definition_id = definition$definition_id,
    definition_slug = definition$slug,
    figure_panel = panel,
    estimand = estimand,
    group_or_contrast = group_or_contrast,
    focus_genus = focus_genus,
    source_file = normalizePath(source_path, mustWork = TRUE),
    source_md5 = unname(tools::md5sum(source_path)),
    source_column = source_column,
    n_items = n,
    n_positive = unname(counts[["positive"]]),
    n_negative = unname(counts[["negative"]]),
    n_zero = unname(counts[["zero"]]),
    median = unname(statistics[["median"]]),
    q25 = unname(statistics[["q25"]]),
    q75 = unname(statistics[["q75"]]),
    minimum = unname(statistics[["minimum"]]),
    maximum = unname(statistics[["maximum"]]),
    median_direction = median_direction,
    direction_definition = direction_definition,
    selection_rule = selection_rule,
    evidence_tier = definition$evidence_tier,
    chains_per_group = 1L,
    retained_draws_per_chain = definition$retained_draws,
    rhat_status = "unavailable_single_chain",
    qpcr_status = "not_measured_relative_abundance_scale_only",
    interpretation_boundary = evidence_boundary,
    stringsAsFactors = FALSE
  )
}

count_summary <- function(
  definition,
  panel,
  estimand,
  group,
  item_count,
  summary_count,
  source_filename,
  selection_rule
) {
  source_path <- file.path(definition$tables_dir, source_filename)
  data.frame(
    definition_id = definition$definition_id,
    definition_slug = definition$slug,
    figure_panel = panel,
    estimand = estimand,
    group_or_contrast = group,
    focus_genus = "",
    source_file = normalizePath(source_path, mustWork = TRUE),
    source_md5 = unname(tools::md5sum(source_path)),
    source_column = "module",
    n_items = as.integer(item_count),
    n_positive = NA_integer_,
    n_negative = NA_integer_,
    n_zero = NA_integer_,
    median = as.numeric(summary_count),
    q25 = NA_real_,
    q75 = NA_real_,
    minimum = NA_real_,
    maximum = NA_real_,
    median_direction = "not_applicable_count",
    direction_definition = "not applicable; median field stores the number of modules",
    selection_rule = selection_rule,
    evidence_tier = definition$evidence_tier,
    chains_per_group = 1L,
    retained_draws_per_chain = definition$retained_draws,
    rhat_status = "unavailable_single_chain",
    qpcr_status = "not_measured_relative_abundance_scale_only",
    interpretation_boundary = evidence_boundary,
    stringsAsFactors = FALSE
  )
}

definition_rows <- function(definition) {
  tables <- definition$tables
  perturbation <- tables[["perturbation_effects_combined.csv"]]
  growth <- tables[["growth_effects_combined.csv"]]
  top25 <- tables[["Fig3C_top25_perturbation_dumbbell.csv"]]
  fig3d <- tables[["Fig3D_growth_vs_perturbation.csv"]]

  rows <- list()
  add <- function(row) {
    rows[[length(rows) + 1L]] <<- row
  }

  for (group in c("LDR", "HDR")) {
    assignments_name <- paste0(group, "_module_assignments.csv")
    assignments <- tables[[assignments_name]]
    add(count_summary(
      definition,
      if (group == "LDR") "A" else "B",
      "posterior_consensus_module_count",
      group,
      nrow(assignments),
      length(unique(assignments$module)),
      assignments_name,
      "All taxa in the finalized average-linkage consensus assignment"
    ))
  }

  top25_difference <- as.numeric(top25$HDR_median) - as.numeric(top25$LDR_median)
  add(numeric_summary(
    definition,
    "C",
    "displayed_top25_antibiotic_perturbation_HDR_minus_LDR",
    "HDR_minus_LDR",
    top25_difference,
    "Fig3C_top25_perturbation_dumbbell.csv",
    "computed from HDR_median - LDR_median",
    "Exact 25 taxa exported by the finalized Fig3C renderer; no new ranking in R",
    "positive means the displayed posterior median is higher in HDR than LDR"
  ))
  add(numeric_summary(
    definition,
    "C/D",
    "all_taxa_antibiotic_perturbation_HDR_minus_LDR",
    "HDR_minus_LDR",
    perturbation$HDR_minus_LDR_median,
    "perturbation_effects_combined.csv",
    "HDR_minus_LDR_median",
    "All taxa in the fixed finalized Fig3 posterior summary",
    "positive means posterior median perturbation is higher in HDR than LDR"
  ))
  add(numeric_summary(
    definition,
    "D",
    "all_taxa_growth_HDR_minus_LDR",
    "HDR_minus_LDR",
    growth$HDR_minus_LDR_median,
    "growth_effects_combined.csv",
    "HDR_minus_LDR_median",
    "All taxa in the fixed finalized Fig3 posterior summary",
    "positive means posterior median growth is higher in HDR than LDR"
  ))

  for (group in c("LDR", "HDR")) {
    group_points <- fig3d[fig3d$group == group, , drop = FALSE]
    add(numeric_summary(
      definition,
      "D",
      "displayed_antibiotic_perturbation_response",
      group,
      group_points$perturbation_response,
      "Fig3D_growth_vs_perturbation.csv",
      "perturbation_response",
      "All displayed Fig3D points for the stated group",
      "positive/negative follows the model-scale perturbation-response sign"
    ))
    add(numeric_summary(
      definition,
      "D",
      "displayed_growth_rate",
      group,
      group_points$growth_rate,
      "Fig3D_growth_vs_perturbation.csv",
      "growth_rate",
      "All displayed Fig3D points for the stated group",
      "positive/negative follows the model-scale growth-parameter sign"
    ))
  }

  # Fig. 3D highlights these genera using the exact genus string produced by
  # render_fig3.py.  Pair each taxon across HDR/LDR, then describe the vector
  # of taxon-wise HDR-minus-LDR changes.  This is descriptive aggregation of
  # already-exported points, not a hypothesis test or posterior refit.
  for (focus_genus in FOCUS_GENERA) {
    focus <- fig3d[fig3d$genus == focus_genus, , drop = FALSE]
    ldr <- focus[
      focus$group == "LDR",
      c("taxon", "perturbation_response", "growth_rate"),
      drop = FALSE
    ]
    hdr <- focus[
      focus$group == "HDR",
      c("taxon", "perturbation_response", "growth_rate"),
      drop = FALSE
    ]
    paired <- merge(
      ldr,
      hdr,
      by = "taxon",
      suffixes = c("_LDR", "_HDR"),
      sort = TRUE
    )
    if (
      nrow(paired) == 0L ||
        nrow(paired) != nrow(ldr) ||
        nrow(paired) != nrow(hdr) ||
        anyDuplicated(paired$taxon)
    ) {
      fail(
        definition$definition_id,
        ": cannot form a complete one-to-one Fig3D pair for ",
        focus_genus
      )
    }

    perturbation_delta <-
      as.numeric(paired$perturbation_response_HDR) -
      as.numeric(paired$perturbation_response_LDR)
    growth_delta <-
      as.numeric(paired$growth_rate_HDR) -
      as.numeric(paired$growth_rate_LDR)

    add(numeric_summary(
      definition,
      "D",
      "focus_genus_taxonwise_perturbation_HDR_minus_LDR",
      "HDR_minus_LDR",
      perturbation_delta,
      "Fig3D_growth_vs_perturbation.csv",
      "taxon-wise perturbation_response_HDR - perturbation_response_LDR",
      paste(
        FOCUS_MATCHING_RULE,
        "All matched taxa in the stated genus; no ranking or inferential test."
      ),
      paste(
        "median_direction is the sign of the across-taxon median HDR-minus-LDR",
        "model-scale perturbation change; minimum and maximum give its descriptive range"
      ),
      focus_genus = focus_genus
    ))
    add(numeric_summary(
      definition,
      "D",
      "focus_genus_taxonwise_growth_HDR_minus_LDR",
      "HDR_minus_LDR",
      growth_delta,
      "Fig3D_growth_vs_perturbation.csv",
      "taxon-wise growth_rate_HDR - growth_rate_LDR",
      paste(
        FOCUS_MATCHING_RULE,
        "All matched taxa in the stated genus; no ranking or inferential test."
      ),
      paste(
        "median_direction is the sign of the across-taxon median HDR-minus-LDR",
        "model-scale growth change; minimum and maximum give its descriptive range"
      ),
      focus_genus = focus_genus
    ))
  }

  network_specs <- list(
    HDR = list(panel = "E", filename = "Fig3E_HDR_module_spearman_edges_drawn.csv"),
    LDR = list(panel = "F", filename = "Fig3F_LDR_module_spearman_edges_drawn.csv")
  )
  for (group in names(network_specs)) {
    spec <- network_specs[[group]]
    network <- tables[[spec$filename]]
    add(numeric_summary(
      definition,
      spec$panel,
      "drawn_module_stage_trajectory_spearman_rho",
      group,
      network$rho,
      spec$filename,
      "rho",
      paste(
        "Exact undirected edges exported as drawn by the finalized renderer;",
        "R performs no edge selection or hypothesis test"
      ),
      "positive/negative is the sign of the exported Spearman rho"
    ))
  }

  do.call(rbind, rows)
}

audit <- do.call(rbind, lapply(validated, definition_rows))
row.names(audit) <- NULL
audit <- audit[order(
  match(audit$definition_id, names(definition_slugs)),
  audit$figure_panel,
  audit$estimand,
  audit$focus_genus,
  audit$group_or_contrast
), , drop = FALSE]

directional <- audit[
  audit$median_direction != "not_applicable_count",
  ,
  drop = FALSE
]
direction_keys <- unique(directional[c(
  "figure_panel",
  "estimand",
  "group_or_contrast",
  "focus_genus"
)])
direction_keys <- direction_keys[order(
  direction_keys$figure_panel,
  direction_keys$estimand,
  direction_keys$focus_genus,
  direction_keys$group_or_contrast
), , drop = FALSE]

consistency_rows <- lapply(seq_len(nrow(direction_keys)), function(index) {
  key <- direction_keys[index, , drop = FALSE]
  selected <- directional[
    directional$figure_panel == key$figure_panel &
      directional$estimand == key$estimand &
      directional$group_or_contrast == key$group_or_contrast &
      directional$focus_genus == key$focus_genus,
    ,
    drop = FALSE
  ]
  observed <- setNames(selected$median_direction, selected$definition_id)
  observed <- observed[names(definition_slugs)]
  if (length(observed) != 6L || anyNA(observed)) {
    fail("Internal audit error: directional summary does not contain D1--D6 exactly once")
  }
  counts <- table(factor(
    observed,
    levels = c("positive", "negative", "zero", "not_applicable_empty")
  ))
  available <- observed[observed != "not_applicable_empty"]
  nonzero <- unique(available[available != "zero"])
  shared <- if (length(unique(observed)) == 1L) {
    unique(observed)
  } else if (length(available) == 0L) {
    "not_applicable_all_empty"
  } else if (length(unique(available)) == 1L && any(observed == "not_applicable_empty")) {
    paste0(unique(available), "_with_empty_definitions")
  } else if (length(nonzero) == 1L && all(available %in% c(nonzero, "zero"))) {
    paste0(nonzero, "_with_zero_ties")
  } else {
    "mixed"
  }
  data.frame(
    figure_panel = key$figure_panel,
    estimand = key$estimand,
    group_or_contrast = key$group_or_contrast,
    focus_genus = key$focus_genus,
    n_definitions = length(observed),
    n_positive_medians = unname(counts[["positive"]]),
    n_negative_medians = unname(counts[["negative"]]),
    n_zero_medians = unname(counts[["zero"]]),
    n_empty_definitions = unname(counts[["not_applicable_empty"]]),
    shared_descriptive_direction = shared,
    D1_direction = unname(observed[["D1"]]),
    D2_direction = unname(observed[["D2"]]),
    D3_direction = unname(observed[["D3"]]),
    D4_direction = unname(observed[["D4"]]),
    D5_direction = unname(observed[["D5"]]),
    D6_direction = unname(observed[["D6"]]),
    inference_status = paste(
      "descriptive direction audit only; no significance, equivalence, causality,",
      "or convergence claim"
    ),
    evidence_boundary = evidence_boundary,
    stringsAsFactors = FALSE
  )
})
consistency <- do.call(rbind, consistency_rows)
row.names(consistency) <- NULL

focus_audit <- audit[
  audit$focus_genus %in% FOCUS_GENERA &
    audit$estimand %in% c(
      "focus_genus_taxonwise_perturbation_HDR_minus_LDR",
      "focus_genus_taxonwise_growth_HDR_minus_LDR"
    ),
  ,
  drop = FALSE
]
expected_focus_rows <- length(definition_slugs) * length(FOCUS_GENERA) * 2L
focus_keys <- paste(
  focus_audit$definition_id,
  focus_audit$focus_genus,
  focus_audit$estimand,
  sep = "::"
)
if (
  nrow(focus_audit) != expected_focus_rows ||
    anyDuplicated(focus_keys) ||
    any(focus_audit$group_or_contrast != "HDR_minus_LDR")
) {
  fail(
    "Focus-genus audit must contain exactly one perturbation and one growth row ",
    "for each D1--D6 x genus combination; no output was written"
  )
}
focus_export <- data.frame(
  definition_id = focus_audit$definition_id,
  definition_slug = focus_audit$definition_slug,
  genus = focus_audit$focus_genus,
  quantity = ifelse(
    grepl("perturbation", focus_audit$estimand, fixed = TRUE),
    "perturbation_response_HDR_minus_LDR",
    "growth_rate_HDR_minus_LDR"
  ),
  n_paired_taxa = focus_audit$n_items,
  n_positive_taxon_differences = focus_audit$n_positive,
  n_negative_taxon_differences = focus_audit$n_negative,
  n_zero_taxon_differences = focus_audit$n_zero,
  median_HDR_minus_LDR = focus_audit$median,
  median_direction = focus_audit$median_direction,
  q25_HDR_minus_LDR = focus_audit$q25,
  q75_HDR_minus_LDR = focus_audit$q75,
  range_min_HDR_minus_LDR = focus_audit$minimum,
  range_max_HDR_minus_LDR = focus_audit$maximum,
  render_fig3_name_matching_rule = focus_audit$selection_rule,
  source_file = focus_audit$source_file,
  source_md5 = focus_audit$source_md5,
  evidence_tier = focus_audit$evidence_tier,
  chains_per_group = focus_audit$chains_per_group,
  retained_draws_per_chain = focus_audit$retained_draws_per_chain,
  rhat_status = focus_audit$rhat_status,
  qpcr_status = focus_audit$qpcr_status,
  inference_status = paste(
    "descriptive taxon-wise HDR-minus-LDR summary only; no new test,",
    "significance, equivalence, causality, or convergence claim"
  ),
  interpretation_boundary = focus_audit$interpretation_boundary,
  stringsAsFactors = FALSE
)
focus_export <- focus_export[order(
  match(focus_export$definition_id, names(definition_slugs)),
  match(focus_export$genus, FOCUS_GENERA),
  focus_export$quantity
), , drop = FALSE]
row.names(focus_export) <- NULL

atomic_write_csv <- function(frame, path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    fail("Cannot create output directory: ", parent)
  }
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = parent,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary), add = TRUE)
  write.csv(frame, temporary, row.names = FALSE, na = "")
  if (!file.exists(temporary) || file.info(temporary)$size <= 0) {
    fail("Temporary CSV was not written: ", temporary)
  }
  if (!file.rename(temporary, path)) {
    fail("Atomic replacement failed for: ", path)
  }
  invisible(path)
}

# All frames are fully constructed before any output is touched.  The final
# renames are the only mutations performed by this script.
atomic_write_csv(audit, output_path)
atomic_write_csv(consistency, consistency_path)
atomic_write_csv(focus_export, focus_path)

cat(
  paste0(
    "PASS: finalized D1--D6 rapid Fig3 exports validated.\n",
    "Audit rows: ", nrow(audit), " -> ", output_path, "\n",
    "Direction rows: ", nrow(consistency), " -> ", consistency_path, "\n",
    "Focus-genus rows: ", nrow(focus_export), " -> ", focus_path, "\n",
    "Evidence boundary: ", evidence_boundary, "\n"
  )
)
