# design_matrix(view_output)
#
# Flatten the output of view() to a modeling table. Rows are indexed by
# study_id, or by (study_id, t) when visits are retained. Baseline colData are
# carried once per subject.

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)

# ---- main -------------------------------------------------------------------

# Merge baseline attributes and experiment matrices into one wide table.
# Redundant ID columns are removed to prevent namespace conflicts.
#
# Inputs:
#   view_output: A MultiAssayExperiment (typically the output of view()).
# Output: A tibble where rows are uniquely identified by (study_id) or
#   (study_id, t), binding baseline characteristics to flattened assay variables.
design_matrix <- function(view_output) {
    # extract tables that will go into the merged design
    baseline_df <- as.data.frame(colData(view_output)) |>
        as_tibble()
    experiment_tbls <- imap(
        experiments(view_output),
        .extract_experiment_features
    )
    key_cols <- .key_columns(experiment_tbls)
    combined_df <- .join_experiments(experiment_tbls, key_cols)

    # Join on key_cols, and remove redundant columns from combined_df to avoid .x/.y suffixes
    shared_cols <- intersect(names(baseline_df), names(combined_df))
    cols_to_drop <- setdiff(shared_cols, key_cols)
    combined_df <- combined_df |> select(-any_of(cols_to_drop))
    baseline_df |>
        right_join(combined_df, by = intersect(names(baseline_df), key_cols))
}

# ---- helpers ----------------------------------------------------------------

.extract_sample_data <- function(se) {
    keep_cols <- intersect(names(colData(se)), c("study_id", "t"))
    as.data.frame(colData(se)[, keep_cols, drop = FALSE]) |>
        rownames_to_column("sample_col") |>
        as_tibble()
}

.extract_assay_data <- function(se) {
    as.data.frame(t(assay(se))) |>
        rownames_to_column("sample_col") |>
        as_tibble()
}

.extract_experiment_features <- function(se, exp_name) {
    # Keep only columns necessary for joining (study_id, and optionally t)
    sample_data <- .extract_sample_data(se)
    assay_data <- .extract_assay_data(se)

    sample_data |>
        left_join(assay_data, by = "sample_col") |>
        select(-sample_col)
}

.key_columns <- function(experiment_tbls) {
    if ("t" %in% names(experiment_tbls[[1]])) {
        c("study_id", "t")
    } else {
        "study_id"
    }
}

.join_experiments <- function(experiment_tbls, key_cols) {
    reduce(experiment_tbls, ~ full_join(.x, .y, by = intersect(names(.x), names(.y))))
}

model_frame <- function(pma, ...) {
    design_matrix(view(pma, ...))
}