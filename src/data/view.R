# view(pma, t_star, ...)
#
# Materialize covariates up to time t_star. Each experiment is column-subset to
# the visits retained by `strategy`. Experiment-level colData are augmented with
# the standardized pair (study_id, t). MAE-level colData remain subject-level.

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)

# ---- helpers ----------------------------------------------------------------

# Keep one column per subject (`last`) or all eligible columns (`all`).
# Internal wrapper around `dplyr::slice_max()`.
.apply_strategy <- function(eligible_samples, strategy) {
    switch(
        strategy,
        last = eligible_samples |>
            group_by(study_id) |>
            slice_max(t, n = 1, with_ties = FALSE) |>
            ungroup() |>
            pull(sample_col),
        all = eligible_samples$sample_col,
        summary = stop("strategy='summary' not implemented"),
        locf = stop("strategy='locf' not implemented yet")
    )
}

# Slice one experiment and attach the standardized time coordinate t.
# Extracts the desired subset and binds `visit_lookup` mappings.
.subset_experiment <- function(pma, exp_name, sample_cols) {
    se <- experiments(pma$mae)[[exp_name]][, sample_cols, drop = FALSE]

    visit_lookup <- pma$visit_table |>
        filter(experiment == exp_name) |>
        select(sample_col, t)

    colData(se)$study_id <- colData(se)[[pma$id_col]]
    colData(se)$t <- visit_lookup$t[match(colnames(se), visit_lookup$sample_col)]
    se
}

# ---- main -------------------------------------------------------------------

view <- function(pma, ...) UseMethod("view")

# Map each subject to their latest valid assay prior to `t_star`.
# Constructs a standard MultiAssayExperiment subset.
#
# Inputs:
#   pma: A PanelMAE object.
#   t_star: Upper bound on valid observation times.
#   strategy: Resolution strategy for multiple observations (e.g., 'last').
# Output: A standard MultiAssayExperiment mapping subjects to resolved assays.
view.PanelMAE <- function(
    pma,
    t_star = Inf,
    strategy = "last",
    summary_fns = NULL,
    assays = NULL,
    drop_post_event = TRUE
) {
    if (is.null(assays)) {
        assays <- names(experiments(pma$mae))
    }

    eligible_visits <- pma$visit_table |>
        filter(experiment %in% assays, t <= t_star)

    if (drop_post_event) {
        eligible_visits <- eligible_visits |>
            left_join(pma$events |>
            select(study_id, t_end), by = "study_id") |>
            filter(is.na(t_end) | t <= t_end) |>
            select(-t_end)
    }

    visits_by_experiment <- split(eligible_visits, eligible_visits$experiment)
    experiments_out <- imap(visits_by_experiment, \(visit_data, exp_name) {
        sample_cols <- .apply_strategy(visit_data, strategy)
        .subset_experiment(pma, exp_name, sample_cols)
    })

    sample_map <- listToMap(lapply(experiments_out, \(se) {
        DataFrame(
            primary = as.character(colData(se)$study_id),
            colname = colnames(se)
        )
    }))

    primary_data <- colData(pma$mae)
    primary_data$study_id <- rownames(primary_data)

    MultiAssayExperiment(
        experiments = experiments_out,
        colData = primary_data,
        sampleMap = sample_map
    )
}
