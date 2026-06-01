# view(emae, t_star, ...)
#
# Filter or summarize EventMAEs according to a reference time t_star. For
# example, when `strategy` is 'last', we just keep the last visit before time
# t_star for each of the MAEs. We add (study_id, t) to the experiment-level
# colData's so that we retain a memory of which rows from the original
# experiments are kept.

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)

# ---- helpers ----------------------------------------------------------------

# Keep one column per subject (`last`) or all eligible columns (`all`).
# Internal wrapper around `dplyr::slice_max()`.
apply_strategy <- function(eligible_samples, strategy) {
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

# Get an experiment and attach the time-since-enrollment t.
#
# Also adds the (study_id, t) information so that we don't loose that context
# after subsetting.
subset_experiment <- function(emae, exp_name, sample_cols) {
    se <- experiments(emae$mae)[[exp_name]][, sample_cols, drop = FALSE]

    visit_lookup <- emae$visit_table |>
        filter(experiment == exp_name) |>
        select(sample_col, t)

    colData(se)$study_id <- colData(se)[[emae$id_col]]
    colData(se)$t <- visit_lookup$t[match(colnames(se), visit_lookup$sample_col)]
    se
}

# ---- main -------------------------------------------------------------------

view <- function(emae, ...) UseMethod("view")

# Map each subject to their latest valid assay prior to `t_star`.
# Returns a subset of an MultiAssayExperiment.
#
# Inputs:
#   emae: A EventMAE object.
#   t_star: Upper bound on allowed visit times.
#   strategy: How to reduce multiple observations to a single summary.
# Output: An MultiAssayExperiment mapping subjects to a filtered assays.
view.EventMAE <- function(
    emae,
    t_star = Inf,
    strategy = "last",
    summary_fns = NULL,
    assays = NULL,
    drop_post_event = TRUE
) {
    if (is.null(assays)) {
        assays <- names(experiments(emae$mae))
    }

    eligible_visits <- emae$visit_table |>
        filter(experiment %in% assays, t <= t_star)

    if (drop_post_event) {
        eligible_visits <- eligible_visits |>
            left_join(emae$events |>
            select(study_id, t_end), by = "study_id") |>
            filter(is.na(t_end) | t <= t_end) |>
            select(-t_end)
    }

    visits_by_experiment <- split(eligible_visits, eligible_visits$experiment)
    experiments_out <- imap(visits_by_experiment, \(visit_data, exp_name) {
        sample_cols <- apply_strategy(visit_data, strategy)
        subset_experiment(emae, exp_name, sample_cols)
    })

    sample_map <- listToMap(lapply(experiments_out, \(se) {
        DataFrame(
            primary = as.character(colData(se)$study_id),
            colname = colnames(se)
        )
    }))

    primary_data <- colData(emae$mae)
    primary_data$study_id <- rownames(primary_data)

    MultiAssayExperiment(
        experiments = experiments_out,
        colData = primary_data,
        sampleMap = sample_map
    )
}
