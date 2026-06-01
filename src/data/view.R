# view()
#
# Filter EventMAE to visits at or before t_star. Adds (study_id, t) to
# experiment colData so we see which intervals are present.

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)

# ---- helpers ----------------------------------------------------------------

#' Keep One Sample per Subject or All Eligible Samples
#'
#' Wraps `slice_max` for the time variable.
#'
#' @param eligible_samples Tibble of candidate sample rows.
#' @param strategy "last" for the most recent visit per subject, or
#'   "all" to retain every eligible visit.
#' @return Character vector of sample_col values to retain.
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

#' Subset an Experiment and Annotate with Relative Time
#'
#' Extracts a SummarizedExperiment by sample columns then annotates colData with
#' study_id and t.
#'
#' @param emae EventMAE object.
#' @param exp_name Experiment name.
#' @param sample_cols Sample column names to retain.
#' @return SummarizedExperiment with study_id and t in colData.
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

#' Subset EventMAE to Valid Assay Observations
#'
#' Filters each experiment to visits at or before t_star, then reduces
#' multiple observations per subject by strategy.
#'
#' @param emae EventMAE object.
#' @param t_star Upper bound on visit times. Default Inf.
#' @param strategy "last" (most recent visit) or "all"
#'   (every eligible visit).
#' @param summary_fns Not implemented yet...
#' @param assays Experiments to include. Default: all.
#' @param drop_post_event Exclude post-event visits? Default TRUE.
#' @return MultiAssayExperiment with filtered assays.
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
