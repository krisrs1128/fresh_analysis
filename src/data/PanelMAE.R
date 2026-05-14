# PanelMAE: a MultiAssayExperiment on a common subject-time axis.
#
# For each sample,
#
#   t = (visit_date - enroll_date) / time_unit.
#
# Events are stored separately as (study_id, t_end, event).

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(tidyverse)
library(glue)

# ---- helpers ----------------------------------------------------------------

.unit_to_days <- function(unit) {
    switch(
        unit,
        days = 1,
        weeks = 7,
        stop("Unsupported time_unit: ", unit)
    )
}

.relative_time <- function(when, origin, unit) {
    as.numeric(difftime(when, origin, units = "days")) / .unit_to_days(unit)
}

# One row per (study_id, experiment, sample_col).
.build_visit_table <- function(mae, id_col, enroll_col, visit_time_col,
                               time_unit) {
    baseline <- as.data.frame(colData(mae))
    enroll_time <- setNames(
        baseline[[enroll_col]],
        as.character(baseline[[id_col]])
    )

    imap_dfr(experiments(mae), \(se, exp_name) {
        sample_data <- as.data.frame(colData(se))
        study_id <- sample_data[[id_col]]

        tibble(
            study_id = study_id,
            experiment = exp_name,
            sample_col = colnames(se),
            visit_date = sample_data[[visit_time_col]],
            t = .relative_time(
                sample_data[[visit_time_col]],
                enroll_time[study_id],
                time_unit
            )
        )
    })
}

# Re-express events on the same subject-time axis.
.normalize_events <- function(
    events,
    baseline,
    id_col,
    enroll_col,
    events_time_col,
    time_unit
) {
    events <- as_tibble(events)
    if (id_col != "study_id") {
        events <- events |>
            rename(study_id = all_of(id_col))
    }

    enroll_time <- baseline |>
        select(all_of(c(id_col, enroll_col))) |>
        rename(study_id = all_of(id_col), enroll_time = all_of(enroll_col))

    events |>
        left_join(enroll_time, by = "study_id") |>
        mutate(
            t_end = .relative_time(
                .data[[events_time_col]],
                enroll_time,
                time_unit
            ),
            event = as.integer(event)
        ) |>
        select(study_id, t_end, event)
}

# Warn if a recorded visit lies beyond t_end.
.validate_visit_times <- function(pma) {
    post_end_visits <- pma$visit_table |>
        left_join(pma$events, by = "study_id") |>
        filter(!is.na(t_end), t > t_end)

    if (nrow(post_end_visits) > 0) {
        warning(
            glue(
                "{nrow(post_end_visits)} visits in ",
                "{n_distinct(post_end_visits$study_id)} subjects ",
                "occur after t_end. See pma$visit_table."
            )
        )
    }
}

# ---- constructor ------------------------------------------------------------

# Map panel data onto a unified relative scale (time since enrollment)
#
# Inputs:
#   mae: MultiAssayExperiment with baseline colData and assay matrices.
#   events: Data frame defining (subject, event_time, status) tuples.
#   *_col, time_unit: Strings specifying column mappings and temporal scaling.
# Output: A PanelMAE object representing the MultiAssayExperiment, event
#   metadata, and visit table.
PanelMAE <- function(
    mae,
    events,
    id_col = "study_id",
    enroll_col = "enroll_date",
    visit_time_col = "visit_date",
    time_unit = "weeks",
    events_time_col = "event_date"
) {
    baseline <- as.data.frame(colData(mae))

    pma <- structure(
        list(
            mae = mae,
            events = .normalize_events(
                events,
                baseline,
                id_col,
                enroll_col,
                events_time_col,
                time_unit
            ),
            visit_table = .build_visit_table(
                mae,
                id_col,
                enroll_col,
                visit_time_col,
                time_unit
            ),
            id_col = id_col,
            enroll_col = enroll_col,
            visit_time_col = visit_time_col,
            time_unit = time_unit
        ),
        class = "PanelMAE"
    )

    .validate_visit_times(pma)
    pma
}

from_mae <- function(mae, events = metadata(mae)$events, ...) {
    if (is.null(events)) events <- metadata(mae)$event
    PanelMAE(mae, events, ...)
}

# ---- accessors --------------------------------------------------------------

events <- function(pma) UseMethod("events")
events.PanelMAE <- function(pma) pma$events
samples <- function(pma) UseMethod("samples")
samples.PanelMAE <- function(pma) pma$visit_table

# ---- subsetting -------------------------------------------------------------

subset_subjects <- function(pma, study_ids) UseMethod("subset_subjects")

subset_subjects.PanelMAE <- function(pma, study_ids) {
    pma$mae <- pma$mae[, study_ids, ]
    pma$events <- pma$events |> filter(study_id %in% study_ids)
    pma$visit_table <- pma$visit_table |> filter(study_id %in% study_ids)
    structure(pma, class = "PanelMAE")
}

subset_cohorts <- function(pma, cohort_ids, cohort_col = "cohort_id") {
    UseMethod("subset_cohorts")
}

subset_cohorts.PanelMAE <- function(pma, cohort_ids, cohort_col = "cohort_id") {
    keep <- colData(pma$mae) |>
        as.data.frame() |>
        filter(.data[[cohort_col]] %in% cohort_ids) |>
        pull(study_id)

    subset_subjects(pma, keep)
}

# ---- print ------------------------------------------------------------------

print.PanelMAE <- function(x, ...) {
    cat("PanelMAE\n")
    cat("  subjects:", nrow(x$events), "\n")
    cat("  events (1):", sum(x$events$event == 1), "\n")
    cat("  experiments:", paste(names(experiments(x$mae)), collapse = ", "), "\n")
    cat("  visits (rows):", nrow(x$visit_table), "\n")
    cat("  time_unit:", x$time_unit, "\n")
    invisible(x)
}
