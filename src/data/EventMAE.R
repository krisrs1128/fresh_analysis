# EventMAE: MultiAssayExperiment for Longitudinal Time Courses
#
# This class,
# - converts dates into relative time from enrollment
# - distinguishes between event descriptors and visit data.
#
# The relative times are standardized according to,
#
#   t = (visit_date - enroll_date) / time_unit.
#
# Events are stored as (study_id, t_end, event_type). The intent is to simplify
# survival analysis on MAE objects.

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(tidyverse)
library(glue)

# ---- helpers ----------------------------------------------------------------

relative_time <- function(when, origin, unit) {
    relative_unit <- switch(
        unit,
        days = 1,
        weeks = 7,
        stop("Unsupported time_unit: ", unit)
    )
    as.numeric(difftime(when, origin, units = "days")) / relative_unit
}

# One row per (study_id, experiment, sample_col).
build_visit_table <- function(mae, id_col, enroll_col, visit_time_col,
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
            t = relative_time(
                sample_data[[visit_time_col]],
                enroll_time[study_id],
                time_unit
            )
        )
    })
}

# Convert event times to standardized relative-to-enrollment units
standardize_event_times <- function(
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

    # find the enrollment time
    enroll_time <- baseline |>
        select(all_of(c(id_col, enroll_col))) |>
        rename(study_id = all_of(id_col), enroll_time = all_of(enroll_col))

    # compute standardized times
    events |>
        left_join(enroll_time, by = "study_id") |>
        mutate(
            t_end = relative_time(
                .data[[events_time_col]],
                enroll_time,
                time_unit
            ),
            event = as.integer(event)
        ) |>
        select(study_id, t_end, event)
}

# Warn if a recorded visit lies beyond t_end.
validate_visit_times <- function(emae) {
    post_end_visits <- emae$visit_table |>
        left_join(emae$events, by = "study_id") |>
        filter(!is.na(t_end), t > t_end)

    if (nrow(post_end_visits) > 0) {
        warning(
            glue(
                "{nrow(post_end_visits)} visits in ",
                "{n_distinct(post_end_visits$study_id)} subjects ",
                "occur after t_end. See emae$visit_table."
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
# Output: A EventMAE object representing the MultiAssayExperiment, event
#   metadata, and visit table.
EventMAE <- function(
    mae,
    events,
    id_col = "study_id",
    enroll_col = "enroll_date",
    visit_time_col = "visit_date",
    time_unit = "days",
    events_time_col = "event_date"
) {
    baseline <- as.data.frame(colData(mae))

    emae <- structure(
        list(
            mae = mae,
            events = standardize_event_times(
                events,
                baseline,
                id_col,
                enroll_col,
                events_time_col,
                time_unit
            ),
            visit_table = build_visit_table(
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
        class = "EventMAE"
    )

    validate_visit_times(emae)
    emae
}

from_mae <- function(mae, events = metadata(mae)$events, ...) {
    if (is.null(events)) events <- metadata(mae)$event
    EventMAE(mae, events, ...)
}

# ---- accessors --------------------------------------------------------------

events <- function(emae) UseMethod("events")
events.EventMAE <- function(emae) emae$events
samples <- function(emae) UseMethod("samples")
samples.EventMAE <- function(emae) emae$visit_table

# ---- subsetting -------------------------------------------------------------

subset_subjects <- function(emae, study_ids) UseMethod("subset_subjects")

subset_subjects.EventMAE <- function(emae, study_ids) {
    emae$mae <- emae$mae[, study_ids, ]
    emae$events <- emae$events |> filter(study_id %in% study_ids)
    emae$visit_table <- emae$visit_table |> filter(study_id %in% study_ids)
    structure(emae, class = "EventMAE")
}

subset_cohorts <- function(emae, cohort_ids, cohort_col = "cohort_id") {
    UseMethod("subset_cohorts")
}

subset_cohorts.EventMAE <- function(emae, cohort_ids, cohort_col = "cohort_id") {
    keep <- colData(emae$mae) |>
        as.data.frame() |>
        filter(.data[[cohort_col]] %in% cohort_ids) |>
        pull(study_id)

    subset_subjects(emae, keep)
}

# ---- print ------------------------------------------------------------------

print.EventMAE <- function(x, ...) {
    cat("EventMAE\n")
    cat("  subjects:", nrow(x$events), "\n")
    cat("  events (1):", sum(x$events$event == 1), "\n")
    cat("  experiments:", paste(names(experiments(x$mae)), collapse = ", "), "\n")
    cat("  visits (rows):", nrow(x$visit_table), "\n")
    cat("  time_unit:", x$time_unit, "\n")
    invisible(x)
}
