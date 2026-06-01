# life_table(emae, ...)
#
# Return a long table of left-open, right-closed intervals:
# (study_id, t_start, t_end, event_status).

library(tidyverse)

###############################################################################
## helper functions
###############################################################################

#' Partition Times per Subject
#'
#' @param events_data Tibble of per-subject event times.
#' @param visit_table Tibble of observed assay visit times.
#' @param cut_strategy "visits": split at each subject's own visit
#'   times. "event_times": split at all observed event times, which is useful
#'   for some semiparametric approaches.
#' @return Tibble mapping study_id to a list-column of cuts.
subject_cuts <- function(events_data, visit_table, cut_strategy) {
    if (cut_strategy == "visits") {
        return(
            visit_table |>
                distinct(study_id, t) |>
                arrange(study_id, t) |>
                summarise(cuts = list(t), .by = study_id)
        )
    } else if (cut_strategy == "event_times") {
        global_cuts <- events_data |>
            filter(event == 1) |>
            pull(t_end) |>
            unique() |>
            sort()

        return (
            tibble(
                study_id = events_data$study_id,
                cuts = rep(list(global_cuts), nrow(events_data))
            )
        )
    }
}

complete_cuts <- function(cut_vector, t_end) {
    interior_cuts <- cut_vector[cut_vector > 0 & cut_vector < t_end]
    unique(c(0, interior_cuts, t_end))
}

cuts_to_intervals <- function(cuts) {
    tibble(t_start = head(cuts, -1), t_end = tail(cuts, -1))
}

###############################################################################
## Define the actual life table constructors
###############################################################################

life_table <- function(emae, ...) UseMethod("life_table")

#' Build a Person-Period Table
#'
#' One row per subject-interval. Classifies each interval as containing an
#' event, at risk with no event, or post-event.
#'
#' @param emae EventMAE object.
#' @param cut_strategy "visits" or "event_times".
#' @param cuts Optional additional cut points.
#' @param include_post_event_rows Retain rows after the terminal event
#'   interval? Default FALSE.
#' @return Tibble with columns study_id, t_start, t_end,
#'   event_status (one of "event_in_interval",
#'   "at_risk_no_event", "past_event").
life_table.EventMAE <- function(
    emae,
    cut_strategy = c("visits", "event_times"),
    cuts = NULL,
    include_post_event_rows = FALSE
) {
    cut_strategy <- match.arg(cut_strategy)
    events_data <- emae$events
    subject_cuts <- subject_cuts(events_data, emae$visit_table, cut_strategy)

    interval_data <- events_data |>
        left_join(subject_cuts, by = "study_id") |>
        mutate(
            cuts = map2(cuts, t_end, complete_cuts),
            intervals = map(cuts, cuts_to_intervals)
        )

    result <- interval_data |>
        select(study_id, subject_t_end = t_end, event, intervals) |>
        unnest(intervals) |>
        mutate(
            event_status = case_when(
                t_start >= subject_t_end ~ "past_event",
                event == 1 & t_end == subject_t_end ~ "event_in_interval",
                TRUE ~ "at_risk_no_event"
            )
        ) |>
        select(study_id, t_start, t_end, event_status)

    if (!include_post_event_rows) {
        result <- result |> filter(event_status != "past_event")
    }

    result
}
