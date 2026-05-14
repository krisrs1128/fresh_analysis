# life_table(pma, ...)
#
# Return a long table of left-open, right-closed intervals:
# (study_id, t_start, t_end, event_status).

library(tidyverse)

###############################################################################
## helper functions
###############################################################################

# Partition timess for each subject
#
# Inputs:
#   events_data: Tibble containing final event times per subject.
#   visit_table: Tibble containing observed visit times for assays.
#   cut_strategy: 'visits' (split a subject according to visit times) or
#     'event_times' (split a at all the observed event times in the data, useful
#     for some types of models).
# Output: A tibble mapping `study_id` to a list-column of temporal `cuts`.
.subject_cuts <- function(events_data, visit_table, cut_strategy) {
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

.complete_cuts <- function(cut_vector, t_end) {
    interior_cuts <- cut_vector[cut_vector > 0 & cut_vector < t_end]
    unique(c(0, interior_cuts, t_end))
}

.cuts_to_intervals <- function(cuts) {
    tibble(
        t_start = head(cuts, -1),
        t_end = tail(cuts, -1)
    )
}

###############################################################################
## Define the actual life table constructors
###############################################################################

life_table <- function(pma, ...) UseMethod("life_table")

# Make a long-format table mapping intervals (either between visits or events)
# to event outcomes.
#
# Inputs:
#   pma: A PanelMAE object defining events and visit times.
#   cut_strategy: Split temporal intervals by 'visits' or 'event_times'.
# Output: A tibble mapping (study_id) to partitions (t_start, t_end) with
#   `event_status` determined by whether an event has happened before the
#   interval started, during the current interval, or whether the subject is
#   still at risk.
life_table.PanelMAE <- function(
    pma,
    cut_strategy = c("visits", "event_times"),
    cuts = NULL,
    include_post_event_rows = FALSE
) {
    cut_strategy <- match.arg(cut_strategy)
    events_data <- pma$events
    subject_cuts <- .subject_cuts(events_data, pma$visit_table, cut_strategy)

    interval_data <- events_data |>
        left_join(subject_cuts, by = "study_id") |>
        mutate(
            cuts = map2(cuts, t_end, .complete_cuts),
            intervals = map(cuts, .cuts_to_intervals)
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
