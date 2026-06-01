# src/models/fit_poisson_interval.R
#
# Piecewise exponential survival model on an EventMAE's life_table
#
# Returns a `fit_pexp` object carrying:
#   glm                       — the stats::glm fit
#   cuts                      — sorted event times across people
#   cut_strategy              — how intervals were partitioned
#   time_varying              — forwarded from spec
#   impute_covariate_trajectory — forwarded from spec
#
# Choices for cut_strategy are:
#    event_times → fine step function
#     visits → coarser grid tied to sampling design. The grid travels with the
#       fit so predict_survival() reconstructs S(t | x) without the original
#       EventMAE.

library(tidyverse)

# ---- public -----------------------------------------------------------------

#' Fit a Piecewise-Exponential (Poisson GLM) Survival Model
#'
#' @param emae EventMAE object.
#' @param formula Poisson formula with event as outcome. Include
#'   interval on RHS for baseline hazard that changes over intervals
#'   Example:
#'   "event ~ age_at_start + interval".
#' @param cut_strategy "event_times" (fine grid) or '   "visits"
#    (coarser, tied to sampling).
#' @param time_varying RHS variable names that change by interval.
#'   Empty for baseline-only models.
#' @param impute_covariate_trajectory Covariate trajectory strategy:
#'   "extrapolate_baseline", "observed_path",
#'   "counterfactual_path", or "landmark".
#' @return S3 object of class fit_pexp.
fit_poisson_interval <- function(
    emae,
    formula,
    cut_strategy = c("event_times", "visits"),
    time_varying = character(),
    impute_covariate_trajectory = "extrapolate_baseline"
) {
    cut_strategy <- match.arg(cut_strategy)
    lt <- pexp_data(emae, cut_strategy)
    frm <- as.formula(formula)

    glm_fit <- glm(
        frm,
        offset = log(lt$t_end - lt$t_start),
        family = poisson,
        data = lt
    )

    structure(
        list(
            glm = glm_fit,
            cuts = sort(unique(c(lt$t_start, lt$t_end))),
            cut_strategy = cut_strategy,
            time_varying = time_varying,
            impute_covariate_trajectory = impute_covariate_trajectory
        ),
        class = "fit_pexp"
    )
}

# ---- print ------------------------------------------------------------------

print.fit_pexp <- function(x, ...) {
    cat("fit_pexp\n")
    cat("  engine: pexp (Poisson GLM)\n")
    cat("  formula: ", deparse(formula(x$glm)), "\n")
    cat("  cut_strategy: ", x$cut_strategy, "\n")
    cat("  n_cuts: ", length(x$cuts), "\n")
    cat("  time_varying: ", paste(x$time_varying, collapse = ", ") %|% "(none)", "\n")
    cat("  impute_covariate_trajectory: ", x$impute_covariate_trajectory, "\n")
    invisible(x)
}

# ---- internal ---------------------------------------------------------------

#' Build Poisson GLM Input Data
#'
#' One row per subject-interval with (study_id, t_start, t_end, event_status,
#' event, interval factor, baseline colData).
#'
#' @param emae EventMAE object.
#' @param cut_strategy Passed to life_table().
#' @return Tibble for glm().
pexp_data <- function(emae, cut_strategy) {
    lt <- life_table(emae, cut_strategy = cut_strategy) |>
        mutate(
            study_id = as.character(study_id),
            event = as.integer(event_status == "event_in_interval")
        )

    # Build interval factor with levels ordered by t_start
    interval_levels <- lt |>
        distinct(t_start, t_end) |>
        arrange(t_start) |>
        mutate(label = paste0("[", t_start, ",", t_end, ")")) |>
        pull(label)

    lt <- lt |>
        mutate(
            interval = factor(
                paste0("[", t_start, ",", t_end, ")"),
                levels = interval_levels
            )
        )

    baseline <- colData(emae$mae) |>
        as.data.frame() |>
        as_tibble() |>
        mutate(study_id = as.character(study_id))

    lt |> left_join(baseline, by = "study_id")
}

# Null-coalescing helper
`%|%` <- function(x, y) if (length(x) == 0 || is.null(x) || identical(x, "")) y else x
