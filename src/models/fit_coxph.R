# src/models/fit_coxph.R
#
# Fit a Cox proportional hazards model on the life_table from an EventMAE.
#
# Returns S3 class `fit_coxph` carrying:
#   coxph                     — the coxph fit
#   baseline_hazard           — cumulative baseline hazard
#   cuts                      — all observed interval boundaries
#   time_varying              — covariates that change over time
#   impute_covariate_trajectory — covariates for predicting future trajectories
#
# We save these fields so that predict_survival() needs no extra arguments
# beyond fit, newdata, and times.

library(survival)
library(tidyverse)

# ---- public -----------------------------------------------------------------

#' Fit a Cox Proportional Hazards Model
#'
#' @param emae EventMAE object.
#' @param formula Counting-process formula string, e.g.
#'   "Surv(t_start, t_end, event) ~ age_at_start".
#' @param time_varying RHS variable names that change by interval.
#'   Empty for baseline-only models.
#' @param impute_covariate_trajectory Covariate trajectory strategy:
#'   "extrapolate_baseline", "observed_path",
#'   "counterfactual_path", or "landmark".
#' @return S3 object of class fit_coxph.
fit_coxph <- function(
    emae,
    formula,
    time_varying = character(),
    impute_covariate_trajectory = "extrapolate_baseline"
) {
    lt <- coxph_data(emae, time_varying)
    frm <- as.formula(formula)
    cox <- coxph(frm, data = lt, ties = "efron")

    structure(
        list(
            coxph = cox,
            baseline_hazard = basehaz(cox, centered = FALSE),
            cuts = sort(unique(c(lt$t_start, lt$t_end))),
            time_varying = time_varying,
            impute_covariate_trajectory = impute_covariate_trajectory
        ),
        class = "fit_coxph"
    )
}

# ---- print ------------------------------------------------------------------

print.fit_coxph <- function(x, ...) {
    cat("fit_coxph\n")
    cat("  engine: coxph\n")
    cat("  formula: ", deparse(formula(x$coxph)), "\n")
    cat("  time_varying: ", paste(x$time_varying, collapse = ", ") %|% "(none)", "\n")
    cat("  impute_covariate_trajectory: ", x$impute_covariate_trajectory, "\n")
    cat("  n_subjects: ", x$coxph$n, "\n")
    cat("  n_events: ", x$coxph$nevent, "\n")
    cat("  support: [", min(x$cuts), ", ", max(x$cuts), "]\n", sep = "")
    invisible(x)
}

# ---- internal ---------------------------------------------------------------

#' Build coxph Input Data
#'
#' Extracts (study_id, t_start, t_end, event) from the life table and joins
#' baseline colData. Time-varying joins are reserved for future extension;
#' all covariates are currently baseline.
#'
#' @param emae EventMAE object.
#' @param time_varying Time-varying covariate names (unused).
#' @return Tibble ready for coxph().
coxph_data <- function(emae, time_varying) {
    lt <- life_table(emae, cut_strategy = "visits") |>
        mutate(
            study_id = as.character(study_id),
            event = as.integer(event_status == "event_in_interval")
        )

    baseline <- colData(emae$mae) |>
        as.data.frame() |>
        as_tibble() |>
        mutate(study_id = as.character(study_id))

    lt |> left_join(baseline, by = "study_id")
}

# Null-coalescing helper (avoids taking rlang as a hard dependency)
`%|%` <- function(x, y) if (length(x) == 0 || is.null(x) || identical(x, "")) y else x
