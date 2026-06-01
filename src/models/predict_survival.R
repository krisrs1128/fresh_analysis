# src/models/predict_survival.R
#
# predict_survival(fit, newdata, times)
#
#   fit     — the fitted object, from whatever model fitting engine of interest
#             (fit_coxph, fit_pexp, …)
#   newdata — covariate trajectory on which to make survival probability
#             predictions: one row per (subject, interval) with study_id,
#             t_start, t_end, and formula RHS variables. Use
#             constant_trajectory() to extrapolate a baseline row.
#   times   — shared engine-independent grid. Every method evaluates
#             S(t | x) at these times; returns n × length(times) matrix.

library(tidyverse)

#' Predict Survival Curves on a Shared Time Grid
#'
#' @param fit S3 fit object (\code{fit_coxph}, \code{fit_pexp}, …).
#' @param newdata Covariate trajectory: one row per (subject, interval) with
#'   \code{study_id}, \code{t_start}, \code{t_end}, and RHS variables.
#'   Use \code{constant_trajectory()} for baseline-only models.
#' @param times Numeric vector of evaluation times.
#' @return Matrix (n_subjects × length(times)) of S(t | x) in [0, 1]. Rows
#'   named by \code{study_id}; columns named by \code{times}.
predict_survival <- function(fit, newdata, times) {
    UseMethod("predict_survival")
}

# ---- fit_coxph method -------------------------------------------------------

#' Predict S(t | x) from a fit_coxph Object
#'
#' Cumulative baseline hazards over trajectory intervals, weighting each
#' increment by the interval's linear predictor,
#'
#'   H(t | x(·)) = Σ_k [H0(min(b_k, t)) − H0(a_k)] · exp(x_k β)
#'   S(t | x(·)) = exp(−H(t | x(·)))
#'
#' where (a_k, b_k] are trajectory intervals up to t. For constant covariates
#' this reduces to S(t) = exp(−H0(t) · exp(xβ)).
#'
#' S(t) held constant past the last observed event time.
#'
#' @param fit \code{fit_coxph} object.
#' @param newdata Covariate trajectory tibble.
#' @param times Evaluation times.
#' @return Matrix of survival probabilities.
predict_survival.fit_coxph <- function(fit, newdata, times) {
    max_support <- max(fit$baseline_hazard$time)

    # Cumulative baseline hazard, from basehaz
    H0 <- stepfun(fit$baseline_hazard$time, c(0, fit$baseline_hazard$hazard))
    subjects <- unique(newdata$study_id)
    n <- length(subjects)

    S_mat <- matrix(
        NA_real_,
        nrow = n,
        ncol = length(times),
        dimnames = list(subjects, as.character(times))
    )

    for (i in seq_along(subjects)) {
        sid <- subjects[i]
        traj <- newdata |> filter(study_id == sid)

        # Linear predictor for each interval row
        xb <- predict(fit$coxph, newdata = traj, type = "lp")

        for (j in seq_along(times)) {
            t_q <- times[j]
            active <- which(traj$t_start < t_q)

            if (length(active) == 0) {
                S_mat[i, j] <- 1
                next
            }

            a_k <- traj$t_start[active]
            b_k <- pmin(traj$t_end[active], t_q)

            # Clamp to the support before evaluating H0
            H_incr <- H0(pmin(b_k, max_support)) - H0(pmin(a_k, max_support))
            H <- sum(H_incr * exp(xb[active]))
            S_mat[i, j] <- exp(-H)
        }
    }

    S_mat
}

# ---- fit_pexp method (stub) -------------------------------------------------

#' Predict S(t | x) from a fit_pexp Object
#'
#' Plan:
#' 1. Per-interval hazard rate:
#'      lambda_k = predict(fit$glm, traj_k, type = "response") / width_k
#' 2. H(t) = Σ_k lambda_k · (min(t_end_k, t) − t_start_k), t_start_k < t
#' 3. S_cuts = exp(−H) at cut boundaries; interpolate onto \code{times}.
#'
#' Stub ensures dispatch chain is complete and fails informatively.
#'
#' @param fit \code{fit_pexp} object.
#' @param newdata Covariate trajectory tibble.
#' @param times Evaluation times.
#' @return Not implemented; signals error.
predict_survival.fit_pexp <- function(fit, newdata, times) {
    stop("predict_survival.fit_pexp is not yet implemented.")
}

# ---- trajectory helpers -----------------------------------------------------

#' Build a Constant Trajectory to Predict Baselines
#'
#' Repeats a baseline row once per interval defined by \code{cuts}, holding
#' all covariates constant at entry values. Input to \code{predict_survival()}
#' when \code{impute_covariate_trajectory = "extrapolate_baseline"}.
#'
#' Pre-select \code{baseline_row} to formula columns only.
#'
#' @param baseline_row One-row tibble with \code{study_id} and formula
#'   covariates.
#' @param cuts Interval boundaries (e.g. \code{fit$cuts}), ≥ 2 distinct
#'   values.
#' @return Tibble with \code{t_start}, \code{t_end}, \code{study_id}, and all
#'   covariate columns. One row per interval; covariates constant.
constant_trajectory <- function(baseline_row, cuts) {
    cuts <- sort(unique(cuts))
    intervals <- tibble(t_start = head(cuts, -1), t_end = tail(cuts, -1))
    covar_rows <- baseline_row[rep(1L, nrow(intervals)), , drop = FALSE]
    bind_cols(intervals, covar_rows)
}
