# src/models/spec.R
#
# A spec is a plain YAML with keys defining the model to fit. The `fit_engine`
# function wraps specific model fitting definitions from the `models` directory.
# Every key except `engine` must be a named argument of the corresponding fit
# function:
#   do.call(fit_fn, c(list(emae), spec[names(spec) != "engine"]))

library(yaml)
library(here)
library(fs)

# ---- read -------------------------------------------------------------------

#' Read a Model Spec from YAML
#'
#' @param path Path to a YAML file (schema in \code{conf/models/README.md}).
#' @return Named list forwarded as arguments to the fit function.
read_spec <- function(path) {
    spec <- yaml::read_yaml(path)
    if (is.null(spec$time_varying)) {
        spec$time_varying <- character(0)
    } else {
        spec$time_varying <- as.character(spec$time_varying)
    }
    if (is.null(spec$impute_covariate_trajectory)) {
        spec$impute_covariate_trajectory <- "extrapolate_baseline"
    }
    spec
}

# ---- dispatch ---------------------------------------------------------------

# Load source files for all fit functions (idempotent via source()).
load_fit_fns <- function() {
    source(here("src", "models", "fit_coxph.R"))
    source(here("src", "models", "fit_poisson_interval.R"))
}

#' Read Spec and Fit
#'
#' Forwards every spec key except \code{engine} as a named argument to the
#' engine-specific fit function. No mapping layer.
#'
#' @param emae EventMAE object.
#' @param spec Named list from \code{read_spec()}.
#' @return S3 fit object from the engine-specific fit function.
fit_engine <- function(emae, spec) {
    load_fit_fns()

    fit_fn <- switch(
        spec$engine,
        coxph = fit_coxph,
        pexp = fit_poisson_interval,
        stop(
            "Engine '", spec$engine, "' is not yet implemented.\n",
            "See src/models/fit_<engine>.R and conf/models/README.md."
        )
    )

    # Drop `engine`; every remaining key is a named argument of fit_fn
    args <- spec[names(spec) != "engine"]
    do.call(fit_fn, c(list(emae), args))
}
