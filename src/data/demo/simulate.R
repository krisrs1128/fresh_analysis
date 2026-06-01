# Simulate a toy FRESH-style dataset for EventMAE development.
#
# Schemas match project notes. No feature–outcome signal; purely a pipeline
# sketch. Call build_example_mae() to assemble.
#
# Usage:
#   mae <- build_example_mae(seed = 202605)
#   saveRDS(mae, "data/example_mae.rds")

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(tidyverse)


# ---- Helper functions ------------------------------------------------------

random_factor <- function(n, levels, prob = NULL) {
    factor(sample(levels, n, replace = TRUE, prob = prob), levels = levels)
}

random_integer <- function(n, lo, hi) {
    as.integer(sample(lo:hi, n, replace = TRUE))
}


# ---- 1. Simulate participants & cohorts ------------------------------------

simulate_participants <- function(
    n_cohorts = 4,
    per_cohort = 25,
    start_date = as.Date("2018-01-01")
) {
    cohort_id   <- rep(sprintf("C%02d", seq_len(n_cohorts)), each = per_cohort)
    n          <- length(cohort_id)
    study_id   <- factor(sprintf("S%04d", seq_len(n)))

    # Each cohort starts ~6 months apart, enroll dates jittered within cohort
    cohort_start <- start_date + (seq_len(n_cohorts) - 1) * 180
    enroll_date  <- cohort_start[as.integer(factor(cohort_id))] +
                    sample(0:60, n, replace = TRUE)

    tibble(
        study_id = study_id,
        cohort_id = factor(cohort_id),
        enroll_date = enroll_date
    )
}


# ---- 2. Simulate baseline survey -------------------------------------------

simulate_baseline <- function(participants) {
    n <- nrow(participants)

    have_children <- random_factor(n, c("no", "yes"))
    n_children <- ifelse(have_children == "yes",
                        random_integer(n, 1L, 4L), 0L) |> as.integer()

    prior_sex <- random_factor(n, c("no", "yes"), prob = c(0.15, 0.85))
    age_first_sex <- ifelse(prior_sex == "yes",
                           random_integer(n, 14L, 22L), NA_integer_) |> as.integer()

    tibble(
        study_id = participants$study_id,
        age_at_start = random_integer(n, 18L, 23L),
        have_children = have_children,
        n_children = n_children,
        how_many_children = factor(
            ifelse(n_children == 0, "no_children",
            ifelse(n_children == 1, "one_child", "two_or_more")),
            levels = c("no_children", "one_child", "two_or_more")
        ),
        age_first_sex = age_first_sex,
        prior_sex = prior_sex,
        age_first_partner = ifelse(
            prior_sex == "yes", random_integer(n, 14L, 24L), NA_integer_)
            |> as.integer(),
        age_current_partner = random_integer(n, 17L, 40L),
        number_sexual_partners = random_integer(n, 0L, 6L),
        how_often_drink_alcohol = random_factor(n,
            c("never", "rarely", "sometimes", "often"),
            prob = c(0.4, 0.3, 0.2, 0.1)),
        drink_alcohol_3months = random_factor(n,
            c("a_yes_alcohol_3months", "b_no_alcohol_3months")),
        regular_partner = random_factor(n, c("no", "yes"), prob = c(0.3, 0.7)),
        regular_bf = random_factor(n, c("no", "yes"), prob = c(0.4, 0.6)),
        often_regular_partner = random_factor(n,
            c("almost_everyday", "few_times_week", "few_times_month", "rarely")),
        casual_partner = random_factor(n, c("no", "yes"), prob = c(0.7, 0.3)),
        often_casual_partner = random_factor(n,
            c("almost_everyday", "few_times_week", "few_times_month", "rarely"))
    )
}


# ---- 3. Simulate checkpoint survey (longitudinal) --------------------------

simulate_checkpoints <- function(
    participants,
    max_visits = 6,
    interval_days = 91
) {
    rows <- map_dfr(seq_len(nrow(participants)), function(i) {
        n_visits      <- sample(4:max_visits, 1)
        visit_number  <- seq_len(n_visits)
        visit_date    <- participants$enroll_date[i] + (visit_number - 1L) * interval_days
        tibble(
            study_id = participants$study_id[i],
            visit_number = as.integer(visit_number),
            visit_date = visit_date,
            visit_date_baseline = participants$enroll_date[i]
        )
    })

    n     <- nrow(rows)
    vag <- pmax(0, round(rnorm(n, 4, 3)))
    anal <- pmax(0, round(rnorm(n, 0.3, 0.7)))
    oral <- pmax(0, round(rnorm(n, 1, 1.5)))
    sex_past3 <- as.integer(vag + anal + oral)
    contraception_lv <- c(
        "depoprovera", "oral_pill", "implant", "iud",
        "condoms_only", "none"
    )

    rows |>
        mutate(
            nbr_people_fulltime = pmax(0, round(rnorm(n, 2, 1.5))),
            nbr_people_parttime = pmax(0, round(rnorm(n, 1, 1.0))),
            education = random_factor(n, c(
                "did_not_finish_matric", "matric", "some_tertiary", "tertiary")),
            age_current_partner = random_integer(n, 17L, 45L),
            sex_past3month = sex_past3,
            vaginal_sex_last30days = as.numeric(vag),
            anal_sex_last30days = as.numeric(anal),
            oral_sex_last30days = as.numeric(oral),
            condom_use = random_factor(n, c("always", "never", "sometimes")),
            female_condom = random_factor(n, c("always", "never", "sometimes"),
                                         prob = c(0.05, 0.85, 0.1)),
            times_pregnant = as.numeric(random_integer(n, 0L, 3L)),
            live_births = as.numeric(random_integer(n, 0L, 2L)),
            attending_clinic = random_factor(n, c("no", "yes"), prob = c(0.2, 0.8)),
            current_contraception = random_factor(n, contraception_lv),
            start_date_contraception = visit_date_baseline - sample(0:1500, n, replace = TRUE)
        ) |>
        select(
            study_id, visit_date_baseline, nbr_people_fulltime,
            nbr_people_parttime, education, age_current_partner,
            sex_past3month, vaginal_sex_last30days, anal_sex_last30days,
            oral_sex_last30days, condom_use, female_condom, times_pregnant,
            live_births, attending_clinic, current_contraception,
            start_date_contraception, visit_number, visit_date
        )
}


# ---- 4. Simulate 16S counts ------------------------------------------------

simulate_16s <- function(checkpoints, n_taxa = 15, lib_size = 5000) {
    n_samp <- nrow(checkpoints)
    taxa   <- sprintf("ASV%03d", seq_len(n_taxa))

    # Taxon-specific mean abundances on the log scale
    log_mu <- rnorm(n_taxa, mean = log(lib_size / n_taxa), sd = 1.2)
    counts <- vapply(seq_len(n_samp), function(j) {
        as.integer(rnbinom(n_taxa, mu = exp(log_mu), size = 2))
    }, integer(n_taxa))

    rownames(counts) <- taxa
    colnames(counts) <- sprintf(
        "%s_v%d",
        as.character(checkpoints$study_id),
        checkpoints$visit_number
    )

    list(
        counts = counts,
        coldata = DataFrame(
            sample_id = colnames(counts),
            study_id = checkpoints$study_id,
            visit_number = checkpoints$visit_number,
            visit_date = checkpoints$visit_date
        ),
        rowdata = DataFrame(
            taxon = taxa,
            phylum = sample(
                c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria"),
                n_taxa, replace = TRUE
            ),
            row.names = taxa
        )
    )
}


# ---- 5. Simulate HIV acquisition events ------------------------------------

simulate_events <- function(
    participants, checkpoints,
    event_rate = 0.08
) {
    last_visit <- checkpoints |>
        group_by(study_id) |>
        summarise(last_date = max(visit_date), .groups = "drop")

    df <- participants |>
        select(study_id, enroll_date) |>
        left_join(last_visit, by = "study_id")

    n     <- nrow(df)
    event <- as.integer(runif(n) < event_rate)

    # Generate event dates:
    # - If event==1: random date between enrollment and last visit
    # - If event==0: last visit date (censoring time)
    event_date <- ifelse(
        event == 1L,
        df$enroll_date + as.integer(runif(n, 1, as.numeric(df$last_date - df$enroll_date))),
        df$last_date
    ) |> as.Date(origin = "1970-01-01")

    tibble(
        study_id = df$study_id,
        event_date = event_date,
        event = event
    )
}


# ---- 6. Assemble MultiAssayExperiment (MAE) --------------------------------

build_example_mae <- function(seed = 202605, ...) {
    set.seed(seed)

    # Generate core tables
    participants <- simulate_participants(...)
    baseline     <- simulate_baseline(participants)
    checkpoints  <- simulate_checkpoints(participants)
    microbiome     <- simulate_16s(checkpoints)
    events       <- simulate_events(participants, checkpoints)

    # Primary colData: one row per study_id (baseline + cohort_id)
    primary <- baseline |>
        left_join(participants |> select(study_id, cohort_id, enroll_date),
            by = "study_id"
        ) |>
        as.data.frame()
    rownames(primary) <- as.character(primary$study_id)

    # "checkpoint" experiment: SE with trivial assay, info in colData
    cp_coldata <- DataFrame(checkpoints,
        row.names = sprintf(
            "%s_v%d",
            as.character(checkpoints$study_id),
            checkpoints$visit_number
        )
    )
    cp_assay <- matrix(checkpoints$visit_number,
        nrow = 1,
        dimnames = list("visit_number", rownames(cp_coldata))
    )
    cp_se <- SummarizedExperiment(
        assays = list(visit = cp_assay),
        colData = cp_coldata
    )

    # "microbiome" experiment
    s16_se <- SummarizedExperiment(
        assays = list(counts = microbiome$counts),
        colData = microbiome$coldata,
        rowData = microbiome$rowdata
    )
    colnames(s16_se) <- microbiome$coldata$sample_id

    # sampleMap
    cp_map <- data.frame(
        assay = "checkpoint",
        primary = as.character(checkpoints$study_id),
        colname = rownames(cp_coldata),
        stringsAsFactors = FALSE
    )
    s16_map <- data.frame(
        assay = "microbiome",
        primary = as.character(microbiome$coldata$study_id),
        colname = microbiome$coldata$sample_id,
        stringsAsFactors = FALSE
    )
    smap <- listToMap(list(
        checkpoint = cp_map[, c("primary", "colname")],
        microbiome = s16_map[, c("primary", "colname")]
    ))

    mae <- MultiAssayExperiment(
        experiments = ExperimentList(checkpoint = cp_se, microbiome = s16_se),
        colData = DataFrame(primary),
        sampleMap = smap,
        metadata = list(
            events = events,
            sim_info = list(
                seed = seed,
                n_subjects = nrow(participants),
                n_cohorts = length(unique(participants$cohort_id)),
                notes = "Toy simulation; no signal between covariates and events."
            )
        )
    )
    return(mae)
}
