# user-inputs.R
# User-facing data specification and canonical input adapter.
# =============================================================================

#' Describe the columns and derived features used by SILK
#'
#' SILK operates internally on a small canonical schema, but users do not need
#' to rename their data.  This constructor records the mapping from a user's
#' subject/visit tables to that schema and, optionally, named functions for
#' engineered covariates.  Comparator-specific model formulas do not belong in
#' this object; they are application or simulation code.
#'
#' @param subject_id Subject identifier column in both tables.
#' @param recorded_age Recorded landmark-age column in the subject table.
#' @param event_time Follow-up time column in the subject table.
#' @param event Event indicator column in the subject table.
#' @param visit_number Visit-position column in the visit table.  If NULL,
#'   positions are assigned within subject after sorting by visit_age.
#' @param visit_age Recorded age at a visit.
#' @param visit_lag Time before the recorded landmark.  If NULL, it is derived
#'   from recorded_age and visit_age.
#' @param biomarker_cols Character vector of visit-table biomarker columns.
#'   NULL preserves the package convention of detecting B1, B2, ... columns.
#' @param covariate_cols Character vector of subject-level covariate columns
#'   used by the SILK survival layer.
#' @param engineered_covariates Named list of functions.  Each function takes a
#'   data frame and returns one vector; the resulting named column is added to
#'   both tables (subject-level values are propagated to visits).
#' @param template_input_covariates Character vector used by the template-input
#'   kernel.  If omitted with covariate_cols supplied, covariates plus lag are
#'   used.
#' @param template_input_bandwidths Optional named positive numeric vector.
#' @param clock_covariates Character vector used by the longitudinal clock
#'   residualizer.  Defaults to covariate_cols when supplied.
#' @param anchor_covariates Subject-level covariates used by the shift anchor.
#' @param anchor_mode Either "intercept" or "rx".  A user-created
#'   specification defaults to an intercept-only anchor unless columns are
#'   explicitly supplied.
#' @return An object of class code{silk_data_spec}.
#' @export
silk_data_spec <- function(
    subject_id = "id", recorded_age = "A_obs", event_time = "U",
    event = "delta", visit_number = "visit", visit_age = "A_obs_il",
    visit_lag = "lag", biomarker_cols = NULL, covariate_cols = NULL,
    engineered_covariates = NULL, template_input_covariates = NULL,
    template_input_bandwidths = NULL, clock_covariates = NULL,
    anchor_covariates = NULL, anchor_mode = NULL) {
  scalar_character <- function(x, name, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(NULL)
    if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
      stop(name, " must be one non-empty character string or NULL.", call. = FALSE)
    }
    x
  }
  character_vector <- function(x, name, allow_null = TRUE) {
    if (is.null(x) && allow_null) return(NULL)
    if (!is.character(x) || any(!nzchar(x)) || anyDuplicated(x)) {
      stop(name, " must be a character vector with unique non-empty names.", call. = FALSE)
    }
    unique(x)
  }

  subject_id <- scalar_character(subject_id, "subject_id")
  recorded_age <- scalar_character(recorded_age, "recorded_age")
  event_time <- scalar_character(event_time, "event_time", allow_null = TRUE)
  event <- scalar_character(event, "event", allow_null = TRUE)
  visit_number <- scalar_character(visit_number, "visit_number", allow_null = TRUE)
  visit_age <- scalar_character(visit_age, "visit_age")
  visit_lag <- scalar_character(visit_lag, "visit_lag", allow_null = TRUE)
  biomarker_cols <- character_vector(biomarker_cols, "biomarker_cols")
  covariate_cols <- character_vector(covariate_cols, "covariate_cols")
  template_input_covariates <- character_vector(
    template_input_covariates, "template_input_covariates"
  )
  clock_covariates <- character_vector(clock_covariates, "clock_covariates")
  anchor_covariates <- character_vector(anchor_covariates, "anchor_covariates")

  if (is.null(engineered_covariates)) engineered_covariates <- list()
  if (!is.list(engineered_covariates) ||
      (length(engineered_covariates) &&
       (is.null(names(engineered_covariates)) ||
        any(!nzchar(names(engineered_covariates))) ||
        anyDuplicated(names(engineered_covariates))))) {
    stop("engineered_covariates must be a named list of functions.", call. = FALSE)
  }
  if (length(engineered_covariates) &&
      !all(vapply(engineered_covariates, is.function, logical(1)))) {
    stop("Each engineered_covariates entry must be a function.", call. = FALSE)
  }

  # A user-created specification is self-contained: omitted survival
  # covariates means an age-only survival layer, not an implicit dependency on
  # the simulation's X1--X4 columns.  The visit lag remains a useful default
  # input for the template kernel.
  if (is.null(covariate_cols)) covariate_cols <- character(0)
  if (is.null(template_input_covariates)) {
    template_input_covariates <- unique(c(covariate_cols, "lag"))
  }
  if (is.null(clock_covariates)) {
    clock_covariates <- covariate_cols
  }
  if (is.null(anchor_mode)) anchor_mode <- if (length(anchor_covariates)) "rx" else "intercept"
  anchor_mode <- match.arg(anchor_mode, c("intercept", "rx"))
  if (identical(anchor_mode, "rx") && !length(anchor_covariates)) {
    stop("anchor_mode='rx' requires anchor_covariates.", call. = FALSE)
  }

  if (!is.null(template_input_bandwidths)) {
    if (is.null(names(template_input_bandwidths)) ||
        any(!nzchar(names(template_input_bandwidths))) ||
        any(!is.finite(as.numeric(template_input_bandwidths))) ||
        any(as.numeric(template_input_bandwidths) <= 0)) {
      stop("template_input_bandwidths must be a named positive numeric vector.", call. = FALSE)
    }
    bw_names <- names(template_input_bandwidths)
    template_input_bandwidths <- as.numeric(template_input_bandwidths)
    names(template_input_bandwidths) <- bw_names
  }

  structure(
    list(
      subject_id = subject_id,
      recorded_age = recorded_age,
      event_time = event_time,
      event = event,
      visit_number = visit_number,
      visit_age = visit_age,
      visit_lag = visit_lag,
      biomarker_cols = biomarker_cols,
      covariate_cols = covariate_cols,
      engineered_covariates = engineered_covariates,
      template_input_covariates = template_input_covariates,
      template_input_bandwidths = template_input_bandwidths,
      clock_covariates = clock_covariates,
      anchor_covariates = anchor_covariates,
      anchor_mode = anchor_mode
    ),
    class = "silk_data_spec"
  )
}

#' @keywords internal
validate_silk_data_spec <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (inherits(spec, "silk_data_spec")) return(spec)
  if (is.list(spec)) return(do.call(silk_data_spec, spec))
  stop("data_spec must be created by silk_data_spec().", call. = FALSE)
}

#' Adapt user data to SILK's internal schema
#'
#' @param subjects Subject-level data frame.
#' @param visits Visit-level data frame.
#' @param spec A code{silk_data_spec} object.
#' @return A list with normalized code{subjects}, code{visits}, and code{spec}.
#' @export
prepare_silk_inputs <- function(subjects, visits, spec = NULL) {
  spec <- validate_silk_data_spec(spec)
  if (is.null(spec)) {
    return(list(subjects = subjects, visits = visits, spec = NULL))
  }
  if (!is.data.frame(subjects) || !is.data.frame(visits)) {
    stop("subjects and visits must both be data frames.", call. = FALSE)
  }
  subject_value <- function(data, source, canonical, required = TRUE) {
    if (is.null(source)) {
      if (required) stop("A source column is required for ", canonical, ".", call. = FALSE)
      return(data)
    }
    if (!source %in% names(data)) {
      if (required) stop("Missing column '", source, "' required for ", canonical, ".", call. = FALSE)
      return(data)
    }
    data[[canonical]] <- data[[source]]
    data
  }
  subjects <- subject_value(subjects, spec$subject_id, "id")
  subjects <- subject_value(subjects, spec$recorded_age, "A_obs")
  subjects <- subject_value(subjects, spec$event_time, "U", required = FALSE)
  subjects <- subject_value(subjects, spec$event, "delta", required = FALSE)
  if (anyDuplicated(as.character(subjects$id))) {
    stop("subject_id must identify exactly one row per subject.", call. = FALSE)
  }

  visits <- subject_value(visits, spec$subject_id, "id")
  visits <- subject_value(visits, spec$visit_age, "A_obs_il")
  if (!is.null(spec$visit_number)) {
    visits <- subject_value(visits, spec$visit_number, "visit")
  } else {
    order_index <- order(visits$id, visits$A_obs_il, seq_len(nrow(visits)))
    sorted <- visits[order_index, , drop = FALSE]
    sorted$visit <- ave(seq_len(nrow(sorted)), sorted$id, FUN = seq_along)
    visits$visit <- sorted$visit[match(seq_len(nrow(visits)), order_index)]
  }
  if (!is.null(spec$visit_lag) && spec$visit_lag %in% names(visits)) {
    visits$lag <- visits[[spec$visit_lag]]
  } else {
    landmark <- stats::setNames(as.numeric(subjects$A_obs), as.character(subjects$id))
    visits$lag <- as.numeric(landmark[as.character(visits$id)]) - as.numeric(visits$A_obs_il)
  }
  if (any(!is.finite(as.numeric(visits$lag)))) {
    stop("visit_lag could not be computed for every visit row.", call. = FALSE)
  }

  apply_engineered <- function(data) {
    if (!length(spec$engineered_covariates)) return(data)
    for (nm in names(spec$engineered_covariates)) {
      value <- tryCatch(
        spec$engineered_covariates[[nm]](data),
        error = function(e) stop("Engineered covariate '", nm, "' failed: ", conditionMessage(e), call. = FALSE)
      )
      if (length(value) == 1L) value <- rep(value, nrow(data))
      if (length(value) != nrow(data)) {
        stop("Engineered covariate '", nm, "' must return one value per row.", call. = FALSE)
      }
      data[[nm]] <- value
    }
    data
  }
  subjects <- apply_engineered(subjects)

  # Raw baseline covariates are propagated before evaluating engineered visit
  # features, so a function such as function(d) d$X3^2 - 1 works on both
  # tables.  The remaining derived columns are propagated after evaluation.
  propagate_columns <- function(data, columns) {
    for (nm in unique(columns)) {
      if (!nzchar(nm) || nm %in% names(data) || !nm %in% names(subjects)) next
      data[[nm]] <- subjects[[nm]][match(as.character(data$id), as.character(subjects$id))]
    }
    data
  }
  visits <- propagate_columns(visits, c(
    spec$covariate_cols, spec$clock_covariates, spec$template_input_covariates,
    spec$anchor_covariates,
    setdiff(names(subjects), c("id", "A_obs", "U", "delta"))
  ))
  visits <- apply_engineered(visits)

  # Baseline covariates and engineered subject features are available at every
  # visit.  Copy them by id only when the visit table does not already carry a
  # column with that name.
  propagate <- unique(c(spec$covariate_cols, spec$clock_covariates,
                        spec$template_input_covariates, spec$anchor_covariates,
                        names(spec$engineered_covariates)))
  for (nm in propagate) {
    if (!nzchar(nm) || nm %in% names(visits) || !nm %in% names(subjects)) next
    visits[[nm]] <- subjects[[nm]][match(as.character(visits$id), as.character(subjects$id))]
  }

  requested_biomarkers <- spec$biomarker_cols
  if (!is.null(requested_biomarkers)) {
    missing <- setdiff(requested_biomarkers, names(visits))
    if (length(missing)) stop("Missing biomarker columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!all(c("id", "visit", "A_obs_il", "lag") %in% names(visits))) {
    stop("The normalized visit table is missing id, visit, A_obs_il, or lag.", call. = FALSE)
  }
  list(subjects = subjects, visits = visits, spec = spec)
}
