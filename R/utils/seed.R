# R/utils/seed.R
# Reproducibility utilities for random seed management

#' @export
box::use(
  ./logging[get_logger]
)

#' Set reproducible seed
#'
#' Sets random seed with logging for reproducibility tracking.
#' This function should be called at the start of any stochastic operation.
#'
#' @param seed Integer seed value
#' @param context Optional context description for logging
#' @return Invisible seed value
#' @export
set_reproducible_seed <- function(seed, context = NULL) {
  log <- get_logger("Seed")

  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
    stop("Seed must be a single numeric value")
  }

  seed <- as.integer(seed)
  set.seed(seed)

  if (!is.null(context)) {
    log$info("Set seed = {seed} for: {context}", seed = seed, context = context)
  } else {
    log$debug("Set seed = {seed}", seed = seed)
  }

  invisible(seed)
}

#' Create a seed generator
#'
#' Creates a function that generates reproducible sub-seeds from a master seed.
#' Useful for generating multiple independent seeds from a single master seed.
#'
#' @param master_seed Master seed value
#' @return Function that generates sequential sub-seeds
#' @export
create_seed_generator <- function(master_seed) {
  counter <- 0L

  function(context = NULL) {
    counter <<- counter + 1L
    sub_seed <- master_seed + counter * 1000L
    set_reproducible_seed(sub_seed, context = context)
    sub_seed
  }
}

#' Get default seed from config or environment
#'
#' Retrieves seed value from configuration, environment variable, or options.
#'
#' @param config Optional configuration object
#' @return Integer seed value
#' @export
get_default_seed <- function(config = NULL) {
  # Priority: config > environment > option > fallback
  seed <- NULL

  if (!is.null(config) && !is.null(config$random$seed)) {
    seed <- config$random$seed
  }

  if (is.null(seed)) {
    env_seed <- Sys.getenv("DEADLIFT_SEED", "")
    if (nzchar(env_seed)) {
      seed <- as.integer(env_seed)
    }
  }

  if (is.null(seed)) {
    seed <- getOption("deadlift.default_seed", 42L)
  }

  as.integer(seed)
}
