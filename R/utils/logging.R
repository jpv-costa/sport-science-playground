# R/utils/logging.R
# Simplified logging utilities using lgr package

box::use(
  R6[R6Class],
  lgr[lgr, get_logger]
)

# Singleton storage
.state <- new.env(parent = emptyenv())
.state$configured <- FALSE

#' Setup logging with configuration
#'
#' @param level Log level (trace, debug, info, warn, error, fatal)
#' @param log_file Optional log file path
#' @return Invisible NULL
#' @export
setup_logging <- function(level = "info", log_file = NULL) {
  if (.state$configured) return(invisible(NULL))

  lgr$set_threshold(level)

  if (!is.null(log_file)) {
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE)
    }
    lgr$add_appender(lgr::AppenderFile$new(file = log_file), name = "file")
  }

  .state$configured <- TRUE
  invisible(NULL)
}

#' Get a logger by name
#'
#' @param name Logger name (typically class or module name)
#' @return lgr Logger instance
#' @export
get_logger <- function(name) {
  lgr::get_logger(name)
}
