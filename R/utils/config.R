# R/utils/config.R
# Configuration management utilities

#' @export
box::use(
  config[get],
  yaml[read_yaml],
  ./logging[get_logger]
)

#' Load Configuration
#'
#' Loads configuration from YAML file with environment override.
#'
#' @param config_file Path to config file (default: "config/default.yml")
#' @param environment Environment name (default: from R_CONFIG_ACTIVE env var)
#' @return Validated configuration list
#' @export
load_config <- function(config_file = "config/default.yml",
                        environment = Sys.getenv("R_CONFIG_ACTIVE", "default")) {
  log <- get_logger("Config")

  if (!file.exists(config_file)) {
    stop(sprintf("Config file not found: %s", config_file))
  }

  log$info("Loading config from {config_file} [env: {environment}]",
           config_file = config_file,
           environment = environment)

  cfg <- config::get(file = config_file, config = environment)

  # Validate required sections

  required <- c("effect_size", "models", "validation")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop(sprintf("Missing required config sections: %s", paste(missing, collapse = ", ")))
  }

  log$debug("Configuration loaded successfully")
  cfg
}

#' Get nested configuration parameter
#'
#' Safely retrieves a nested parameter from configuration.
#'
#' @param cfg Configuration object
#' @param ... Path to parameter (e.g., "models", "frequentist", "method")
#' @param default Default value if not found (default: NULL, which throws error)
#' @return Parameter value
#' @export
get_param <- function(cfg, ..., default = NULL) {
  path <- c(...)
  value <- cfg

  for (key in path) {
    if (is.null(value[[key]])) {
      if (is.null(default)) {
        stop(sprintf("Config key not found: %s", paste(path, collapse = ".")))
      }
      return(default)
    }
    value <- value[[key]]
  }

  value
}

#' Merge configurations
#'
#' Deep merge of two configuration lists, with override taking precedence.
#'
#' @param base Base configuration
#' @param override Override configuration
#' @return Merged configuration
merge_config <- function(base, override) {
  if (!is.list(base) || !is.list(override)) {
    return(override)
  }

  for (name in names(override)) {
    if (name %in% names(base) && is.list(base[[name]]) && is.list(override[[name]])) {
      base[[name]] <- merge_config(base[[name]], override[[name]])
    } else {
      base[[name]] <- override[[name]]
    }
  }

  base
}
