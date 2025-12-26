# R/__init__.R
# Package-level exports for box module system

#' @export
box::use(
  ./domain/study[Study],
  ./domain/effect_size[EffectSize],
  ./calculators/smcr_calculator[SMCRCalculator],
  ./calculators/romc_calculator[ROMCCalculator],
  ./models/frequentist_model[FrequentistModel],
  ./validators/study_validator[StudyValidator],
  ./utils/config[load_config, get_param],
  ./utils/logging[get_logger, setup_logging],
  ./utils/seed[set_reproducible_seed]
)
