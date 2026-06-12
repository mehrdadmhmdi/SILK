# =============================================================================
# zzz.R
# Package load hook
# =============================================================================

.onLoad <- function(libname, pkgname) {
  defaults <- silk_default_options()
  for (nm in names(defaults)) {
    opt_name <- paste0("silk.", nm)
    if (is.null(getOption(opt_name))) {
      options(stats::setNames(list(defaults[[nm]]), opt_name))
    }
  }
}
