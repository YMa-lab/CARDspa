.onLoad <- function(libname, pkgname) {
    # CRAN Note avoidance
    utils::globalVariables(
        # variable names for ggplot
        c("x",
          "y",
          "value")
    )
    # verify_NA_ALIAS(0L)
    invisible(NULL)
}
