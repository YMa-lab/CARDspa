# CARDspa News

## Version 0.99.5
- Bioconductor-ready.

## Version 0.99.4
- Updated example code in function documentation.

## Version 0.99.3
- Added a README file.
- Improved unit test coverage.
- Addressed build check notes:
  - Replaced `sapply` with `vapply` where appropriate.
  - Used `seq()` functions instead of colon `:` for sequence generation.
  - Replaced generic variable names with more descriptive ones.
  - Added stop checks to validate input parameters.
  - Removed commented-out code or moved it to a separate section if necessary.
  - Adjusted parallelization approach to ensure compatibility across different operating systems.
  - Moved class definitions to `AllClasses.R` for better organization.
- Improved code readability and consistency
- Improved handling of `SpatialExperiment` object interoperability.

## Version 0.99.2
- Added support for `SingleCellExperiment (sce)` and `SpatialExperiment (spe)` as input objects.
- Implemented `show()` method for better object representation.

## Version 0.99.1
- Removed internal `set.seed()` calls to ensure reproducibility.


## Version 0.99.0
- Initial bioconductor package

