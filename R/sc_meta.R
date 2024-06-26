#' scRNAseq meta data
#'
#' The scRNAseq meta data must be in the format of data frame while each row
#' represents a cell. The rownames of the scRNAseq meta data should match
#' exactly with the column names of the scRNAseq count data. The sc_meta data
#' must contain the column indicating the cell type assignment for each cell
#' (e.g., “cellType” column in the example sc_meta data). Sample/subject information
#' should be provided, if there is only one sample, we can add a column by
#' sc_meta$sampleInfo = "sample1".
#'
#'
#' @docType data
#' @keywords datasets
#' @usage data(sc_meta)
#' @name sc_meta
"sc_meta"
