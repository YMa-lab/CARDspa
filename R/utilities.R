#' Quality control of scRNA-seq count data
#'
#' @param counts_in Raw scRNAseq count data, each column is a cell and each
#' row is a gene.
#' @param metadata data frame, metadata with "ct_varname" specify the cell type
#' annotation information and "sample_varname" specify the sample information
#' @param ct_varname character, the name of the column in metadata that
#' specifies the cell type annotation information
#' @param ct_select vector of cell type names that you are interested in to
#' deconvolute, default as NULL. If NULL, then use all cell types provided by
#' single cell dataset;
#' @param sample_varname character,the name of the column in metadata that
#' specifies the sample information. If NULL, we just use the whole as one
#' sample.
#' @param min_cells numeric, we filtered out the non-expressed cells.
#' @param min_genes numeric we filtered out the non-expressed genes
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return Return the filtered scRNA-seq data and meta data stored in a S4 class
#' (SingleCellExperiment)
#'

sc_QC <- function(
        counts_in, 
        metadata, 
        ct_varname, 
        ct_select,
        sample_varname = NULL, 
        min_cells = 0, 
        min_genes = 0) {
    ##### Filter based on min.features
    coldf <- metadata
    counts <- counts_in
    if (min_genes >= 0) {
        nfeatures <- colSums(x = counts)
        counts <- counts[, which(x = nfeatures > min_genes)]
        coldf <- coldf[which(x = nfeatures > min_genes), ]
    }
    
    ##### filter genes on the number of cells expressing
    if (min_cells >= 0) {
        num.cells <- rowSums(x = counts > 0)
        counts <- counts[which(x = num.cells > min_cells), ]
    }
    fdata <- as.data.frame(rownames(counts))
    rownames(fdata) <- rownames(counts)
    keep_cell <- as.character(coldf[, ct_varname]) %in% ct_select
    counts <- counts[, keep_cell]
    coldf <- coldf[keep_cell, ]
    keep_gene <- rowSums(counts) > 0
    fdata <- as.data.frame(fdata[keep_gene, ])
    counts <- counts[keep_gene, ]
    sce <- SingleCellExperiment(
        list(counts = counts),
        colData = as.data.frame(coldf),
        rowData = as.data.frame(fdata)
    )
    return(sce)
}

#' Create the CARD object
#'
#' @param sc_count Raw scRNA-seq count data, each column is a cell and each row
#' is a gene.
#' @param sc_meta data frame, with each row representing the cell type and/or
#' sample information of a specific cell. The row names of this data frame
#' should match exactly with the column names of the sc_count data
#' @param spatial_count Raw spatial resolved transcriptomics data, each column
#' is a spatial location, and each row is a gene.
#' @param spatial_location data frame, with two columns representing the x and
#' y coordinates of the spatial location. The rownames of this data frame
#' should match eaxctly with the columns of the spatial_count.
#' @param sce a \code{SingleCellExperiment} object containing scRNA-seq count 
#' data in the \code{counts} assay, and cell types and sample information in 
#' the colData.
#' @param spe a \code{SpatialExperiment} object containing spatial data in 
#' the \code{counts} assay, and spatial coordinates in the spatialCoords.
#' @param ct_varname character, the name of the column in metadata that
#' specifies the cell type annotation information
#' @param ct_select vector of cell type names that you are interested in to
#' deconvolute, default as NULL. If NULL, then use all cell types provided by
#' single cell dataset;
#' @param sample_varname character,the name of the column in metadata that
#' specifies the sample information. If NULL, we just use the whole as one
#' sample.
#' @param mincountgene Minimum counts for each gene
#' @param mincountspot Minimum counts for each spatial location
#'
#' @importFrom SummarizedExperiment assays
#' @import methods
#' @import SingleCellExperiment
#' @import SpatialExperiment
#' @return Returns CARD object with filtered spatial count and single cell
#' RNA-seq dataset.
#'
createCARDObject <- function(
        sc_count, 
        sc_meta, 
        spatial_count, 
        spatial_location,
        ct_varname, 
        ct_select, 
        sample_varname,
        mincountgene = 100, 
        mincountspot = 5,
        sce = NULL, 
        spe = NULL) {
    ##### Extract count and spatial information
    if(is.null(sc_count) || 
        is.null(sc_meta) || 
        is.null(spatial_count) || 
        is.null(spatial_location)){
        if(!is.null(sce) && !is.null(spe)){
            sc_count <- assays(sce)[['counts']]
            sc_meta <- as.data.frame(sce@colData)
            spatial_count <- assays(spe)[['counts']]
            spatial_location <- as.data.frame(spatialCoords(spe))
        } else{
            stop("Please provide SingleCellExperiment object and 
            SpatialExperiment object")
        }
    }
    
    ##### QC on scRNASeq dataset
    message("## QC on scRNASeq dataset! ...\n")
    if (is(sc_count, "matrix")) {
        sc_countmat <- as(as.matrix(sc_count), "sparseMatrix")
    } else if (is(sc_count, "vector")) {
        sc_countmat <- as(t(as.matrix(sc_count)), "sparseMatrix")
    } else if (is(sc_count, "sparseMatrix")) {
        sc_countmat <- sc_count
    } else {
        stop("scRNASeq counts has to be of following forms: vector,
        matrix or sparseMatrix")
    }
    if (missing(x = sc_countmat)) {
        stop("Please provide scRNASeq count data")
    } else if (is.null(sample_varname) || missing(sample_varname)) {
        sample_varname <- "sampleID"
        sc_meta <- as.data.frame(sc_meta)
        sc_meta$sampleID <- "Sample"
    } else if (any(rownames(x = sc_countmat) == "")) {
        stop("Feature names of sc_count matrix cannot be empty",
            call. = FALSE
        )
    } else if (sum(rownames(sc_meta) == colnames(sc_countmat)) !=
            ncol(sc_countmat)) {
        stop("Cell name in scRNAseq count data does not match with
        the rownames of metadata")
    } else if (ncol(sc_countmat) != nrow(sc_meta)) {
        stop("The number of cells in scRNA-seq counts and sc_meta
            should be consistent!
            (sc_count -- p x c; sc_meta -- c x 2)")
    }
    if (is.null(ct_varname)) {
        stop("Please provide the column name indicating the cell type 
        information in the meta data of scRNA-seq")
    } else if (is.null(ct_select)) {
        message("No cell types selected, we will use all the cell types in the
                scRNA-seq data\n")
        ct_select <- unique(sc_meta[, ct_varname])
    }
    ct_select <- as.character(ct_select[!is.na(ct_select)])
    sc_eset <- sc_QC(
        sc_countmat, 
        sc_meta, 
        ct_varname, 
        ct_select, 
        sample_varname)
    
    ##### Check the spatial count dataset
    #### QC on spatial dataset
    message("## QC on spatially-resolved dataset! ...\n")
    if (is(spatial_count, "matrix")) {
        spatial_countmat <- as(as.matrix(spatial_count), "sparseMatrix")
    } else if (is(spatial_count, "vector")) {
        spatial_countmat <- as(t(as.matrix(spatial_count)), "sparseMatrix")
    } else if (is(spatial_count, "sparseMatrix")) {
        spatial_countmat <- spatial_count
    } else {
        stop("spatial resolved transcriptomic counts has to be of
                    following forms: vector,matrix or sparseMatrix")
    }
    if (any(rownames(x = spatial_countmat) == "")) {
        stop("Gene names of spatial count matrix cannot be empty",
            call. = FALSE
        )
    }
    commongene <- intersect(
        rownames(spatial_countmat),
        rownames(assays(sc_eset)$counts)
    )
    if (length(commongene) == 0) {
        stop("There are no common gene names in spatial count data and
                single cell RNAseq count data", call. = FALSE)
    }
    if (is.null(spatial_location)) {
        stop("Please provide the matched spatial location data frame")
    }
    if (ncol(spatial_countmat) != nrow(spatial_location)) {
        stop("The number of spatial locations in spatial_count and 
        spatial_location should be consistent!
        (spatial_count -- p x n; spatial_location -- n x 2)")
    } # end fi
    ## check data order should consistent
    if (!identical(colnames(spatial_countmat), rownames(spatial_location))) {
        stop("The column names of spatial_count and row names of 
        spatial_location should be should be matched each other!
        (spatial_count -- p x n; spatial_location -- n x 2)")
    } # end fi
    
    ##### QC on spatial dataset
    spatial_countmat <- spatial_countmat[rowSums(spatial_countmat > 0) >
        mincountspot, ]
    spatial_countmat <- 
        spatial_countmat[, (colSums(spatial_countmat) >= mincountgene &
                                colSums(spatial_countmat) <= 1e6)]
    spatial_location <- spatial_location[rownames(spatial_location) %in%
        colnames(spatial_countmat), ]
    spatial_location <- spatial_location[
        match(colnames(spatial_countmat),rownames(spatial_location)), ]
    object <- new(
        Class = "CARD",
        sc_eset = sc_eset,
        spatial_countMat = spatial_countmat,
        spatial_location = spatial_location,
        project = "Deconvolution",
        info_parameters = list(
            ct.varname = ct_varname, 
            ct.select = ct_select,
            sample.varname = sample_varname
        )
    )
    return(object)
}

#' Create the CARD object
#'
#' @param markerlist a list of marker genes, with each element of the list
#' being the vector of cell type specific marker genes
#' @param spatial_count Raw spatial resolved transcriptomics data, each column
#' is a spatial location, and each row is a gene.
#' @param spatial_location data frame, with two columns representing the x and
#' y coordinates of the spatial location. The rownames of this data frame
#' should match eaxctly with the columns of the spatial_count.
#' @param spe a \code{SpatialExperiment} object containing spatial data in 
#' the \code{counts} assay, and spatial coordinates in the spatialCoords.
#' @param mincountgene Minimum counts for each gene
#' @param mincountspot Minimum counts for each spatial location
#'
#' @import methods
#' @import SpatialExperiment
#' @return Returns CARDfree object with filtered spatial count and marker gene
#' list.
#'

createCARDfreeObject <- function(
        markerlist, 
        spatial_count, 
        spatial_location,
        mincountgene = 100, 
        mincountspot = 5,
        spe = NULL) {
    ##### Extract count and spatial information
    if(is.null(spatial_count) || is.null(spatial_location)){
        if(!is.null(spe)){
            spatial_count <- assays(spe)[['counts']]
            spatial_location <- as.data.frame(spatialCoords(spe))
        }else{
            stop("Please provide SpatialExperiment object")
        }
    }
    if (is(spatial_count, "matrix")) {
        spatial_countmat <- as(as.matrix(spatial_count), "sparseMatrix")
    } else if (is(spatial_count, "vector")) {
        spatial_countmat <- as(t(as.matrix(spatial_count)), "sparseMatrix")
    } else if (is(spatial_count, "sparseMatrix")) {
        spatial_countmat <- spatial_count
    } else {
        stop("spatial resolved transcriptomic counts has to be of
        following forms: vector,matrix or sparseMatrix")
    } # end fi
    if (any(rownames(x = spatial_countmat) == "")) {
        stop("Gene names of spatial count matrix cannot be empty",
            call. = FALSE
        )
    } # end fi
    if (is.null(spatial_location)) {
        stop("Please provide the matched spatial location data frame")
    } # end fi
    if (ncol(spatial_countmat) != nrow(spatial_location)) {
        stop("The number of spatial locations in spatial_count and 
        spatial_location should be consistent!
        (spatial_count -- p x n; spatial_location -- n x 2)")
    } # end fi
    ## check data order should consistent
    if (!identical(colnames(spatial_countmat), rownames(spatial_location))) {
        stop("The column names of spatial_count and row names of 
        spatial_location should be should be matched each other!
        (spatial_count -- p x n; spatial_location -- n x 2)")
    } # end fi
    
    ##### QC on spatial dataset
    spatial_countmat <- spatial_countmat[rowSums(spatial_countmat > 0) >
        mincountspot, ]
    spatial_countmat <- spatial_countmat[, (colSums(spatial_countmat) >= 
                                                mincountgene &
        colSums(spatial_countmat) <= 1e6)]
    spatial_location <- spatial_location[rownames(spatial_location) %in%
        colnames(spatial_countmat), ]
    spatial_location <- spatial_location[match(
        colnames(spatial_countmat),
        rownames(spatial_location)
    ), ]
    #### check marker gene list
    if (missing(x = markerlist)) {
        stop("Please provide the marker list for CARDfree")
    }
    object <- new(
        Class = "CARDfree",
        spatial_countMat = spatial_countmat,
        spatial_location = spatial_location,
        project = "Deconvolution (reference-free)",
        markerList = markerlist,
        info_parameters = list()
    )
    return(object)
}


