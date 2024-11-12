#' The function to estimate the cell type composition signature for each
#' single cell in the scRNaseq reference data
#'
#' @param sc_eset the sc_eset stored in the CARD object
#' @param ct_varname character, the name of the column in metaData that
#' specifies the cell type annotation information, stored in the CARD object
#' @param ct_select vector of cell type names that you are interested in to
#' deconvolute, default as NULL. stored in the CARD object
#' @param sample_varname character,the name of the column in metaData that
#' specifies the sample information. stored in the CARD object
#' @param B reference basis matrix stored in the CARD object.

#' @importFrom SummarizedExperiment assays
#' @importFrom nnls nnls

#' @return Returns a matrix of the cell type composition signature for each
#' single cell in the scRNaseq reference
#'
get_weight_for_cell <- function(
        sc_eset, 
        ct_varname, 
        ct_select, 
        sample_varname, 
        B){
    ##### Calculate cell type composition signature matrix
    count <- assays(sc_eset)$counts
    count <- count[, colSums(count) > 0]
    count <- sweep(count, 2, colSums(count), "/")
    count <- count[rownames(count) %in% rownames(B), ]
    count <- count[match(rownames(B), rownames(count)), ]
    mean_cell <- vapply(seq_len(ncol(count)), function(icell) {
        mod1 <- nnls(as.matrix(B), as.matrix(count[, icell]))
        mod1$x
    }, FUN.VALUE = numeric(ncol(B)))
    mean_cell <- t(mean_cell)
    rownames(mean_cell) <- colnames(count)
    colnames(mean_cell) <- colnames(B)
    return(mean_cell)
    }
#' The function to sample the spatial location information for each single
#' cell
#'
#' @param cords The spatial location information in the measure spatial
#' locations, with the first and second columns represent the 2-D 
#' x-y coordinate system
#' @param numcell a numeric value indicating the number of single cells in each
#' measured location, we suggest 20 for ST technology, 7 for 10x Viisum and 2
#' for Slide-seq
#' @param shape a character indicating whether the sampled spatial coordinates
#' for single cells locating in a Square-like region or a Circle-like region.
#' The center of this region is the measured spatial location in the non-single
#' cell resolution spatial transcriptomics data. The default is 'Square', the
#' other shape is 'Circle'
#' @importFrom spatstat.random runifdisc
#' @importFrom fields rdist
#' @importFrom stats runif median
#' @return Returns a dataframe with the sampled spatial location information
#' for each single cell
#'
get_high_res_cords <- function(cords, numcell, shape = "Square") {
    ##### Return sampled spatial location information
    ed <- rdist(as.matrix(cords))
    n <- dim(cords)[1]
    dis <- c()
    for (i in seq_len(dim(ed)[1])) {
        dis[i] <- min(ed[i, -i])
    }
    min_distance <- median(dis) / 2
    cords_new <- NULL
    for (i in seq_len(nrow(cords))) {
        get_points_within_circle <- function(ed, i, numcell, min_distance) {
            circle <- runifdisc(
                numcell,
                radius = min_distance,
                centre = c(cords[i, 1], cords[i, 2]),
                nsim = 1,
                drop = TRUE
            )
            df <- data.frame(x = circle$x, y = circle$y)
            return(df)
        }
        get_points_within_square <- function(cords, i, numcell, min_distance) {
            minX <- min_distance
            minY <- min_distance
            rectangular_x <- runif(
                numcell,
                min = cords[i, 1] - minX,
                max = cords[i, 1] + minX
            )
            rectangular_y <- runif(
                numcell,
                min = cords[i, 2] - minY,
                max = cords[i, 2] + minY
            )
            df <- data.frame(x = rectangular_x, y = rectangular_y)
            return(df)
        }
        if (shape == "Square") {
            df <- get_points_within_square(cords, i, numcell, min_distance)
        } else if (shape == "Circle") {
            df <- get_points_within_circle(ed, i, numcell, min_distance)
        }
        colnames(df) <- c("x", "y")
        df$centerSPOT <- paste0(cords[i, 1], "x", cords[i, 2])
        df$centerx <- cords[i, 1]
        df$centery <- cords[i, 2]
        cords_new <- rbind(cords_new, df)
    }
    cords_new <- cords_new[!duplicated(paste0(cords_new$x, "x", cords_new$y)), ]
    rownames(cords_new) <- paste0(cords_new$x, "x", cords_new$y)
    return(cords_new)
}

#' The function to assign the spatial location information for
#' each single cell
#'
#' @param mappint_spot_cell_cor a mapped correlation matrix indicating the
#' relashionship between each measured spatial location and the single cell in
#' the scRNAseq reference
#' @param cords_new output from the function get_high_res_cords
#' @param numcell a numeric value indicating the number of single cells in each
#' measured location, we suggest 20 for ST technology, 7 for 10x Viisum and 2
#' for Slide-seq
#' @param sc_eset a single cell experiment object stored in CARD object
#' @param ct_varname character, the name of the column in metaData that
#' specifies the cell type annotation information, stroed in CARD object

#' @importFrom SingleCellExperiment colData
#' @return Return the assigned spatial location information for the mapped
#' single cell
#'
assign_sc_cords <- function(
        mappint_spot_cell_cor, 
        cords_new, 
        numcell,
        sc_eset, 
        ct_varname) {
    ##### Assign spatial information for mapping single cell
    mapcellcords <- NULL
    for (ispot in seq_len(nrow(mappint_spot_cell_cor))) {
        mapcell <- mappint_spot_cell_cor[ispot, ]
        mapcell <- mapcell[order(mapcell, decreasing = TRUE)]
        centerspot <- rownames(mappint_spot_cell_cor)[ispot]
        ispot_cords <- data.frame(
            x = as.numeric(
                vapply(
                    strsplit(centerspot, split = "x"), 
                    `[`, 
                    character(1),
                    1)
                ),
            y = as.numeric(
                vapply(
                    strsplit(centerspot, split = "x"), 
                    `[`, 
                    character(1), 
                    2)
                )
        )
        subcordscell <- cords_new[cords_new$centerSPOT == centerspot, ]
        
        ##### calculate Euclidiean distance with this
        ed_with_center <- rdist(
            as.matrix(subcordscell[, seq_len(2)]),
            as.matrix(ispot_cords)
        )
        subcordscell$ed_with_center <- ed_with_center
        subcordscell <- subcordscell[
            order(subcordscell$ed_with_center,decreasing = FALSE), ]
        map_cell_cords_temp <- data.frame(
            CorwithSpot = mapcell[seq_len(numcell)])
        map_cell_cords_temp <- cbind(map_cell_cords_temp, subcordscell)
        map_cell_cords_temp$CT <- 
            colData(sc_eset)[rownames(map_cell_cords_temp),ct_varname]
        map_cell_cords_temp$Cell <- rownames(map_cell_cords_temp)
        mapcellcords <- rbind(mapcellcords, map_cell_cords_temp)
    }
    return(mapcellcords)
}

#' Extension of CARD into performing single cell Mapping from non-single
#' cell spatial transcriptomics dataset.
#'
#' @param CARD_object CARD object create by the CARD_deconvolution function.
#' @param shapeSpot a character indicating whether the sampled spatial
#' coordinates for single cells locating in a Square-like region or a
#' Circle-like region. The center of this region is the measured spatial
#' location in the non-single cell resolution spatial transcriptomics data.
#' The default is 'Square', the other shape is 'Circle'
#' @param numcell a numeric value indicating the number of single cells in
#' each measured location, we suggest 20 for ST technology, 7 for 10x Viisum
#' and 2 for Slide-seq
#' @param ncore a numeric value indicating the number of cores used to
#' accelerating the procedure

#' @importFrom SummarizedExperiment assays colData
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam
#' @return Returns a SingleCellExperiment SCE object with the mapped expression
#' at single cell resolution and the spatial location information of each
#' single cell
#'
#' @export
#' @examples
#' library(SingleCellExperiment)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
#' scMapping <- CARD_scmapping(
#' CARD_obj, 
#' shapeSpot = "Square", 
#' numcell = 20, 
#' ncore = 2)
#' print(scMapping)
#'
CARD_scmapping <- function(
        CARD_object, 
        shapeSpot = "Square",
        numcell, 
        ncore = 10) {
    #####Check
    if(shapeSpot != 'Square' && shapeSpot != 'Circle'){
        stop("Please use 'Square' or 'Circle'")
    }
    if (!inherits(CARD_object, "SpatialExperiment")) {
        stop("CARD_object must be a SpatialExperiment object is created by
        CARD functions")
    }
    if (!is.numeric(ncore) || 
        length(ncore) != 1 || 
        ncore <= 0 || 
        ncore != as.integer(ncore)) {
        stop("ncore must be a single positive integer")
    }
    if (!is.numeric(numcell) || 
        length(numcell) != 1 || 
        numcell <= 0 || 
        numcell != as.integer(numcell)) {
        stop("numcell must be a single positive integer")
    }
    
    ##### load in spatial transcriptomics data stored in CARD_object
    sc_eset <- CARD_object@metadata$sc_eset
    B <- CARD_object@metadata$algorithm_matrix$B
    ct_select <- CARD_object@metadata$info_parameters$ct.select
    ct_varname <- CARD_object@metadata$info_parameters$ct.varname
    sample_varname <- CARD_object@metadata$info_parameters$sample.varname
    res_CARD <- if (!is.null(CARD_object$Proportion_CARD)) {
        CARD_object$Proportion_CARD
    } else {
        CARD_object@metadata$Proportion_CARD
    }
    res_CARD <- res_CARD[, order(colnames(res_CARD))]
    mean_cell <- get_weight_for_cell(
        sc_eset, 
        ct_varname, 
        ct_select,
        sample_varname, 
        B
    )
    mean_cell <- mean_cell[, order(colnames(mean_cell))]
    cords <- if (dim(spatialCoords(CARD_object))[2] > 0) {
        as.data.frame(spatialCoords(CARD_object)[rownames(res_CARD),])
    } else {
        as.data.frame(
            CARD_object@metadata$spatial_location[rownames(res_CARD), ])
    }
    rownames(res_CARD) <- paste0(cords$x, "x", cords$y)
    cords_new <- get_high_res_cords(
        cords,
        numcell = numcell,
        shape = shapeSpot
    )
    mappint_spot_cell_cor <- cor(t(res_CARD), t(mean_cell))
    rownames(mappint_spot_cell_cor) <- 
        paste0(
            cords$x, 
            "x",
            cords$y
            )
    
    ##### high resolution coordinates
    mapcellcords <- assign_sc_cords(
        mappint_spot_cell_cor, 
        cords_new,
        numcell, 
        sc_eset, 
        ct_varname
    )
    
    ##### map scRNAseq count data
    count_sc <- assays(sc_eset)$counts
    count_sc <- count_sc[, colSums(count_sc) > 0]
    count_ct <- NULL
    
    ##### parallelization parameter
    param <- if (.Platform$OS.type == "windows") {
        SnowParam(workers = ncore)
    } else {
        MulticoreParam(workers = ncore)
    }
    count_ct <- bplapply(seq_len(nrow(res_CARD)), function(ispot) {
        spot <- rownames(res_CARD)[ispot]
        map_cell_cords_spot <- mapcellcords[mapcellcords$centerSPOT == spot, ]
        df <- as(count_sc[, map_cell_cords_spot$Cell], "sparseMatrix")
        colnames(df) <- paste0(map_cell_cords_spot$Cell, ":", spot)
        colnames(df) <- paste0(
            colnames(df), 
            ":", 
            map_cell_cords_spot$x,
            "x", 
            map_cell_cords_spot$y
        )
        df
        }, BPPARAM = param)
    count_ct <- do.call("cbind", count_ct)
    rownames(mapcellcords) <- colnames(count_ct)
    sce <- SingleCellExperiment(
        list(counts = count_ct),
        colData = as.data.frame(mapcellcords[, c(
            "x", 
            "y", 
            "centerSPOT",
            "centerx", 
            "centery",
            "CT", 
            "Cell")]),
        rowData = as.data.frame(rownames(count_ct))
    )
    return(sce)
}
