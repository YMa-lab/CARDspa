#' Make new spatial locations on unmeasured tissue through grids.
#'
#' @param location Data frame, spatial location data frame of the original
#' spatial resolved transcriptomics dataset, stored in the
#' spatialCoords(CARD_object)
#' @param num_sample Numeric, approximate number of cells in grid within the
#' shape of the spatial location data frame
#' @param concavity Numeric, a relative measure of concavity. The default is
#' 2.0, which can prodecure detailed enough shapes. Infinity results in a
#' convex hull while 1 results in a more detailed shape.
#'
#' @import concaveman
#' @import sp
#' @importFrom dplyr '%>%'
#' @importFrom sf st_as_sf
#' @return Return a list of data frame with newly grided points
#'
#'

sample_grid_within <- function(location, num_sample, concavity = 2) {
    #### recognise the outline of the shape and extract the coordinates of
    ##### the polygons
    df <- location[, seq_len(2)]
    pnts <- df %>%
        st_as_sf(coords = c("x", "y"))
    polygon <- concaveman(pnts, concavity) #### tune the number, default is 2.0
    poly_coords <- as.data.frame(as.matrix(polygon$polygons[[1]]))
    colnames(poly_coords) <- c("x", "y")
    coordinates(poly_coords) <- ~ x + y
    
    ##### create the spatial polygon data frame
    crds <- poly_coords
    pl <- Polygon(crds)
    id <- "21*15"
    pls <- Polygons(list(pl), ID = id)
    spls <- SpatialPolygons(list(pls))
    df <- data.frame(value = 1, row.names = id)
    spdf <- SpatialPolygonsDataFrame(spls, df)
    
    ##### random make the grids within the polygon
    pts <- makegrid(spdf, num_sample)
    pts1 <- SpatialPoints(pts)
    
    ##### extract the points within the locations
    spgrd_within <- SpatialPixels(pts1[spdf, ])
    data <- as.data.frame(spgrd_within@coords)
    data$id <- seq_len(nrow(data))
    colnames(data) <- c("x", "y")
    return(data)
}

#' Normalize the new spatial locations without changing the shape and
#' relative positions
#'
#' @param location_orig Data frame, spatial location data frame of the original
#' spatial resolved transcriptomics dataset, stored in the
#' spatialCoords(CARD_object)
#' @param train_ind Vector, Index of the original spatial locations
#' @param test_ind Vector, Index of the newly grided spatial locations
#'
#' @return Return the normalized spatial location data frame
#'

norm_coords_train_test <- function(location_orig, train_ind, test_ind) {
    ##### normalize to 0-1 scale
    norm_coords_train <- as.data.frame(location_orig[train_ind, seq_len(2)])
    location_factor_x <- min(norm_coords_train$x)
    location_factor_y <- min(norm_coords_train$y)
    norm_coords_train$x <- norm_coords_train$x - location_factor_x
    norm_coords_train$y <- norm_coords_train$y - location_factor_y
    scale_factor <- max(norm_coords_train$x, norm_coords_train$y)
    norm_coords_train$x <- norm_coords_train$x / scale_factor
    norm_coords_train$y <- norm_coords_train$y / scale_factor
    norm_coords_test <- as.data.frame(location_orig[test_ind, seq_len(2)])
    norm_coords_test[, 1] <- (norm_coords_test[, 1] -
        location_factor_x) / scale_factor
    norm_coords_test[, 2] <- (norm_coords_test[, 2] -
        location_factor_y) / scale_factor
    locationTemp <- rbind(norm_coords_test, norm_coords_train)
    return(locationTemp)
}

#' Calculate the variance covariance matrix used in the imputation of the
#' new grided locations
#'
#' @param location_orig Data frame, spatial location data frame of the original
#' spatial resolved transcriptomics dataset, stored in the
#' spatialCoords(CARD_object)
#' @param train_ind Vector, index of the original spatial locations
#' @param test_ind Vector, index of the newly grided spatial locations
#' @param optimal_phi Numeric, the optimal phi value stored in CARD_object
#' @param ineibor Numeric, number of neighbors used in the imputation on newly
#' grided spatial locations, default is 10.
#'
#' @import Matrix
#' @importFrom RANN nn2
#' @return Return a list with the imputed Cell type composition Vtest matrix on
#' the newly grided spatial locations and predicted normalized gene expression
#'

Sigma <- function(location_orig, train_ind, test_ind, optimal_phi, ineibor) {
    norm_cords_temp <- norm_coords_train_test(
        location_orig, 
        train_ind, 
        test_ind)
    ##### find neiibors
    near_data <- nn2(norm_cords_temp[, seq_len(2)], k = ineibor + 1)
    neibors <- near_data$nn.idx
    neibors <- neibors[, -1] ##### delete the location itself as the neighbor
    nmat <- Matrix(0,
        nrow = nrow(neibors),
        ncol = nrow(neibors), 
        sparse = TRUE
    )
    
    ##### speed up the matrices when it is large
    for (icol in seq_len(ncol(neibors))) {
        edges <- data.frame(
            i = seq_len(nrow(neibors)),
            j = neibors[, icol]
        )
        adjacency <- sparseMatrix(
            i = as.integer(edges$i),
            j = as.integer(edges$j),
            x = 1,
            dims = rep(nrow(neibors), 2),
            use.last.ij = TRUE
        )
        nmat <- nmat + adjacency
    }
    
    ##### find mutual neibors, for non-mutual neighbors, the value will be zero
    nmat <- nmat * t(nmat)
    
    ##### create euclidiean distance
    dist_neibors <- near_data$nn.dists
    isigma <- 0.1
    kernel_neibors <- exp(-dist_neibors^2 / (2 * isigma^2))
    wtemp <- Matrix(0,
        nrow = nrow(neibors),
        ncol = nrow(neibors),
        sparse = TRUE
    )
    
    ##### speed up the matrices when it is large
    for (icol in seq_len(ncol(kernel_neibors))) {
        temp <- near_data$nn.idx
        edges <- data.frame(
            i = seq_len(nrow(temp)),
            j = temp[, icol]
        )
        adjacency <- sparseMatrix(
            i = as.integer(edges$i),
            j = as.integer(edges$j),
            x = kernel_neibors[, icol],
            dims = rep(nrow(neibors), 2),
            use.last.ij = TRUE
        )
        wtemp <- wtemp + adjacency
    }
    
    #### create euclidiean distance for mutual neibors, for non-mutual
    ##### neighbors, the value will be zero
    wtemp <- wtemp * nmat
    rownames(wtemp) <- rownames(norm_cords_temp)
    colnames(wtemp) <- rownames(norm_cords_temp)
    
    ##### calculate variance matrix that I need
    dtemp <- Diagonal(x = colSums(wtemp))
    In <- Diagonal(nrow(wtemp))
    w22 <- wtemp[
        (length(test_ind) + 1):nrow(wtemp),
        (length(test_ind) + 1):ncol(wtemp)
    ]
    sigmatemp <- dtemp - 
        optimal_phi * wtemp
    sigma11 <- sigmatemp[
        seq_len(length(test_ind)),
        seq_len(length(test_ind))
    ]
    sigma12 <- sigmatemp[
        seq_len(length(test_ind)),
        (length(test_ind) + 1):ncol(sigmatemp)
    ]
    sigma21 <- sigmatemp[
        (length(test_ind) + 1):nrow(sigmatemp),
        seq_len(length(test_ind))
    ]
    sigma22 <- sigmatemp[
        (length(test_ind) + 1):nrow(sigmatemp),
        (length(test_ind) + 1):ncol(sigmatemp)
    ]
    return(list(
        SigmaTemp = sigmatemp,
        Sigma11 = sigma11, 
        Sigma12 = sigma12,
        Sigma21 = sigma21, 
        Sigma22 = sigma22,
        W22 = w22
    ))
}

#' Imputation and Construction of High-Resolution Spatial Maps for Cell Type
#' Composition and Gene Expression by the spatial correlation structure between
#' original spatial locations and new grided spatial locations
#'
#' @param vtrain Matrix, estimated V matrix from CARD
#' @param location_orig Data frame, spatial location data frame of the original
#' spatial resolved transcriptomics dataset, stored in the
#' spatialCoords(CARD_object)
#' @param train_ind Vector, index of the original spatial locations
#' @param test_ind Vector, index of the newly grided spatial locations
#' @param B Matrix, used in the deconvolution as the reference basis matrix
#' @param xinput_norm Matrix, used in the deconvolution as the normalized
#' spatial count data
#' @param optimal_b Vector, vector of the intercept for each cel type estimated
#' based on the original spatial resolution
#' @param optimal_phi Numeric, the optimal phi value stored in CARD_object
#' @param lambda Vector, vector of cell type specific scalar in the CAR model
#' @param ineibor Numeric, number of neighbors used in the imputation on newly
#' grided spatial locations, default is 10.
#'
#' @return Return a list with the imputed Cell type composition Vtest matrix on
#' the newly grided spatial locations and predicted normalized gene expression
#'


mvn_cv <- function(
        vtrain, 
        location_orig, 
        train_ind, 
        test_ind, 
        B, 
        xinput_norm,
        optimal_b, 
        optimal_phi, 
        lambda, 
        ineibor) {
    ##### Calculate imputed Cell type composition Vtest matrix
    sigma_list <- Sigma(
        location_orig, 
        train_ind, 
        test_ind, 
        optimal_phi, 
        ineibor)
    sigma11 <- sigma_list$Sigma11
    sigma21 <- sigma_list$Sigma21
    sigma12 <- sigma_list$Sigma12
    sigma22 <- sigma_list$Sigma22
    w22 <- sigma_list$W22
    mu_cond <- NULL
    for (ict in seq_len(length(lambda))) {
        temp <- 
            sigma12 %*% 
            (vtrain[, ict] - optimal_b[seq_len(nrow(vtrain)), ict])
        solver <- solve(sigma11, temp)
        Mu <- optimal_b[seq_len(length(test_ind)), ict] - solver ### q * k
        mu_cond <- cbind(mu_cond, Mu)
    }
    vtest <- mu_cond
    colnames(vtest) <- colnames(vtrain)
    
    ##### Calculate the predicted normalized gene expression
    x_test_hat <- B %*% t(vtest)
    rownames(vtest) <- colnames(x_test_hat)
    vtest <- sweep(vtest, 1, rowSums(vtest), "/")
    return(list(Vtest = vtest, XtestHat = x_test_hat))
}

#' Construct an enhanced spatial expression map on the unmeasured tissue
#' locations
#'
#' @param CARD_object SpatialExperiment Object created by 
#' CARD_deconvolution with estimated cell type compositions on the
#' original spatial resolved transcriptomics data.
#' @param num_grids Initial number of newly grided spatial locations. The final
#' number of newly grided spatial locations will be lower than this value since
#' the newly grided locations outside the shape of the tissue will be filtered
#' @param ineibor Numeric, number of neighbors used in the imputation on newly
#' grided spatial locations, default is 10.
#' @param exclude Vector, the rownames of spatial location data on the original
#' resolution that you want to exclude. This is to avoid the weird detection of
#' the shape.
#'
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#'
#' @return Return a SpatialExperiment object with the refined cell type 
#' compositions estimated for newly grided spots and the refined predicted 
#' gene expression (normalized).
#'
#' @export
#' @examples
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
#' CARD_obj <- CARD_imputation(
#'     CARD_obj,
#'     num_grids = 200,
#'     ineibor = 10,
#'     exclude = NULL
#' )
#'
CARD_imputation <- function(
        CARD_object, 
        num_grids, 
        ineibor = 10, 
        exclude = NULL){
    ##### Check inputs
    if (!inherits(CARD_object, "SpatialExperiment")) {
        stop("CARD_object must be a SpatialExperiment object is created by
        CARD functions")
    }
    if (!is.numeric(num_grids) || 
        length(num_grids) != 1 || 
        num_grids <= 0 || 
        num_grids != as.integer(num_grids)) {
        stop("num_grids must be a single positive integer")
    }
    if (!is.numeric(ineibor) ||
        length(ineibor) != 1 ||
        ineibor <= 0 ||
        ineibor != as.integer(ineibor)) {
        stop("ineibor must be a single positive integer")
    }
    if (!is.null(exclude) && !is.character(exclude)) {
        stop("exclude must be NULL or a character vector")
    }
    
    ##### extract results from CARD object
    B <- CARD_object@metadata$algorithm_matrix$B
    xinput_norm <- CARD_object@metadata$algorithm_matrix$Xinput_norm
    vtrain <- CARD_object@metadata$algorithm_matrix$Res$V
    location <- as.data.frame(spatialCoords(CARD_object))
    
    ##### check
    if (sum(rownames(location) == rownames(vtrain)) == nrow(vtrain)) {
        message("## The rownames of locations are matched ...\n")
    }
    
    ##### Make the new grids
    message("## Make grids on new spatial locations ...\n")
    if (!is.null(exclude)) {
        location_use_to_sample <- 
            location[!(rownames(location) %in% exclude), ]
    } else {
        location_use_to_sample <- location
    }
    data <- sample_grid_within(location_use_to_sample, num_grids, concavity = 2)
    rownames(data) <- paste0(data$x, "x", data$y)
    
    ##### delete the newly grided locations that are the same as the original
    ##### ones
    data <- data[!(rownames(data) %in% rownames(location)), ]
    location_train <- location[, seq_len(2)]
    location_test <- data[, seq_len(2)]
    location_combine <- rbind(data[, seq_len(2)], location[, seq_len(2)])
    test_ind <- seq_len(nrow(data))
    train_ind <- (nrow(data) + 1):nrow(location_combine)
    location_orig <- location_combine
    
    ##### results from the CARD
    optimal_b <- CARD_object@metadata$algorithm_matrix$Res$b
    vec_one <- rep(1, nrow(location_combine))
    optimal_b_sum <- vec_one %*% t(optimal_b)
    optimal_phi <- CARD_object@metadata$info_parameters$phi
    lambda <- CARD_object@metadata$algorithm_matrix$Res$lambda
    imputation <- 
        mvn_cv(
            vtrain, 
            location_orig, 
            train_ind, 
            test_ind,
            B, 
            xinput_norm, 
            optimal_b_sum, 
            optimal_phi, 
            lambda, 
            ineibor
        )
    spe_out <- SpatialExperiment(
        assay = list(refined_expression = as.matrix(imputation$XtestHat)))
    colData(spe_out)$refined_prop <- as.matrix(imputation$Vtest)
    
    ## remove sample_id 
    spe_out@colData@listData[["sample_id"]] <- NULL
    metadata <- metadata(CARD_object)
    metadata(spe_out) <- metadata
    metadata(spe_out)$Proportion_CARD <- CARD_object$Proportion_CARD
    metadata(spe_out)$spatial_countMat <- 
        CARD_object@assays@data$spatial_countMat
    metadata(spe_out)$spatial_location <- spatialCoords(CARD_object)
    
    return(spe_out)
}
