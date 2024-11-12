#' Each CARD object has a number of slots which store information. Key
#' slots to access are listed below.
#'
#' @slot sc_eset The filtered scRNA-seq data along with meta data stored
#' in the format of SingleCellExperiment.
#' @slot spatial_countMat The filtered spatial count data.
#' @slot spatial_location The weights for combining p-values from multiple
#' kernels.
#' @slot Proportion_CARD The estimated cell type proportion by
#' CARD with each row is a spatial location and each column is a
#' cell type.
#' @slot project The name of the project, default is deconvolution.
#' @slot info_parameters The paramters that are used in model fitting.
#' @slot algorithm_matrix The intermediate matrices that are used in
#' the model fitting step.
#' @slot refined_prop The refined cell type proportion matrix estimated
#' by CARD for the newly grided spatial locations. The number of
#' initial grids are defined by the user.
#' @slot refined_expression The refined predicted expression matrix
#' (normalized) estimated by CARD for the newly grided spatial locations.
#' The number of initial grids are defined by the user.
#'
#' @return Return an object of CARD class
setClass("CARD",
        slots = list(
            sc_eset = "ANY",
            spatial_countMat = "ANY",
            spatial_location = "data.frame",
            Proportion_CARD = "matrix",
            project = "character",
            info_parameters = "list",
            algorithm_matrix = "list",
            refined_prop = "matrix",
            refined_expression = "matrix"
        ))


#' Each CARDfree object has a number of slots which store information.
#' Key slots to access are listed below.
#'
#' @slot spatial_countMat The filtered spatial count data.
#' @slot spatial_location The weights for combining p-values from multiple
#' kernels.
#' @slot Proportion_CARD The estimated cell type proportion by CARD
#' with each row is a spatial location and each column is a cell type.
#' @slot estimated_refMatrix The estimated reference matrix by CARDfree
#' with each row represents a gene and each column represents a cell type
#' cluster.
#' @slot project The name of the project, default is deconvolution.
#' @slot markerList The nlist of cell type specific markers, with
#' each element represents the vector of cell type specific markers
#' @slot info_parameters The paramters that are used in model fitting.
#' @slot algorithm_matrix The intermediate matrices that are used in
#' the model fitting step.
#' @slot refined_prop The refined cell type proportion matrix estimated
#' by CARD for the newly grided spatial locations. The number of initial
#' grids are defined by the user.
#' @slot refined_expression The refined predicted expression matrix
#' (normalized) estimated by CARD for the newly grided spatial locations.
#' The number of initial grids are defined by the user.
#'
#' @return Return an object of CARDfree class
setClass("CARDfree",
        slots = list(
            spatial_countMat = "ANY",
            spatial_location = "data.frame",
            Proportion_CARD = "matrix",
            estimated_refMatrix = "matrix",
            project = "character",
            markerList = "list",
            info_parameters = "list",
            algorithm_matrix = "list",
            refined_prop = "matrix",
            refined_expression = "matrix"
        ))


#' Show method for the CARD class
#'
#' This method provides a concise summary of an object of class \code{CARD}, 
#' displaying key information including the project name, the number of spots,
#' the number of cell types, and a sample of the 
#' \code{Proportion_CARD} matrix.
#'
#' @param object An object of class \code{CARD}.
#' 
#' @return A concise summary of the \code{CARD} object is printed to the 
#' console.
#' 


setMethod("show", "CARD", function(object) {
    cat("An object of class 'CARD'\n")
    cat("\nProject: ", object@project, "\n")
    cat("Number of spots: ", dim(object@spatial_countMat)[2], "\n")
    cat("Number of cell types: ",length(unique(object@sc_eset$cellType)), "\n")
})



#' Show method for the CARDfree class
#'
#' This method provides a concise summary of an object of 
#' class \code{CARDfree}, displaying key information including the project 
#' name, the number of spots, the number of cell types, and a sample of the 
#' \code{Proportion_CARD} matrix.
#'
#' @param object An object of class \code{CARDfree}.
#' 
#' @return A concise summary of the \code{CARDfree} object is 
#' printed to the console.
#' 

setMethod("show", "CARDfree", function(object) {
    cat("An object of class 'CARDfree'\n")
    cat("\nProject: ", object@project, "\n")
    cat("Number of spots: ", dim(object@spatial_countMat)[2], "\n")
})


