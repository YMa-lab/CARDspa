% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{createCARDfreeObject}
\alias{createCARDfreeObject}
\title{Create the CARD object}
\usage{
createCARDfreeObject(
  markerList,
  spatial_count,
  spatial_location,
  minCountGene = 100,
  minCountSpot = 5,
  spe = NULL
)
}
\arguments{
\item{markerList}{a list of marker genes, with each element of the list
being the vector of cell type specific marker genes}

\item{spatial_count}{Raw spatial resolved transcriptomics data, each column
is a spatial location, and each row is a gene.}

\item{spatial_location}{data frame, with two columns representing the x and
y coordinates of the spatial location. The rownames of this data frame
should match eaxctly with the columns of the spatial_count.}

\item{minCountGene}{Minimum counts for each gene}

\item{minCountSpot}{Minimum counts for each spatial location}

\item{spe}{a \code{SpatialExperiment} object containing spatial data in
the \code{counts} assay, and spatial coordinates in the spatialCoords.}
}
\value{
Returns CARDfree object with filtered spatial count and marker gene
list.
}
\description{
Create the CARD object
}
\examples{
data(markerList)
data(spatial_count)
data(spatial_location)
CARDfree_obj <- createCARDfreeObject(
    markerList = markerList,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    minCountGene = 100,
    minCountSpot = 5,
    spe = NULL
)

}
