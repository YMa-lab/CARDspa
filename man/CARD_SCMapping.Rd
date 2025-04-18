% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CARD.SCMapping.R
\name{CARD_scmapping}
\alias{CARD_scmapping}
\title{Extension of CARD into performing single cell Mapping from non-single
cell spatial transcriptomics dataset.}
\usage{
CARD_scmapping(CARD_object, shapeSpot = "Square", numcell, ncore = 10)
}
\arguments{
\item{CARD_object}{CARD object create by the CARD_deconvolution function.}

\item{shapeSpot}{a character indicating whether the sampled spatial
coordinates for single cells locating in a Square-like region or a
Circle-like region. The center of this region is the measured spatial
location in the non-single cell resolution spatial transcriptomics data.
The default is 'Square', the other shape is 'Circle'}

\item{numcell}{a numeric value indicating the number of single cells in
each measured location, we suggest 20 for ST technology, 7 for 10x Viisum
and 2 for Slide-seq}

\item{ncore}{a numeric value indicating the number of cores used to
accelerating the procedure}
}
\value{
Returns a SingleCellExperiment SCE object with the mapped expression
at single cell resolution and the spatial location information of each
single cell
}
\description{
Extension of CARD into performing single cell Mapping from non-single
cell spatial transcriptomics dataset.
}
\examples{
library(SingleCellExperiment)
data(spatial_count)
data(spatial_location)
data(sc_count)
data(sc_meta)
CARD_obj <- CARD_deconvolution(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct_varname = "cellType",
    ct_select = unique(sc_meta$cellType),
    sample_varname = "sampleInfo",
    mincountgene = 100,
    mincountspot = 5
)
scMapping <- CARD_scmapping(
CARD_obj, 
shapeSpot = "Square", 
numcell = 20, 
ncore = 2)
print(scMapping)

}
