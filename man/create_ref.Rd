% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CARD.prop.R
\name{create_ref}
\alias{create_ref}
\title{Construct the mean gene expression basis matrix (B), this is the
faster version}
\usage{
create_ref(sc_eset, ct_select = NULL, ct_varname, sample_varname = NULL)
}
\arguments{
\item{sc_eset}{S4 class for storing data from single-cell experiments.
This format is usually created by the package SingleCellExperiment with
stored counts, along with the usual metadata for genes and cells.}

\item{ct_select}{vector of cell type names that you are interested in to
deconvolute, default as NULL. If NULL, then use all cell types provided by
single cell dataset;}

\item{ct_varname}{character, the name of the column in metaData that
specifies the cell type annotation information}

\item{sample_varname}{character,the name of the column in metaData that
specifies the sample information. If NULL, we just use the whole as one
sample.}
}
\value{
Return a list of basis (B) matrix
}
\description{
Construct the mean gene expression basis matrix (B), this is the
faster version
}
