% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{sc_QC}
\alias{sc_QC}
\title{Quality control of scRNA-seq count data}
\usage{
sc_QC(
  counts_in,
  metadata,
  ct_varname,
  ct_select,
  sample_varname = NULL,
  min_cells = 0,
  min_genes = 0
)
}
\arguments{
\item{counts_in}{Raw scRNAseq count data, each column is a cell and each
row is a gene.}

\item{metadata}{data frame, metadata with "ct_varname" specify the cell type
annotation information and "sample_varname" specify the sample information}

\item{ct_varname}{character, the name of the column in metadata that
specifies the cell type annotation information}

\item{ct_select}{vector of cell type names that you are interested in to
deconvolute, default as NULL. If NULL, then use all cell types provided by
single cell dataset;}

\item{sample_varname}{character,the name of the column in metadata that
specifies the sample information. If NULL, we just use the whole as one
sample.}

\item{min_cells}{numeric, we filtered out the non-expressed cells.}

\item{min_genes}{numeric we filtered out the non-expressed genes}
}
\value{
Return the filtered scRNA-seq data and meta data stored in a S4 class
(SingleCellExperiment)
}
\description{
Quality control of scRNA-seq count data
}
