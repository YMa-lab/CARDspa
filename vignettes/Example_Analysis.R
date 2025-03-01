## ----include = FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE--------------------------------------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)) {
#      install.packages("BiocManager")
#  }
#  BiocManager::install("CARDspa")

## ----------------------------------------------------------------------------------------------------------
library(CARDspa)
library(RcppML)
library(NMF)
library(RcppArmadillo)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
#### load the example spatial transcriptomics count data,
data(spatial_count)
spatial_count[1:4, 1:4]

## ----------------------------------------------------------------------------------------------------------
#### load the example spatial location data,
data(spatial_location)
spatial_location[1:4, ]

## ----------------------------------------------------------------------------------------------------------
data(sc_count)
sc_count[1:4, 1:4]

## ----------------------------------------------------------------------------------------------------------
data(sc_meta)
sc_meta[1:4, ]

