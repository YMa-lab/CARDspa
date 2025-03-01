# CARDspa

## Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics 
We developed CARDspa, a Bioconductor version of the original CARD method 
[https://github.com/YMa-lab/CARD], for spatially informed cell type 
deconvolution in spatial transcriptomics. CARDspa is a reference-based 
deconvolution tool that estimates the cell type composition in spatial 
transcriptomics data by utilizing cell-type-specific expression information 
obtained from reference single-cell RNA sequencing (scRNA-seq) data.

A key feature of CARDspa is its ability to accommodate spatial correlation in 
the cell type composition across different tissue locations, which ensures 
accurate and spatially informed deconvolution, as well as refined spatial map 
construction. This is achieved through an efficient optimization algorithm for 
constrained maximum likelihood estimation. CARDspa is scalable to handle spatial 
transcriptomics data with tens of thousands of spatial locations and genes.

As a Bioconductor package, CARDspa is designed to integrate seamlessly with 
other Bioconductor tools for analyzing and visualizing high-dimensional 
biological data. It leverages Bioconductor's robust framework for 
reproducibility and ease of use within the R ecosystem.

CARDspa is implemented as an open-source R package, and is freely available on 
Bioconductor at [https://bioconductor.org/packages/CARDspa/].

For further details, please visit our project page at 
[www.xzlab.org/software.html]. 

Installation
------------
You can install the released version of CARDspa using BioConductor.

``` r
# install BiocManager if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CARDspa")

```

Or you can install using devtools in R.
``` r
# install devtools if necessary
install.packages('devtools')

# install the CARDspa package
devtools::install_github('https://github.com/YMa-lab/CARDspa')

# load package
library(CARDspa)

```

# Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure 
to raise issues with a detailed and reproducible exmple and also please 
provide the output of your sessionInfo() in R! 

How to cite `CARDspa`
-------------------
Ma, Y., Zhou, X. Spatially informed cell-type deconvolution for spatial 
transcriptomics. Nat Biotechnol 40, 1349–1359 (2022). 
https://doi.org/10.1038/s41587-022-01273-7

How to use `CARDspa`
-------------------
Please see vignettes.
``` r
browseVignettes("CARDspa")
```

