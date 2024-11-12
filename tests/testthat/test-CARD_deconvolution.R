test_that("Reference-Based Deconvolution works well", {
    expect_true(exists("CARD_obj"))
    expect_true("Proportion_CARD" %in% names(colData(CARD_obj)))
})


test_that("Reference-free Deconvolution works well", {
    library(NMF)
    data(spatial_count)
    data(spatial_location)
    data(markerList)
    markerList2 <- markerList[6:10]
    CARDfree_obj = CARD_refFree(
        markerlist = markerList2,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        mincountgene = 10,
        mincountspot =5
    ) 
    res <- rowSums(CARDfree_obj$Proportion_CARD)
    expect_equal(as.numeric(res), as.numeric(rep(1, 428)))
})

test_that("Refined spatial map works well", {
    ## use global reference-based deconvolution result
    expect_no_error({
        CARD_obj <- CARD_imputation(CARD_obj, num_grids = 2000, ineibor = 10, 
                                    exclude = NULL)
    })
})

test_that("Single cell resolution mapping works well", {
    ## use global reference-based deconvolution result
    expect_no_error({
        scMapping <- CARD_scmapping(CARD_obj,shapeSpot="Square",
                                    numcell=20,ncore=1)
    })
})


