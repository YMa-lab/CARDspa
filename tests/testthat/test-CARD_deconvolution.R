test_that("Deconvolution works well", {
    data(spatial_count)
    data(spatial_location)
    data(sc_count)
    data(sc_meta)
    CARD_obj <- createCARDObject(
        sc_count = sc_count,
        sc_meta = sc_meta,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        ct.varname = "cellType",
        ct.select = unique(sc_meta$cellType),
        sample.varname = "sampleInfo",
        minCountGene = 100,
        minCountSpot = 5
    )
    CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
    res <- rowSums(slot(CARD_obj, "Proportion_CARD"))
    expect_equal(as.numeric(res), as.numeric(rep(1, 423)))
})
