# setup-card.R
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
    mincountgene = 10,
    mincountspot = 5
)
assign("CARD_obj", CARD_obj, envir = .GlobalEnv)
