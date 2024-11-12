test_that("Visulization for each cell type works well", {
    ## use global reference-based deconvolution result
    expect_no_error({
        p1 <- CARD_visualize_pie(
            proportion = CARD_obj$Proportion_CARD,
            spatial_location = spatialCoords(CARD_obj), 
            colors = NULL, 
            radius = NULL)
    })
})


test_that("Visulization foe selected cell types works well", {
    ## use global reference-based deconvolution result
    ct.visualize = c("Acinar_cells","Cancer_clone_B",
                     "Ductal_terminal_ductal_like",
                     "Ductal_CRISP3_high-centroacinar_like",
                     "Ductal_MHC_Class_II","Ductal_APOL1_high-hypoxic")
    ## visualize the spatial distribution of the cell type proportion
    expect_no_error({
        p2 <- CARD_visualize_prop(
            proportion = CARD_obj$Proportion_CARD,        
            spatial_location = spatialCoords(CARD_obj), 
            ct_visualize = ct.visualize,                 
            colors = c("lightblue","lightyellow","red"), 
            NumCols = 4,                                 
            pointSize = 3.0) 
    })
})

test_that("Visualization for two cell types works well", {
    ## use global reference-based deconvolution result
    ## visualize the spatial distribution of the cell type proportion
    expect_no_error({
        p3 <- CARD_visualize_prop_2CT(
            proportion = CARD_obj$Proportion_CARD,                          
            spatial_location = spatialCoords(CARD_obj),                     
            ct2_visualize = c("Acinar_cells","Cancer_clone_B"),              
            colors = list(c("lightblue","lightyellow","red"),
                          c("lightblue","lightyellow","black")))
    })
})

test_that("Visualization for correlation works well", {
    ## use global reference-based deconvolution result
    expect_no_error({
        p4 <- CARD_visualize_Cor(CARD_obj$Proportion_CARD,colors = NULL)
    })
})


test_that("Visualization for genes works well", {
    ## use global reference-based deconvolution result
    expect_no_error({
        p7 <- CARD_visualize_gene(
            spatial_expression = assays(CARD_obj)$spatial_countMat,
            spatial_location = spatialCoords(CARD_obj),
            gene_visualize = c("A4GNT", "AAMDC", "CD248"),
            colors = NULL,
            NumCols = 3)
    })
})




