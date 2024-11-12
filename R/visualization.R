#' Visualize the spatial distribution of cell type proportion
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in
#' either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param ct_visualize Vector of selected cell type names that are
#' interested to visualize
#' @param colors Vector of color names that you want to use, if NULL, we will
#' use the default color scale c("lightblue","lightyellow","red")
#' @param NumCols Numeric, number of columns in the figure panel, it depends on
#' the number of cell types you want to visualize.
#' @param pointSize Size of each point used for plotting
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return Returns a ggplot2 figure.
#'
#' @export
#' @examples
#' library(ggplot2)
#' library(SpatialExperiment)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
#' ct_visualize <- c(
#'     "Acinar_cells", "Cancer_clone_A", "Cancer_clone_B",
#'     "Ductal_terminal_ductal_like", "Ductal_CRISP3_high-centroacinar_like",
#'     "Ductal_MHC_Class_II", "Ductal_APOL1_high-hypoxic", "Fibroblasts"
#' )
#' CARD_visualize_prop(
#'     proportion = CARD_obj$Proportion_CARD,
#'     spatial_location = spatialCoords(CARD_obj),
#'     ct_visualize = ct_visualize,
#'     colors = c("lightblue", "lightyellow", "red"),
#'     NumCols = 4,
#'     pointSize = 3.0
#' )
#'
CARD_visualize_prop <- function(
        proportion, 
        spatial_location,
        ct_visualize = ct_visualize,
        colors = c("lightblue", "lightyellow", "red"),
        NumCols, 
        pointSize = 3.0) {
    #####Preapare dataframe for plotting
    if (is.null(colors)) {
        colors <- c("lightblue", "lightyellow", "red")
    } else {
        colors <- colors
    }
    res_CARD <- as.data.frame(proportion)
    res_CARD <- res_CARD[, order(colnames(res_CARD))]
    location <- as.data.frame(spatial_location)
    if (sum(rownames(res_CARD) == rownames(location)) != nrow(res_CARD)) {
        stop("The rownames of proportion data does not match with the
            rownames of spatial location data")
    }
    ct_select <- ct_visualize
    res_CARD <- res_CARD[, ct_select]
    if (!is.null(ncol(res_CARD))) {
        res_CARD_scale <- as.data.frame(apply(res_CARD, 2, function(x) {
            (x - min(x)) / (max(x) - min(x))}))
    } else {
        res_CARD_scale <-
            as.data.frame((res_CARD - min(res_CARD)) /
                        (max(res_CARD) - min(res_CARD)))
        colnames(res_CARD_scale) <- ct_visualize
    }
    res_CARD_scale$x <- as.numeric(location$x)
    res_CARD_scale$y <- as.numeric(location$y)
    mdata <- melt(res_CARD_scale, id.vars = c("x", "y"))
    colnames(mdata)[3] <- "Cell_Type"
    b <- c(0, 1)
    
    #####Plot
    p <- ggplot(mdata, aes(x, y)) +
        geom_point(aes(colour = value), size = pointSize) +
        scale_color_gradientn(colours = colors) +
        scale_x_discrete(expand = c(0, 1)) +
        scale_y_discrete(expand = c(0, 1)) +
        facet_wrap(~Cell_Type, ncol = NumCols) +
        coord_fixed() +
        theme(
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = "grey89", fill = NA, 
                                        linewidth = 0.5),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 11),
            strip.text = element_text(size = 12, face = "bold"),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.key.size = unit(0.45, "cm")
        )
    return(p)
}

#' Visualize the spatial distribution of two cell type proportions on
#' the same plot
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in
#' either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param ct2_visualize Vector of selected two cell type names that are
#' interested to visualize, here we only focus on two cell types
#' @param colors list of color names that you want to use for each cell type,
#' if NULL, we will use the default color scale list
#' list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return Returns a ggplot2 figure.
#'
#' @export
#' @examples
#' library(ggplot2)
#' library(SpatialExperiment)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
#' CARD_visualize_prop_2CT(
#'     proportion = CARD_obj$Proportion_CARD,
#'     spatial_location = spatialCoords(CARD_obj),
#'     ct2_visualize = c("Cancer_clone_A", "Cancer_clone_B"),
#'     colors = list(c("lightblue", "lightyellow", "red"), c(
#'         "lightblue", "lightyellow",
#'         "black"
#'     ))
#' )
#'
CARD_visualize_prop_2CT <- function(
        proportion, 
        spatial_location,
        ct2_visualize = ct2_visualize,
        colors = NULL) {
    #####Preapare dataframe for plotting
    if (is.null(colors)) {
        colors <- list(
            c("lightblue", "lightyellow", "red"),
            c("lightblue", "lightyellow", "black")
        )
    } else {
        colors <- colors
    }
    res_CARD <- as.data.frame(proportion)
    res_CARD <- res_CARD[, order(colnames(res_CARD))]
    location <- as.data.frame(spatial_location)
    if (sum(rownames(res_CARD) == rownames(location)) != nrow(res_CARD)) {
        stop("The rownames of proportion data does not match with the
            rownames of spatial location data")
    }
    ct_select <- ct2_visualize
    res_CARD <- res_CARD[, ct_select]
    res_CARD_scale <- as.data.frame(apply(res_CARD, 2, function(x) {
        (x - min(x)) / (max(x) - min(x))
    }))
    res_CARD_scale$x <- as.numeric(location$x)
    res_CARD_scale$y <- as.numeric(location$y)
    mdata <- melt(res_CARD_scale, id.vars = c("x", "y"))
    colnames(mdata)[3] <- "Cell_Type"
    b <- c(0, 1)
    
    #####Plot
    p <- ggplot() +
        geom_point(
            data = mdata[mdata$Cell_Type == ct2_visualize[1], ],
            aes(x = x,y = y,color = value), 
            shape = 21, 
            size = 5
        ) +
        scale_color_gradientn(colours = colors[[1]]) +
        geom_point(
            data = mdata[mdata$Cell_Type == ct2_visualize[2], ],
            aes(x = x,y = y,fill = value), 
            color = "white", 
            shape = 22, 
            size = 2
        ) +
        scale_fill_gradientn(colours = colors[[2]]) +
        scale_x_discrete(expand = c(0, 1)) +
        scale_y_discrete(expand = c(0, 1)) +
        coord_fixed() +
        theme(
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(
                colour = "grey89", 
                fill = NA, 
                linewidth = 0.5),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 11),
            strip.text = element_text(size = 12, face = "bold"),
            legend.key.size = unit(0.45, "cm")
        ) +
        labs(colour = ct2_visualize[1], fill = ct2_visualize[2])
    return(p)
}

#' Visualize the spatial distribution of cell type proportion in a geom
#' scatterpie plot
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in
#' either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param colors Vector of color names that you want to use, if NULL, we will
#' use the color palette "Spectral" from RColorBrewer package.
#' @param radius Numeric value about the radius of each pie chart, if NULL, we
#' will calculate it inside the function.
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scatterpie geom_scatterpie
#' @importFrom grDevices colorRampPalette
#' @importFrom gtools mixedsort
#' @return Returns a ggplot2 figure.
#'
#' @export
#' @examples
#' library(ggplot2)
#' library(SpatialExperiment)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
#' colors <- c(
#'     "#FFD92F", "#4DAF4A", "#FCCDE5", "#D9D9D9", "#377EB8", "#7FC97F",
#'     "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", 
#'     "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
#'     "#E6AB02", "#A6761D"
#' )
#' CARD_visualize_pie(
#'     proportion = CARD_obj$Proportion_CARD,
#'     spatial_location = spatialCoords(CARD_obj),
#'     colors = colors,
#'     radius = 0.52
#' )
CARD_visualize_pie <- function(
        proportion, 
        spatial_location, 
        colors = NULL,
        radius = NULL) {
    #####Preapare dataframe for plotting
    res_CARD <- as.data.frame(proportion)
    res_CARD <- res_CARD[, mixedsort(colnames(res_CARD))]
    location <- as.data.frame(spatial_location)
    if (sum(rownames(res_CARD) == rownames(location)) != nrow(res_CARD)) {
        stop("The rownames of proportion data does not match with the rownames
        of spatial location data")
    }
    colorCandidate <- c(
        "#1e77b4", "#ff7d0b", "#ceaaa3", "#2c9f2c", "#babc22",
        "#d52828", "#9267bc",
        "#8b544c", "#e277c1", "#d42728", "#adc6e8", "#97df89", "#fe9795", 
        "#4381bd",
        "#f2941f", "#5aa43a", "#cc4d2e", "#9f83c8", "#91675a",
        "#da8ec8", "#929292", "#c3c237", "#b4e0ea", "#bacceb", "#f7c685",
        "#dcf0d0", "#f4a99f", "#c8bad8",
        "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
        "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
        "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
        "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785",
        "#f4f1de", "#e07a5f", "#3d405b", "#81b29a", "#f2cc8f", "#a8dadc",
        "#f1faee", "#f08080"
    )
    if (is.null(colors)) {
        if (ncol(res_CARD) > length(colorCandidate)) {
            colors <- colorRampPalette(colorCandidate)(ncol(res_CARD))
        } else {
            colors <-
                colorCandidate[sample(
                    seq_len(length(colorCandidate)),
                    ncol(res_CARD)
                )]
        }
    } else {
        colors <- colors
    }
    data <- cbind(res_CARD, location)
    ct_select <- colnames(res_CARD)
    if (is.null(radius)) {
        radius <- (max(data$x) - min(data$x)) * (max(data$y) - min(data$y))
        radius <- radius / nrow(data)
        radius <- radius / pi
        radius <- sqrt(radius) * 0.85
    } else {
        #### avoid the situation when the radius does not generate the correct
        #### figure
        radius <- radius
    }
    
    #####Plot
    p <- 
        ggplot() +
            geom_scatterpie(
                data = data,
                aes(x = x, y = y, r = radius),
                cols = ct_select,
                color = NA
            ) +
            coord_fixed(ratio = 1 * max(data$x) / max(data$y)) +
            scale_fill_manual(values = colors) +
            theme(
                plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                panel.background = element_blank(),
                plot.background = element_blank(),
                panel.border = element_rect(
                    colour = "grey89", 
                    fill = NA, 
                    linewidth = 0.5),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 15),
                legend.key = element_rect(
                    colour = "transparent", 
                    fill = "white"),
                legend.key.size = unit(0.45, "cm"),
                strip.text = element_text(size = 16, face = "bold"),
                legend.position = "bottom"
            ) +
            guides(fill = guide_legend(title = "Cell Type"))
    
    return(p)
}

#' Visualize the cell type proportion correlation
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in
#' either original resolution or enhanced resolution.
#' @param colors Vector of color names that you want to use, if NULL, we will
#' use the default color scale c("#91a28c","white","#8f2c37")
#'
#' @import ggcorrplot
#' @importFrom stats cor
#' @return Returns a ggcorrplot figure.
#'
#' @export
#' @examples
#' library(ggplot2)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
#' CARD_visualize_Cor(CARD_obj$Proportion_CARD, colors = NULL)
CARD_visualize_Cor <- function(proportion, colors = colors) {
    #####Preapare dataframe for plotting
    proportion <- proportion[, order(colnames(proportion))]
    cor_CARD <- cor(as.matrix(proportion))
    if (is.null(colors)) {
        colors <- c("#91a28c", "white", "#8f2c37")
    } else {
        colors <- colors
    }
    
    #####Plot
    p <- ggcorrplot(
        cor_CARD,
        hc.order = FALSE,
        outline.color = "white",
        tl.srt = 60,
        tl.cex = 18,
        lab_size = 7,
        colors = colors
    ) +
        theme(
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(
                colour = "grey89", 
                fill = NA, 
                linewidth = 0.5),
            axis.text = element_text(size = 12),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 16),
            legend.key = element_rect(
                colour = "transparent", 
                fill = "white"),
            legend.key.size = unit(0.45, "cm")
        ) +
        coord_fixed() +
        ggtitle("Correlation") +
        theme(plot.title = element_text(size = 22, face = "bold"))
    return(p)
}

#' Visualize the spatial distribution of cell type proportion
#'
#' @param spatial_expression Data frame, spatial gene expression in either
#' original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param gene_visualize Vector of selected gene names that are interested to
#' visualize
#' @param colors Vector of color names that you want to use, if NULL, we will
#' use the default color scale in virdis palette
#' @param NumCols Numeric, number of columns in the figure panel, it depends on
#' the number of cell types you want to visualize.
#'
#' @import ggplot2
#' @return Returns a ggplot2 figure.
#'
#' @export
#' @examples
#' library(ggplot2)
#' library(SummarizedExperiment)
#' library(SpatialExperiment)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
#' CARD_visualize_gene(
#'     spatial_expression = assays(CARD_obj)$spatial_countMat,
#'     spatial_location = spatialCoords(CARD_obj),
#'     gene_visualize = c("A4GNT", "AAMDC", "CD248"),
#'     colors = NULL,
#'     NumCols = 3
#' )
#'
CARD_visualize_gene <- function(
        spatial_expression, 
        spatial_location,
        gene_visualize, 
        colors = colors, 
        NumCols) {
    #####Preapare dataframe for plotting
    expression <- as.data.frame(as.matrix(spatial_expression))
    expression <- sweep(expression, 2, colSums(expression), "/")
    location <- as.data.frame(spatial_location)
    if (sum(colnames(expression) == rownames(location)) != nrow(location)) {
        stop("The colnames of expression data does not match with the rownames
        of spatial location data")
    }
    gene_select <- gene_visualize
    if (sum(toupper(gene_select) %in% toupper(rownames(spatial_expression))) !=
        length(gene_select)) {
        stop("There existing selected genes that are not in the 
            expression data!")
    }
    plotdata <- NULL
    for (i in seq_len(length(gene_select))) {
        #### load spatial dataset
        gene <- gene_select[i]
        ind <- which(toupper(rownames(expression)) == toupper(gene))
        df <- as.numeric(expression[ind, ])
        names(df) <- colnames(expression)
        df <- (df - min(df)) / (max(df) - min(df))
        d <- data.frame(
            value = df,
            x = as.numeric(location$x),
            y = as.numeric(location$y)
        )
        d$gene <- gene
        plotdata <- rbind(plotdata, d)
    }
    plotdata$gene <- factor(plotdata$gene, levels = gene_select)
    
    ##### Plot
    p <- ggplot(plotdata, aes(x, y)) +
        geom_point(
            aes(color = value),
            size = 1.5,
            shape = 15,
            position = position_dodge2(width = 0.05)
        ) +
        scale_x_discrete(expand = c(0, 1)) +
        scale_y_discrete(expand = c(0, 1)) +
        coord_equal() +
        facet_wrap(~gene, ncol = NumCols) +
        theme(
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.position = "bottom",
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(
                colour = "grey89", 
                fill = NA, 
                linewidth = 0.5),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 14),
            strip.text = element_text(size = 18, face = "bold"),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.key.size = unit(1.0, "cm")
        ) +
        guides(color = guide_colourbar(title = "Expression"))
    if (is.null(colors)) {
        p <- p + scale_color_viridis_c(
            labels = c("0", "0.25", "0.5", "0.75", "1.0"))
    } else {
        p <- p + scale_color_gradientn(colours = colors)
    }
    return(p)
}
