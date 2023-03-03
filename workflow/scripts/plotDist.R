# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

plotDist <- function(object, ...) {

    UseMethod("plotDist")

}

plotDist.DGEList <- function(object, group) {

    require(RColorBrewer)

    expr <- edgeR::cpm(object, log = TRUE)

    dist <- dist(t(expr))

    brew <- RColorBrewer::brewer.pal(5, "Reds")

    cols <- colorRampPalette(rev(brew))(255)

    anno <- object$samples[, group, drop = FALSE]

    pheatmap(
        mat = as.matrix(dist),
        color = cols,
        clustering_distance_rows = dist,
        clustering_distance_cols = dist,
        treeheight_row = 0,
        annotation_row = anno,
        show_colnames = FALSE,
        silent = TRUE
    )

}

main <- function(input, output, params, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    library(pheatmap)

    object <- readRDS(input$rds)

    graphic <- plotDist(object, group = params$group)

    png(output$png, width = 8, height = 4, units = "in", res = 300)

    grid::grid.newpage()

    grid::grid.draw(graphic$gtable)

    dev.off()

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
