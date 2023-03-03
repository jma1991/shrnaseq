# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

plotMDS <- function(object, ...) {

    UseMethod("plotMDS")

}

plotMDS.DGEList <- function(object, group) {

    mat <- cpm(object, log = TRUE)

    var <- matrixStats::rowVars(mat)

    num <- min(500, length(var))

    ind <- order(var, decreasing = TRUE)[seq_len(num)]

    dst <- dist(t(mat[ind, ]))

    mds <- cmdscale(as.matrix(dst))

    dat <- data.frame(
        MD1 = mds[, 1],
        MD2 = mds[, 2],
        group = object$samples[, group]
    )

    ggplot(dat, aes(MD1, MD2, colour = group)) +
        geom_point(size = 3) +
        labs(x = "MDS 1", y = "MDS 2", colour = "Group")

}

main <- function(input, output, params, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    library(ggplot2)

    object <- readRDS(input$rds)

    graphic <- plotMDS(object, group = params$group)

    ggsave(output$png, plot = graphic, width = 8, height = 6)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
