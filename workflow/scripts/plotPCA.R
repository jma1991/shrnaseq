# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

plotPCA <- function(object, ...) {

    UseMethod("plotPCA")

}

plotPCA.DGEList <- function(object, group) {

    mat <- cpm(object, log = TRUE)

    var <- matrixStats::rowVars(mat)

    num <- min(500, length(var))

    ind <- order(var, decreasing = TRUE)[seq_len(num)]

    pca <- prcomp(t(mat[ind, ]))

    pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)

    dat <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = object$samples[, group]
    )

    ggplot(dat, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", round(pct[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(pct[2] * 100), "% variance")) +
        coord_fixed() +
        labs(colour = "Group")

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

    graphic <- plotPCA(object, group = params$group)

    ggsave(output$png, plot = graphic, width = 8, height = 6)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
