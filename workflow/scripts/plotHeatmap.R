# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: james.ashmore@zifornd.com
# License: MIT

pheatmap.mat <- function(x) {

    # Scale rows by 'variance-aware' Z-transformation

    M <- rowMeans(x, na.rm = TRUE)

    DF <- ncol(x) - 1

    isNA <- is.na(x)

    anyNA <- any(isNA)

    if (anyNA) {

        mode(isNA) <- "integer"

        DF <-  DF - rowSums(isNA)

        DF[DF == 0] <- 1

    }

    x <- x - M

    V <- rowSums(x ^ 2, na.rm = TRUE) / DF

    x <- x / sqrt(V + 0.01)

}

pheatmap.color <- function(x) {

    # Return color vector

    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = x)))(100)

}

pheatmap.breaks <- function(x) {

    # Return breaks vector

    abs <- max(abs(x))

    abs <- min(abs, 5)

    seq(-abs, +abs, length.out = 101)

}

pheatmap.cluster_rows <- function(x) {

    # Return hclust object for rows

    hclust(dist(x, method = "euclidean"), method = "complete")

}

pheatmap.cluster_cols <- function(x) {

    # Return hclust object for columns

    hclust(dist(t(x), method = "euclidean"), method = "complete")

}

pheatmap.annotation_col <- function(x) {

    # Return column annotations

    data.frame(
        Condition = x$condition,
        row.names = colnames(x)
    )

}

main <- function(input, output, params, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(edgeR)

    library(pheatmap)

    obj <- readRDS(input$rds)

    res <- read.delim(input$tsv)

    ##

    ind <- order(res$P.Value)[seq_len(params$ntop)]

    ids <- as.character(res$PROBEID[ind])

    obj <- obj[ids, ]

    ##

    cpm <- exprs(obj)

    mat <- pheatmap.mat(cpm)

    lab <- fData(obj)$SYMBOL

    pheatmap(
        mat = mat,
        color = pheatmap.color("RdBu"),
        breaks = pheatmap.breaks(mat),
        cellwidth = 10,
        cellheight = 10,
        cluster_rows = pheatmap.cluster_rows(mat), 
        cluster_cols = pheatmap.cluster_cols(cpm),
        annotation_col = pheatmap.annotation_col(obj),
        show_colnames = FALSE,
        labels_row = lab,
        filename = output$pdf
    )

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
