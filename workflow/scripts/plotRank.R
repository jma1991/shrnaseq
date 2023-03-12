# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

plotRank <- function(x, FDR = 0.05) {

    #

    x$rank.logFC <- rank(x$rep.logFC)

    #

    index.min <- sort(x$rep.logFC, index.return = TRUE, decreasing = FALSE)

    index.min <- head(index.min$ix, n = 10)

    index.max <- sort(x$rep.logFC, index.return = TRUE, decreasing = TRUE)

    index.max <- head(index.max$ix, n = 10)

    #

    x$label.logFC <- ""

    x$label.logFC[index.min] <- rownames(x)[index.min]

    x$label.logFC[index.max] <- rownames(x)[index.max]

    #

    x$color.logFC <- x$direction

    x$color.logFC[x$FDR > params$FDR] <- "ns"

    #

    pal <- c(
        "up"    = "#4daf4a",
        "down"  = "#e41a1c",
        "mixed" = "#ff7f00",
        "ns"    = "#999999"
    )

    #

    ggplot(x, aes(rep.logFC, rank.logFC, label = label.logFC)) +
        geom_point(aes(colour = color.logFC)) +
        geom_label_repel(aes(fill = color.logFC), colour = "#ffffff", segment.colour = "#000000", max.overlaps = Inf) +
        scale_colour_manual(values = pal) +
        scale_fill_manual(values = pal) +
        labs(x = "Score", y = "Rank") +
        theme_classic() +
        theme(aspect.ratio = 1, legend.position = "none")

}

main <- function(input, output, params, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(ggplot2)

    library(ggrepel)

    dat <- read.delim(input$tsv, row.names = 1)

    plt <- plotRank(dat, FDR = params$FDR)

    ggsave(output$png, plot = plt)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
