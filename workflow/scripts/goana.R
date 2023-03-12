# Author: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

main <- function(input, output, params, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(limma)

    library(
        package = params$organism,
        character.only = TRUE,
        lib.loc = "resources/bioconductor/organism/lib/R/library"
    )

    res <- read.delim(input$tsv)

    sig <- subset(res, FDR < params$FDR)

    org <- getFromNamespace(params$organism, params$organism)

    dir.create(path = output$dir)

    # Up

    sig.up <- subset(sig, direction == "up")

    if (nrow(sig.up) > 0) {

        ids.up <- mapIds(
            x = org,
            keys = rownames(sig.up),
            column = "ENTREZID",
            keytype = "SYMBOL"
        )

        out.up <- goana(
            de = ids.up,
            species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2]
        )

        write.table(
            x = out.up,
            file = file.path(output$dir, "up.tsv"),
            quote = FALSE,
            sep = "\t",
            col.names = NA
        )

    }

    # Down

    sig.down <- subset(sig, direction == "down")

    if (nrow(sig.down) > 0) {

        ids.down <- mapIds(
            x = org,
            keys = rownames(sig.down),
            column = "ENTREZID",
            keytype = "SYMBOL"
        )

        out.down <- goana(
            de = ids.down,
            species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2]
        )

        write.table(
            x = out.down,
            file = file.path(output$dir, "down.tsv"),
            quote = FALSE,
            sep = "\t",
            col.names = NA
        )

    }

    # Mixed

    sig.mixed <- subset(sig, direction == "mixed")

    if (nrow(sig.mixed) > 0) {

        ids.mixed <- mapIds(
            x = org,
            keys = rownames(sig.mixed),
            column = "ENTREZID",
            keytype = "SYMBOL"
        )

        out.mixed <- goana(
            de = ids.mixed,
            species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2]
        )

        write.table(
            x = out.mixed,
            file = file.path(output$dir, "mixed.tsv"),
            quote = FALSE,
            sep = "\t",
            col.names = NA
        )

    }

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
