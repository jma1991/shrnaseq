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
        params$organism,
        character.only = TRUE,
        lib.loc = "resources/bioconductor/organism/lib/R/library"
    )

    res <- read.delim(input$tsv)

    sig <- subset(res, FDR < params$FDR)

    org <- getFromNamespace(params$organism, params$organism)

    ids <- mapIds(
        x = org,
        keys = rownames(sig),
        column = "ENTREZID",
        keytype = "SYMBOL"
    )

    out <- kegga(
        de = ids,
        species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2]
    )

    write.table(
        x = out,
        file = output$tsv,
        quote = FALSE,
        sep = "\t",
        col.names = NA
    )

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
