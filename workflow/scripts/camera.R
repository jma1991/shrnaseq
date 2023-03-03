# Authors: James Ashmore, Claire Prince
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

    library(edgeR)

    library(limma)

    object <- readRDS(input$rds[1])

    ids <- split(object$genes$ID, object$genes$Gene)

    indices <- ids2indices(ids, identifiers = object$genes$ID)

    design <- readRDS(input$rds[2])

    contrasts <- readRDS(input$rds[3])

    results <- camera(
        y = object,
        index = indices,
        design = design,
        contrast = contrasts[, params$contrast],
        sort = FALSE
    )

    write.table(results, file = output$tsv, quote = FALSE, sep = "\t")

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
