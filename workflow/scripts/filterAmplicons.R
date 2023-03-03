#!/usr/bin/env Rscript

main <- function(input, output, params, log) {

    # Log 

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    object <- readRDS(input$rds)

    keep <- filterByExpr(object, group = object$samples$Condition)

    object <- object[keep, , keep.lib.sizes = FALSE]

    saveRDS(object, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
