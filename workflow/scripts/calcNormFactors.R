#!/usr/bin/env Rscript

main <- function(input, output, log) {

    # Log 

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    object <- readRDS(input$rds)

    object <- calcNormFactors(object)

    saveRDS(object, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
