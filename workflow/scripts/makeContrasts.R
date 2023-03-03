#!/usr/bin/env Rscript

makeContrasts <- function(object, config) {

    data <- object$samples

    condition <- factor(data$Condition)

    groups <- levels(condition)

    contrasts <- sapply(config$contrasts, paste, collapse = "-")

    contrasts <- limma::makeContrasts(contrasts = contrasts, levels = groups)

}

main <- function(input, output, log, config) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(limma)

    object <- readRDS(input$rds)

    contrasts <- makeContrasts(object, config)

    saveRDS(contrasts, file = output$rds)
}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@config)
