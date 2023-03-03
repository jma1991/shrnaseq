#!/usr/bin/env Rscript

main <- function(input, output, params, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    library(limma)

    object <- readRDS(input$rds)

    logcounts <- cpm(object, log = TRUE)

    batch <- factor(object$samples$Batch)

    batch2 <- factor(object$samples$Batch2)

    design <- model.matrix(~ 0 + object$samples$Condition)

    if (nlevels(batch) > 1 | nlevels(batch2) > 1) {

        corrected <- removeBatchEffect(
            x = logcounts,
            batch = batch,
            batch2 = batch2,
            design = design
        )

    } else {

        corrected <- logcounts

    }

    saveRDS(corrected, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
