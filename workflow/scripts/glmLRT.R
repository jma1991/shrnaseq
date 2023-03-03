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

    object <- readRDS(input$rds[1])

    contrasts <- readRDS(input$rds[2])

    results <- glmLRT(object, contrast = contrasts[, params$contrast])

    saveRDS(results, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
