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

    library(csaw)

    res <- read.delim(input$tsv)

    out <- combineTests(
        ids = res$Gene,
        tab = res,
        weights = NULL,
        pval.col = "PValue",
        fc.col = "logFC",
        fc.threshold = params$fdr
    )

    write.table(
        x = out,
        file = output$tsv,
        quote = FALSE,
        row.names = TRUE,
        sep = "\t"
    )

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
