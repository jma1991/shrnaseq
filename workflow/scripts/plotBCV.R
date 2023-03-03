# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

main <- function(input, output, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    object <- readRDS(input$rds)

    png(output$png, width = 8, height = 4, units = "in", res = 300)

    plotBCV(object)

    dev.off()

}

main(snakemake@input, snakemake@output, snakemake@log)
