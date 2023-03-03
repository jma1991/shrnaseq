main <- function(input, output, log) {

    # Log

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script

    library(edgeR)

    object <- readRDS(input$rds)

    results <- topTags(object, n = Inf, sort.by = "none")

    write.table(results, file = output$tsv, quote = FALSE, sep = "\t", row.names = FALSE)

}

main(snakemake@input, snakemake@output, snakemake@log)
