analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(limma)
    go=readRDS(input$rds)
    topgo <- topGO(go, sort="up")

    write.table(go, file = output$tsv, quote = FALSE, sep = '\t', col.names = NA)
}

analysis(snakemake@input, snakemake@output, snakemake@log)