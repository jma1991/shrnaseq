analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    lrt=readRDS(input$rds)
    write.table(topTags(lrt), output$tsv, row.names=F, quote=F)
}

analysis(snakemake@input, snakemake@output, snakemake@log)