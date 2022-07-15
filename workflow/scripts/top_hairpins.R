analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    lrt=readRDS(input$rds[1])
    write.table(topTags(lrt), output$tsv[1], row.names=F, quote=F)
    corrected_lrt=readRDS(input$rds[2])
    write.table(topTags(corrected_lrt), output$tsv[2], row.names=F, quote=F)
}

analysis(snakemake@input, snakemake@output, snakemake@log)