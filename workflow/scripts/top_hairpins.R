analysis=function(input, output) {
    library(edgeR)
    lrt=readRDS(input$rds)
    write.table(topTags(lrt), output$txt, row.names=F, quote=F)
}

analysis(snakemake@input, snakemake@output)